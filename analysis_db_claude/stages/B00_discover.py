# %%
"""
B00: discover the tracks of a batch.

    python stages/B00_discover.py <batch_key>

Reads  batches/<batch_key>.toml
Writes work/<batch>/tracks.csv               one row per expected track (ID from the granule name)
       work/<batch>/chunks.csv               time chunks for the B01 download jobs
       work/<batch>/rgt_start_points.geojson lowest-latitude point of every RGT crossing the box
                                             (the only step that touches the 510 MB RGT shapefile)
       work/<batch>/batch.json               copy of the batch definition
       plots/<hemis>/<batch>/_batch/B00_overview.png  polar map + estimates
       status/B00/<batch_key>.json
"""
import re
import sys
import json
import datetime as dt
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # analysis_db_claude/

from pipeline_config import HERE, mconfig, np, pd, plt, MT, paths_for, save_fig
from pipeline_status import StageRun, SkipTrack
from pipeline_params import load_batch, load_params, materialize, region_polygon
import ICEsat2_SI_tools.sliderule_converter_tools as sct

STAGE = 'B00'

# ATL03_YYYYMMDDHHMMSS_RRRRCCSS_RL_VV.h5
GRANULE_RX = re.compile(r'ATL03_(\d{4})(\d{2})(\d{2})(\d{2})(\d{2})(\d{2})_(\d{4})(\d{2})(\d{2})_(\d{3})_(\d{2})\.h5')

# rough data-volume estimates, from the test track SH_20190502_05180312 (2026-09-21)
BYTES_PER_POINT_B01 = 120          # B01_binned.h5, all columns, one beam point
MB_INTERMEDIATE_PER_TRACK = 150    # B02..B06 products per track (gFT_x/k.nc dominate)


def parse_granule(name):
    m = GRANULE_RX.search(name)
    if not m:
        return None
    y, mo, d, H, M, S, rgt, cyc, seg, rl, vv = m.groups()
    return dict(granule=name, date=f'{y}{mo}{d}', time=dt.datetime(int(y), int(mo), int(d), int(H), int(M), int(S)),
                rgt=int(rgt), cycle=int(cyc), segment=int(seg), release=rl, version=vv)


def make_chunks(t0, t1, chunk_days):
    """[(chunk_index, t0, t1), ...] consecutive windows covering [t0, t1)"""
    if not chunk_days or chunk_days <= 0:
        return [(0, t0, t1)]
    out, i, a = [], 0, t0
    while a < t1:
        b = min(a + dt.timedelta(days=chunk_days), t1)
        out.append((i, a, b))
        i, a = i + 1, b
    return out


def cmr_granules(part, t0, t1, depth=0):
    """ATL03 granules in one polygon part and time window; SlideRule's CMR proxy refuses more than
    300 hits per query, so the window is halved recursively until every query fits"""
    from sliderule import earthdata
    fmt = '%Y-%m-%dT%H:%M:%S'
    try:
        return list(earthdata.cmr(short_name='ATL03', polygon=part, time_start=t0.strftime(fmt), time_end=t1.strftime(fmt)))
    except Exception as e:
        if 'exceeded maximum' in str(e) and depth < 12 and (t1 - t0) > dt.timedelta(hours=6):
            tm = t0 + (t1 - t0) / 2
            print(f'  CMR: > 300 hits in {t0:%Y-%m-%d}..{t1:%Y-%m-%d}, splitting')
            return cmr_granules(part, t0, tm, depth + 1) + cmr_granules(part, tm, t1, depth + 1)
        raise


def polar_overview(batch, poly, Gtrack, Gtrack_lowest, tracks, chunks, est, path):
    """polar map (r = 90 - |lat|) with the box, RGT points in the box and the estimates"""
    hemis = batch['batch']['hemis']
    sign = -1 if hemis == 'SH' else 1
    r_of = lambda lat: 90 - np.abs(np.asarray(lat))
    th_of = lambda lon: np.deg2rad(np.asarray(lon))

    fig = plt.figure(figsize=(7.5, 8))
    ax = fig.add_subplot(111, projection='polar')
    ax.set_theta_zero_location('S' if hemis == 'SH' else 'N')
    ax.set_theta_direction(-1 if hemis == 'SH' else 1)
    r_max = 40 if hemis == 'SH' else 40
    ax.set_rlim(0, r_max)
    lat_ticks = np.arange(sign * 90, sign * (90 - r_max) + (0 if sign < 0 else 1), sign * -10)
    ax.set_rgrids(r_of(lat_ticks), labels=[f'{int(l)}°' for l in lat_ticks], angle=0, fontsize=7)
    ax.grid(alpha=0.4)

    # coastline (Natural Earth 50 m, analysis_db_claude/support/), clipped to the hemisphere
    try:
        import geopandas as gpd
        coast = gpd.read_file(str(HERE / 'support' / 'ne_50m_coastline.zip'))
        for geom in coast.geometry:
            lines = geom.geoms if geom.geom_type == 'MultiLineString' else [geom]
            for ln in lines:
                lon, lat = np.asarray(ln.coords).T
                keep = (sign * lat) > (90 - r_max)
                if keep.any():
                    lat = np.where(keep, lat, np.nan)          # break lines at the map edge
                    ax.plot(th_of(lon), r_of(lat), '-', color='0.35', linewidth=0.5, zorder=1)
    except Exception as e:
        print('coastline not drawn:', repr(e))

    # RGT ground-track points inside the box (thinned) and their start points
    if len(Gtrack):
        g = Gtrack.iloc[::20]
        ax.plot(th_of(g.geometry.x), r_of(g.geometry.y), '.', color='0.6', markersize=1, label='RGT points in box')
    if len(Gtrack_lowest):
        ax.plot(th_of(Gtrack_lowest.geometry.x), r_of(Gtrack_lowest.geometry.y), 'o', color='tab:blue',
                markersize=3, label=f'RGT start points ({len(Gtrack_lowest)})')

    # the region, densified along the edges (the unwrapped ring; theta wraps naturally on the polar axes)
    lons = [p['lon'] for p in poly['list']]
    lats = [p['lat'] for p in poly['list']]
    bl, bb = [], []
    for i in range(len(lons) - 1):
        bl += list(np.linspace(lons[i], lons[i + 1], 50))
        bb += list(np.linspace(lats[i], lats[i + 1], 50))
    ax.plot(th_of(bl), r_of(bb), '-', color='tab:green', linewidth=2,
            label='batch polygon' + (' (crosses 180°, %d parts)' % len(poly['parts']) if poly['crosses_antimeridian'] else ''))

    ax.set_title(f"{batch['batch']['key']}  ({hemis}, polar view, r = 90-|lat|)\n"
                 f"{poly.get('kind', 'box')}: lat {poly['lats'][0]:.2f}..{poly['lats'][1]:.2f}  lon {poly['lons'][0]:.2f}..{poly['lons'][1]:.2f}",
                 fontsize=9, loc='left')
    ax.legend(loc='lower right', fontsize=7, bbox_to_anchor=(1.15, -0.12))

    sel = tracks[tracks.selected]
    txt = (f"time      : {batch['time']['t0']} .. {batch['time']['t1']}  ({len(chunks)} chunk(s))\n"
           f"granules  : {len(tracks)} in CMR, {len(sel)} selected\n"
           f"tracks    : {sel.ID.nunique()} expected (unique rgt,cycle), {sel.rgt.nunique()} RGTs\n"
           f"RGTs in box (shapefile): {Gtrack_lowest.RGT.nunique() if len(Gtrack_lowest) else 0}\n"
           f"est. B01 download : {est['B01_MB']:.0f} MB  ({est['points_per_track']/1e3:.0f}k points/track)\n"
           f"est. B02-B06 data : {est['intermediate_MB']/1e3:.1f} GB\n"
           f"est. B02+B04 time : {est['cpu_hours']:.1f} CPU h (serial ~6 min/track)")
    fig.text(0.02, 0.01, txt, fontsize=8, family='monospace', va='bottom',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
    save_fig(fig, path, 'B00_overview', pdf=False)


def run_stage(batch_key, run):
    import geopandas as gpd
    batch = load_batch(batch_key)
    P = paths_for(batch_key)
    MT.mkdirs_r(P.batch_work)
    hemis = batch['batch']['hemis']
    t0 = dt.datetime.fromisoformat(str(batch['time']['t0']))
    t1 = dt.datetime.fromisoformat(str(batch['time']['t1']))

    # %% polygon and CMR granule list
    poly = region_polygon(batch)
    print('polygon:', poly['list'], '| parts:', len(poly['parts']), '| crosses antimeridian:', poly['crosses_antimeridian'])
    granules = []
    for part in poly['parts']:                       # one CMR query per part (cut at +-180) ...
        for _, ca, cb in make_chunks(t0, t1, batch['time'].get('chunk_days', 0)):   # ... and per time chunk
            granules += cmr_granules(part, ca, cb)
    granules = sorted(set(granules))
    print(f'CMR: {len(granules)} ATL03 granules in box and time window')
    if len(granules) == 0:
        raise SkipTrack('no ATL03 granules in box/time window')

    rows = [parse_granule(g) for g in granules]
    bad = [g for g, r in zip(granules, rows) if r is None]
    if bad:
        print('unparsable granule names (ignored):', bad)
    tracks = pd.DataFrame([r for r in rows if r]).sort_values('time').reset_index(drop=True)

    # %% RGT start points from the shapefile subset in the box
    shp = P.groundtracks + f'IS2_mission_points_{hemis}_RGT_all.shp'
    print('reading RGT shapefile subset:', shp)
    Gtrack = gpd.read_file(shp, mask=poly['shapely'])
    Gtrack_lowest = sct.get_RGT_start_points(Gtrack) if len(Gtrack) else gpd.GeoDataFrame(columns=['RGT', 'geometry'])
    if len(Gtrack_lowest):
        Gtrack_lowest['RGT'] = Gtrack_lowest['RGT'].astype(int)
        Gtrack_lowest.to_file(P.batch_work + 'rgt_start_points.geojson', driver='GeoJSON')
    print(f'RGTs in box: {Gtrack_lowest.RGT.nunique() if len(Gtrack_lowest) else 0}')

    # %% expected track IDs: one per (rgt, cycle); several granules (segments) may fall in the box
    tracks['ID'] = [f"{hemis}_{r.date}_{r.rgt:04d}{r.cycle:02d}{r.segment:02d}" for r in tracks.itertuples()]
    grp = tracks.groupby(['rgt', 'cycle'], sort=False)
    first = grp.head(1).set_index(['rgt', 'cycle'])['ID']
    tracks['ID'] = [first.loc[(r.rgt, r.cycle)] for r in tracks.itertuples()]
    tracks['n_granules'] = grp['granule'].transform('count')
    tracks['rgt_in_shapefile'] = tracks.rgt.isin(set(Gtrack_lowest.RGT)) if len(Gtrack_lowest) else False
    lowest = Gtrack_lowest.set_index('RGT').geometry if len(Gtrack_lowest) else None
    tracks['rgt_start_lon'] = [lowest.loc[r].x if lowest is not None and r in lowest.index else np.nan for r in tracks.rgt]
    tracks['rgt_start_lat'] = [lowest.loc[r].y if lowest is not None and r in lowest.index else np.nan for r in tracks.rgt]

    # %% chunks
    chunks = make_chunks(t0, t1, batch['time'].get('chunk_days', 0))
    tracks['chunk'] = [next(i for i, a, b in chunks if a <= t < b or (t >= b and i == chunks[-1][0])) for t in tracks.time]
    chunks_df = pd.DataFrame([{'chunk': i, 't0': a.isoformat(), 't1': b.isoformat(),
                               'n_granules': int((tracks.chunk == i).sum())} for i, a, b in chunks])

    # %% selection
    sel = batch['selection']
    selected = pd.Series(True, index=tracks.index)
    if sel.get('require_rgt', True):
        selected &= tracks.rgt_in_shapefile
    if sel.get('include_ids'):
        selected &= tracks.ID.isin(sel['include_ids'])
    if sel.get('exclude_ids'):
        selected &= ~tracks.ID.isin(sel['exclude_ids'])
    if sel.get('max_tracks', 0):
        keep_ids = tracks[selected].drop_duplicates('ID').ID.head(int(sel['max_tracks']))
        selected &= tracks.ID.isin(keep_ids)
    tracks['selected'] = selected
    tracks['time'] = tracks.time.dt.strftime('%Y-%m-%dT%H:%M:%S')
    cols = ['ID', 'granule', 'rgt', 'cycle', 'segment', 'date', 'time', 'release', 'version', 'n_granules',
            'chunk', 'rgt_in_shapefile', 'rgt_start_lon', 'rgt_start_lat', 'selected']
    tracks = tracks[cols]
    tracks.to_csv(P.batch_work + 'tracks.csv', index=False)
    chunks_df.to_csv(P.batch_work + 'chunks.csv', index=False)
    with open(P.batch_work + 'batch.json', 'w') as f:
        json.dump(batch, f, indent=2, default=str)
    materialize(batch_key)

    # %% estimates + overview figure
    res = load_params(batch_key)['B01']['sliderule']['res']
    lat_extent_m = sct.haversine(0, poly['lats'][0], 0, poly['lats'][1]) * 1e3
    n_tracks = int(tracks[tracks.selected].ID.nunique())
    points_per_track = 6 * lat_extent_m / res
    est = {'n_granules': int(len(tracks)), 'n_selected_granules': int(tracks.selected.sum()),
           'n_tracks': n_tracks, 'n_rgts_in_box': int(Gtrack_lowest.RGT.nunique()) if len(Gtrack_lowest) else 0,
           'points_per_track': points_per_track,
           'B01_MB': n_tracks * points_per_track * BYTES_PER_POINT_B01 / 1e6,
           'intermediate_MB': n_tracks * MB_INTERMEDIATE_PER_TRACK,
           'cpu_hours': n_tracks * 6 / 60, 'n_chunks': len(chunks)}
    polar_overview(batch, poly, Gtrack, Gtrack_lowest, tracks, chunks, est, P.plot_batch + '_batch/')
    run.info(**est)
    print(tracks[tracks.selected][['ID', 'granule', 'chunk', 'n_granules']].drop_duplicates('ID').to_string())
    print(f"selected {n_tracks} tracks; est. B01 {est['B01_MB']:.0f} MB, intermediates {est['intermediate_MB']/1e3:.1f} GB")


def main(batch_key):
    prm = load_params(batch_key)
    with StageRun(STAGE, batch_key, batch_key, params={'version': prm['version']}, script=__file__) as run:
        run_stage(batch_key, run)
    return run.status


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('-')]
    batch_key = args[0] if args else 'SH_dev_small'
    sys.exit(0 if main(batch_key) != 'fail' else 1)
