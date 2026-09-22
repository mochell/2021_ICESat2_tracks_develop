# %%
"""
B01: download ATL06-type segments with SlideRule for one time chunk of a batch and split them
into tracks with an absolute along-track coordinate.

    python stages/B01_sliderule_load.py <batch_key> <chunk>

Reads  work/<batch>/tracks.csv, chunks.csv, rgt_start_points.geojson   (from B00)
Writes work/<batch>/B01_regrid/<ID>_B01_binned.h5       one h5py group per beam
       work/<batch>/A01b_ID/A01b_ID_<ID>.json           track metadata used by A02
       plots/<hemis>/<batch>/<ID>/B01b_ATL06_corrected.png, B01b_beam_statistics.png, B01_track.png
       work/<batch>/A01b_ID/<batch>_chunk<n>_point_density.html
       status/B01/<ID>.json for every expected ID of the chunk (success, or skip 'no data in box')
       status/B01/chunk<n>.json for the download job itself
"""
import sys
import copy
import json
import datetime
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline_config import mconfig, np, pd, plt, M, MT, col, paths_for, save_fig, font_for_pres
from pipeline_status import StageRun, SkipTrack, write_status, read_status
from pipeline_params import load_batch, load_params
import ICEsat2_SI_tools.sliderule_converter_tools as sct
import ICEsat2_SI_tools.io as io
import ICEsat2_SI_tools.beam_stats as beam_stats

STAGE = 'B01'
BEAMS = ['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']


def make_B01_dict(table_data, split_by_beam=True, to_hdf5=False):
    """
    converts a GeoDataFrame from SlideRule to one DataFrame per beam with the column names the
    downstream stages expect (from analysis_db/B01_SL_load_single_file.py)
    """
    table_data.rename(columns={'n_fit_photons': 'N_photos',
                               'w_surface_window_final': 'signal_confidence',
                               'y_atc': 'y',
                               'x_atc': 'distance'}, inplace=True)
    table_data['lons'] = table_data['geometry'].x
    table_data['lats'] = table_data['geometry'].y
    drop_columns = [c for c in ['cycle', 'gt', 'rgt', 'pflags'] if c in table_data.columns]
    if to_hdf5:
        drop_columns.append('geometry')
    table_data.drop(columns=drop_columns, inplace=True)
    if not split_by_beam:
        return table_data
    B01b = dict()
    for spot, beam in zip([1, 2, 3, 4, 5, 6], BEAMS):
        B01b[beam] = table_data[table_data.spot == spot]
    return B01b


def track_products(ID, gdf_track, granules, P, prm, Gtrack_lowest):
    """
    absolute x coordinate, per-beam tables, h5 + json + figures for one (rgt, cycle) track.
    returns info dict
    """
    plot_path = P.plot_batch + ID + '/'
    MT.mkdirs_r(plot_path)
    save_path = P.stage_dir('B01_regrid')
    save_path_json = P.stage_dir('A01b_ID')

    table_data = copy.copy(gdf_track)
    ascending = sct.ascending_test_distance(table_data)
    # absolute x: distance from the equator via the RGT start point in the box (sct.define_x_coordinate_with_RGT)
    table_data = sct.define_x_coordinate_with_RGT(table_data, Gtrack_lowest)
    table_data.sort_values(by='x', inplace=True)
    table_data.reset_index(inplace=True)           # 'time' index -> column
    x_ref = float(abs(table_data['x_atc'].iloc[0] - table_data['x'].iloc[0]))   # reference distance from the equator
    x_extent = float(table_data['x'].max() - table_data['x'].min())
    if x_extent < 2 * prm['beam_stats_Lmeter']:
        raise SkipTrack(f'track too short in box ({x_extent/1e3:.1f} km)', x_extent_km=x_extent / 1e3)
    table_time = table_data['time']
    table_data.drop(columns=['time'], inplace=True)

    Ti = make_B01_dict(table_data, split_by_beam=True, to_hdf5=True)
    n_beam = {b: int(Ti[b].shape[0]) for b in BEAMS}
    thin = [b for b, n in n_beam.items() if n < prm['min_points_per_beam']]
    if thin:
        raise SkipTrack(f'beams with < {prm["min_points_per_beam"]} points: {thin}', n_points=n_beam)
    for kk in Ti.keys():
        Ti[kk]['dist'] = Ti[kk]['x'].copy()
        Ti[kk]['heights_c_weighted_mean'] = Ti[kk]['h_mean'].copy()
        Ti[kk]['heights_c_std'] = Ti[kk]['h_sigma'].copy()
    io.write_track_to_HDF5(Ti, ID + '_B01_binned', save_path)

    # figures
    cdict = {s: col.rels[b] for s, b in zip([1, 2, 3, 4, 5, 6], BEAMS)}
    font_for_pres()
    F_atl06 = M.figure_axis_xy(6.5, 5, view_scale=0.6)
    F_atl06.fig.suptitle(ID)
    beam_stats.plot_ATL06_track_data(gdf_track, cdict)
    save_fig(F_atl06, plot_path, 'B01b_ATL06_corrected', pdf=False)

    beam_stats_ok = True
    try:
        D = beam_stats.derive_beam_statistics(Ti, BEAMS, Lmeter=prm['beam_stats_Lmeter'], dx=prm['beam_stats_dx'])
        F = M.figure_axis_xy(8, 4.3, view_scale=0.6)
        beam_stats.plot_beam_statistics(D, mconfig['beams']['high_beams'], mconfig['beams']['low_beams'], col.rels,
                                        track_name=ID + ' |  ascending =' + str(ascending))
        save_fig(F, plot_path, 'B01b_beam_statistics', pdf=False)
    except Exception as e:                  # diagnostics only; a short/sparse beam must not fail the track
        print(f'  beam statistics figure failed for {ID}: {e!r}')
        beam_stats_ok = False
        plt.close('all')

    gdf_track[::100].plot(markersize=0.1, figsize=(4, 6))
    plt.title(ID + '\nascending =' + str(ascending), loc='left')
    save_fig(plt.gcf(), plot_path, 'B01_track', pdf=False)

    # A01b json (read by A02 via io.init_data)
    start_pos = abs(table_data.lats).argmin()
    end_pos = abs(table_data.lats).argmax()
    DD = {'case_ID': ID, 'tracks': {'ATL03': granules},
          'pars': {'poleward': bool(sct.ascending_test(gdf_track)), 'ascending': bool(ascending), 'region': '0',
                   'x_reference_m': x_ref,
                   'start': {'longitude': float(table_data.lons[start_pos]), 'latitude': float(table_data.lats[start_pos]),
                             'seg_dist_x': float(table_data.x[start_pos]),
                             'delta_time': datetime.datetime.timestamp(table_time[start_pos])},
                   'end': {'longitude': float(table_data.lons[end_pos]), 'latitude': float(table_data.lats[end_pos]),
                           'seg_dist_x': float(table_data.x[end_pos]),
                           'delta_time': datetime.datetime.timestamp(table_time[end_pos])}}}
    MT.json_save2(name='A01b_ID_' + ID, path=save_path_json, data=DD)

    n_points = {b: int(Ti[b].shape[0]) for b in BEAMS}
    return {'n_points': n_points, 'n_points_total': int(sum(n_points.values())), 'ascending': bool(ascending),
            'beam_stats_ok': beam_stats_ok,
            'x_reference_m': x_ref, 'x_min_km': float(table_data.x.min() / 1e3), 'x_max_km': float(table_data.x.max() / 1e3),
            'N_photos_median': {b: float(Ti[b].N_photos.median()) if len(Ti[b]) else np.nan for b in BEAMS}}


def run_stage(batch_key, chunk, prm, run):
    import geopandas as gpd
    from sliderule import sliderule, icesat2

    batch = load_batch(batch_key)
    P = paths_for(batch_key)
    tracks = pd.read_csv(P.batch_work + 'tracks.csv')
    chunks = pd.read_csv(P.batch_work + 'chunks.csv').set_index('chunk')
    ch = chunks.loc[chunk]
    todo = tracks[(tracks.chunk == chunk) & tracks.selected]
    expected_ids = list(todo.ID.unique())
    granules = list(todo.granule)
    print(f'chunk {chunk}: {ch.t0} .. {ch.t1}, {len(granules)} granules, {len(expected_ids)} expected tracks')
    if not granules:
        raise SkipTrack('no selected granules in chunk')

    Gtrack_lowest = gpd.read_file(P.batch_work + 'rgt_start_points.geojson')
    poly = sct.create_polygons(list(batch['region']['lat']), list(batch['region']['lon']))

    # %% SlideRule request: one call for all granules of the chunk, clipped to the polygon
    sl = batch['sliderule']
    print(f"sliderule.init(desired_nodes={sl['desired_nodes']}, time_to_live={sl['time_to_live']}) ... nodes take ~2-3 min")
    sliderule.init(desired_nodes=sl['desired_nodes'], time_to_live=sl['time_to_live'], verbose=True, user_service=True)
    params = dict(prm['sliderule'])
    params['poly'] = poly['list']
    params['t0'], params['t1'] = str(ch.t0), str(ch.t1)
    t_req = datetime.datetime.now()
    gdf = icesat2.atl06p(params, resources=granules)
    print(f'atl06p returned {len(gdf)} rows in {(datetime.datetime.now() - t_req).seconds}s')
    run.info(n_rows_raw=int(len(gdf)), request_s=(datetime.datetime.now() - t_req).seconds)
    if len(gdf) == 0:
        for ID in expected_ids:
            _write_track_status(run, ID, 'skip', 'no data returned for chunk')
        raise SkipTrack('atl06p returned no data')
    gdf = sct.correct_and_remove_height(gdf, prm['maximum_height'])

    # %% split per (rgt, cycle) -> tracks
    id_of = todo.drop_duplicates('ID').set_index(['rgt', 'cycle'])['ID']
    granules_of = todo.groupby('ID')['granule'].apply(list)
    produced, density = {}, {}
    for (rgt, cyc), tmp in gdf.groupby(['rgt', 'cycle']):
        key = (int(rgt), int(cyc))
        if key not in id_of.index:
            print(f'  rgt {rgt} cycle {cyc}: {len(tmp)} rows but not an expected track, skipped')
            continue
        ID = id_of.loc[key]
        print(f'  {ID}: {len(tmp)} rows')
        if len(tmp) < prm['min_points_per_track']:
            _write_track_status(run, ID, 'skip', f'too few points in box ({len(tmp)} < {prm["min_points_per_track"]})')
            continue
        try:
            info = track_products(ID, tmp, granules_of.loc[ID], P, prm, Gtrack_lowest)
            _write_track_status(run, ID, 'success', None, info=info)
            produced[ID] = info
            density[ID] = pd.Series(info['n_points'], name=ID)
        except SkipTrack as e:
            print(f'  {ID} skipped: {e.reason}')
            _write_track_status(run, ID, 'skip', e.reason, info=e.info)
        except Exception as e:            # one bad track must not sink the chunk
            import traceback
            print(f'  {ID} FAILED: {e!r}')
            traceback.print_exc()
            _write_track_status(run, ID, 'fail', None, error=e)
        plt.close('all')
    for ID in expected_ids:
        if read_status(batch_key, STAGE, ID) is None or read_status(batch_key, STAGE, ID)['t_start'] < run.t_start.isoformat():
            _write_track_status(run, ID, 'skip', 'no data in box for this track')

    # %% chunk summary
    if density:
        D = pd.concat(density, axis=1).T
        D['total'] = D.sum(axis=1)
        D.to_html(P.stage_dir('A01b_ID') + f'{batch_key}_chunk{chunk}_point_density.html')
        print(D.to_string())
    run.info(n_expected=len(expected_ids), n_success=len(produced),
             n_skip=sum(1 for ID in expected_ids if (read_status(batch_key, STAGE, ID) or {}).get('status') == 'skip'),
             n_fail=sum(1 for ID in expected_ids if (read_status(batch_key, STAGE, ID) or {}).get('status') == 'fail'))


def _write_track_status(run, ID, status, reason, info=None, error=None):
    """per-track status record derived from the chunk run"""
    rec = run._record(status)
    rec.update({'ID': ID, 'status': status, 'reason': reason, 'info': info or {},
                't_end': datetime.datetime.now().replace(microsecond=0).isoformat(),
                'runtime_s': None, 'log': f'logs/{STAGE}/chunk{run.ID[5:]}.log' if run.ID.startswith('chunk') else None,
                'outputs': [f'B01_regrid/{ID}_B01_binned.h5', f'A01b_ID/A01b_ID_{ID}.json'] if status == 'success' else [],
                'figures': ['B01b_ATL06_corrected.png', 'B01b_beam_statistics.png', 'B01_track.png'] if status == 'success' else []})
    if error is not None:
        import traceback
        rec.update({'error_type': type(error).__name__, 'error_msg': str(error)[:500],
                    'traceback_tail': traceback.format_exc().splitlines()[-30:]})
    write_status(rec)


def main(batch_key, chunk):
    prm = load_params(batch_key)
    section = dict(prm['B01'], version=prm['version'])
    with StageRun(STAGE, f'chunk{chunk}', batch_key, params=section, script=__file__) as run:
        run_stage(batch_key, int(chunk), prm['B01'], run)
    return run.status


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('-')]
    if len(args) < 2:
        args = ['SH_dev_small', '0']
        print('no arguments, using', args)
    sys.exit(0 if main(args[0], args[1]) != 'fail' else 1)
