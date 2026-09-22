# %%
"""
A02: WW3 (IOWAGA THREDDS) wave prior around the equatorward end of the track.

    python stages/A02_ww3_prior.py <ID> <batch_key>

Reads  work/<batch>/A01b_ID/A01b_ID_<ID>.json                  (track time)
       work/<batch>/B01_regrid/<ID>_B01_binned.h5              (B01, first point per beam)
       IOWAGA-WW3-FORECAST GLOB-30M via OPeNDAP (siphon)
Writes work/<batch>/A02_prior/A02_<ID>.h5   table 'priors_hindcast'  (read by B04/B05 with
       MT.load_pandas_table_dict('/A02_'+ID, path)['priors_hindcast'])
       plots/<hemis>/<batch>/<ID>/A02_hindcast_data, A02_hindcast_prior
       status/A02/<ID>.json

Port of analysis_db/A02c_IOWAGA_thredds_prior.py: same algorithm, parameters from params/<v>.toml
[A02] and [A02.<hemis>], no bare try/except (errors propagate to StageRun), the latitude-shifting
loop of the prior box is capped (-> SkipTrack), no _hindcast_success/_fail.json markers.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline_config import mconfig, np, pd, xr, plt, GridSpec, M, MT, col, paths_for, save_fig, font_for_print, cli_args
from pipeline_status import StageRun, SkipTrack, require_upstream
from pipeline_params import load_params

import h5py
from siphon.catalog import TDSCatalog
import ICEsat2_SI_tools.io as io
import ICEsat2_SI_tools.wave_tools as waves

STAGE = 'A02'
DATA_LABEL = 'LOPS WW3-GLOB-30M'     # figure label of the hindcast product

VAR_LIST = ['dir', 'dp', 'fp', 'hs', 'ice', 'spr',
            't01', 't02',
            'plp0',
            'pdir0', 'pdir1', 'pdir2', 'pdir3', 'pdir4', 'pdir5',
            'pspr0', 'pspr1', 'pspr2', 'pspr3', 'pspr4', 'pspr5',
            'ptp0', 'ptp1', 'ptp2', 'ptp3', 'ptp4', 'ptp5',
            'phs0', 'phs1', 'phs2', 'phs3', 'phs4', 'phs5']

# (amplitude, angle) pairs that are averaged in vector space
KEY_LIST_PAIRS = {
    'mean': ('hs', 'dir'),
    'peak': ('hs', 'dp'),
    'partion0': ('phs0', 'pdp0'),
    'partion1': ('phs1', 'pdp1'),
    'partion2': ('phs2', 'pdp2'),
    'partion3': ('phs3', 'pdp3'),
    'partion4': ('phs4', 'pdp4')}


def sel_data(I, lon_range, lat_range, timestamp=None, time_range=None):
    """
    this method returns the selected data in the lon-lat box at an interpolated timestamp
    """
    lon_flag = (lon_range[0] < I.longitude.data) & (I.longitude.data < lon_range[1])
    lat_flag = (lat_range[0] < I.latitude.data) & (I.latitude.data < lat_range[1])

    if timestamp is None:
        I = I.isel(latitude=lat_flag, longitude=lon_flag)
    else:
        target_time = np.datetime64(timestamp, 'ns')
        time_flag = (time_range[0] < I.time.data) & (I.time.data < time_range[1])
        I = I.isel(latitude=lat_flag, longitude=lon_flag, time=time_flag).sortby('time')
        # Force time coordinate to same precision as target_time
        I = I.assign_coords(time=I.time.values.astype('datetime64[ns]'))
        I = I.interp(time=target_time)
    return I


def draw_range(lon_range, lat_range, *args, **kargs):
    plt.plot([lon_range[0], lon_range[1], lon_range[1], lon_range[0], lon_range[0]],
             [lat_range[0], lat_range[0], lat_range[1], lat_range[1], lat_range[0]], *args, **kargs)


def test_nan_frac(imask, nan_frac_max):
    "True if the fraction of ice-free (False) cells in the mask is below nan_frac_max"
    return ((~imask).sum() / imask.size).data < nan_frac_max


def plot_prior(Prior, axx, lon_range):
    angle = Prior['incident_angle']['value']  # incident direction in degrees from North clockwise (Meteorological convention)
    angle_plot = - angle - 90
    axx.quiver(Prior['center_lon']['value'], Prior['center_lat']['value'],
               - np.cos(angle_plot * np.pi / 180), - np.sin(angle_plot * np.pi / 180),
               scale=4.5, zorder=12, width=0.1, headlength=4.5, minshaft=2, alpha=0.6, color='black')
    axx.plot(Prior['center_lon']['value'], Prior['center_lat']['value'], '.', markersize=6, zorder=12, alpha=1, color='black')
    tstring = (' ' + str(np.round(Prior['peak_period']['value'], 1)) + 'sec \n ' + str(np.round(Prior['Hs']['value'], 1))
               + 'm\n ' + str(np.round(angle, 1)) + 'deg')
    plt.text(lon_range[1], Prior['center_lat']['value'], tstring)


def run_stage(ID, batch_key, prm, run):
    P = paths_for(batch_key, ID)
    require_upstream(run, ['B01'])
    track_name = ID
    hemis = P.hemis
    hprm = prm[hemis]
    all_beams = mconfig['beams']['all_beams']

    save_path = P.stage_dir('A02_prior')
    plot_path = P.track_plot_dir()
    save_name = 'A02_' + track_name
    col.colormaps2(21)

    # %% track time from the A01b ID file, first point per beam from B01
    ID_json = MT.json_load('A01b_ID_' + track_name, P.stage_dir('A01b_ID', mkdir=False))

    load_path = P.stage_dir('B01_regrid', mkdir=False)
    Gd = h5py.File(load_path + track_name + '_B01_binned.h5', 'r')
    G1 = dict()
    for b in all_beams:
        Gi = io.get_beam_hdf_store(Gd[b])
        G1[b] = Gi.iloc[abs(Gi['lats']).argmin()]      # equatorward end of the beam
    Gd.close()
    G1 = pd.DataFrame.from_dict(G1).T

    # %% DEFINE SEARCH REGION AND SIZE OF BOXES FOR AVERAGES
    dlon_deg = hprm['dlon_deg']               # lon degree range around 1st point
    dlat_deg = hprm['dlat_deg']               # lat degree range around 1st point
    dlat_deg_prior = hprm['dlat_deg_prior']   # degree range around 1st point of the prior box
    dtime = prm['dtime_h']                    # in hours

    lon_range = G1['lons'].min() - dlon_deg, G1['lons'].max() + dlon_deg
    if hemis == 'SH':
        lat_range = np.sign(G1['lats'].min()) * 78, G1['lats'].max() + dlat_deg[1]
    else:
        lat_range = G1['lats'].min() - dlat_deg[0], G1['lats'].max() + dlat_deg[1]
    lat_range_prior = G1['lats'].min() - dlat_deg_prior[0], G1['lats'].max() + dlat_deg_prior[1]

    timestamp = pd.to_datetime(ID_json['pars']['start']['delta_time'], unit='s')
    time_range = np.datetime64(timestamp) - np.timedelta64(dtime, 'h'), np.datetime64(timestamp) + np.timedelta64(dtime, 'h')
    print('timestamp', timestamp, 'lon_range', lon_range, 'lat_range', lat_range, 'lat_range_prior', lat_range_prior)

    # %% load WW3 data (ECMWF forecast) from the IOWAGA THREDDS server
    cat = TDSCatalog(prm['catalog_url'])
    ncss = cat.datasets[prm['dataset']].remote_access(use_xarray=True)

    IOWAGA = ncss[VAR_LIST]
    IOWAGA['time'] = np.array([np.datetime64(k0) for k0 in IOWAGA.time.data]).astype('M8[h]')
    IOWAGA = IOWAGA.rename(name_dict={'pdir0': 'pdp0', 'pdir1': 'pdp1', 'pdir2': 'pdp2',
                                      'pdir3': 'pdp3', 'pdir4': 'pdp4', 'pdir5': 'pdp5'})

    G_beam = sel_data(IOWAGA, lon_range, lat_range, timestamp, time_range).load()
    G_prior = sel_data(G_beam, lon_range, lat_range_prior)

    # %% ice mask
    if hemis == 'SH':
        ice_mask = (G_beam.ice > 0) | np.isnan(G_beam.ice)

        # mask all latitudes that are completely full with sea ice.
        lats = list(ice_mask.latitude.data)
        lats.sort(reverse=True)
        # find 1st latitude that is completely full with sea ice.
        ice_lat_pos = next((i for i, j in enumerate((ice_mask.sum('longitude') == ice_mask.longitude.size).sel(latitude=lats)) if j), None)
        # recreate lat mask based on this criteria; no fully ice-covered row -> nothing masked
        # (the old code crashed here with lats[None])
        lats = np.array(lats)
        lat_mask = lats < lats[ice_lat_pos] if ice_lat_pos is not None else np.zeros(lats.size, dtype=bool)
        run.info(ice_edge_row_found=ice_lat_pos is not None)
        lat_mask = xr.DataArray(lat_mask.repeat(ice_mask.longitude.size).reshape(ice_mask.shape), dims=ice_mask.dims, coords=ice_mask.coords)
        lat_mask['latitude'] = lats

        # combine ice mask and new lat mask
        ice_mask = ice_mask + lat_mask
    else:
        ice_mask = np.isnan(G_beam.ice)
        lats = ice_mask.latitude
        # find closest latitude with non-nan data
        ice_lat_pos = abs(lats.where(ice_mask.sum('longitude') > 4, np.nan) - np.array(lat_range).mean()).argmin().data
        # redefine lat-range
        lat_range = lats[ice_lat_pos].data - 2, lats[ice_lat_pos].data + 2

    # %% plot 1st figure: hindcast fields in the search box
    dir_clev = np.arange(0, 380, 20)
    f_clev = np.arange(1 / 40, 1 / 5, 0.01)
    fvar = ['ice', 'dir', 'dp', 'spr', 'fp', 'hs']
    fcmap = [plt.cm.Blues_r, col.circle_medium_triple, col.circle_medium_triple, plt.cm.Blues, plt.cm.Blues, plt.cm.Blues]
    fpos = [0, 1, 2, 3, 4, 5]
    clevs = [np.arange(0, 1, 0.2), dir_clev, dir_clev, np.arange(0, 90, 10), f_clev, np.arange(.5, 9, 0.5)]

    font_for_print()
    F = M.figure_axis_xy(4, 3.5, view_scale=0.9, container=True)
    plt.suptitle(track_name + ' | ' + DATA_LABEL, y=1.3)
    lon, lat = G_beam.longitude, G_beam.latitude

    gs = GridSpec(9, 6, wspace=0.1, hspace=0.4)
    for fv, fp, fc, cl in zip(fvar, fpos, fcmap, clevs):
        ax1 = F.fig.add_subplot(gs[0:7, fp])
        if fp == 0:
            ax1.spines['bottom'].set_visible(False)
            ax1.spines['left'].set_visible(False)
            ax1.tick_params(labelbottom=True, bottom=True)
        else:
            ax1.axis('off')

        plt.plot(G1['lons'], G1['lats'], '.r', markersize=5)
        draw_range(lon_range, lat_range_prior, c='red', linewidth=1, zorder=12)
        draw_range(lon_range, lat_range, c='blue', linewidth=0.7, zorder=10)
        if fv != 'ice':
            cm = plt.pcolor(lon, lat, G_beam[fv], vmin=cl[0], vmax=cl[-1], cmap=fc)
            if G_beam.ice.shape[0] > 1:
                plt.contour(lon, lat, G_beam.ice, colors='black', linewidths=0.6)
        else:
            cm = plt.pcolor(lon, lat, G_beam[fv], vmin=cl[0], vmax=cl[-1], cmap=fc)

        plt.title(G_beam[fv].long_name.replace(' ', '\n') + '\n' + fv, loc='left')
        ax1.axis('equal')

        ax2 = F.fig.add_subplot(gs[-1, fp])
        cbar = plt.colorbar(cm, cax=ax2, orientation='horizontal', aspect=1, fraction=1)
        cl_ticks = np.linspace(cl[0], cl[-1], 3)
        cbar.set_ticks(np.round(cl_ticks, 3))
        cbar.set_ticklabels(np.round(cl_ticks, 2))

    save_fig(F, plot_path, 'A02_hindcast_data')

    # %% derive prior: shift the prior box poleward until enough ice-free cells are inside
    ice_mask_prior = ice_mask.sel(latitude=G_prior.latitude)

    n_shifts = 0
    while test_nan_frac(ice_mask_prior, prm['nan_frac_max']):
        if n_shifts >= prm['lat_shift_max_iter']:
            raise SkipTrack(f'prior box is all ice/nan after {n_shifts} shifts',
                            n_lat_shifts=n_shifts, lat_range_prior=[float(lat_range_prior[0]), float(lat_range_prior[1])])
        print(lat_range_prior)
        lat_range_prior = lat_range_prior[0] + 0.5, lat_range_prior[1] + 0.5
        G_prior = sel_data(G_beam, lon_range, lat_range_prior)
        ice_mask_prior = ice_mask.sel(latitude=G_prior.latitude)
        n_shifts += 1

    G_prior_masked = G_prior.where(~ice_mask_prior, np.nan)
    ice_fraction = float(ice_mask_prior.sum() / ice_mask_prior.size)

    # %% make pandas table with obs track end positions
    key_list = list(G_prior_masked.keys())
    key_list_pairs2 = list()
    for k in KEY_LIST_PAIRS.values():
        key_list_pairs2.append(k[0])
        key_list_pairs2.append(k[1])
    key_list_scaler = set(key_list) - set(key_list_pairs2)

    # derive angle average
    Tend = pd.DataFrame(index=key_list, columns=['mean', 'std', 'name'])
    for k, pair in KEY_LIST_PAIRS.items():
        ave_amp, ave_deg, std_amp, std_deg = waves.get_ave_amp_angle(G_prior_masked[pair[0]].data, G_prior_masked[pair[1]].data)
        Tend.loc[pair[0]] = ave_amp, std_amp, G_prior_masked[pair[0]].long_name
        Tend.loc[pair[1]] = ave_deg, std_deg, G_prior_masked[pair[1]].long_name

    for k in key_list_scaler:
        Tend.loc[k] = G_prior_masked[k].mean().data, G_prior_masked[k].std().data, G_prior_masked[k].long_name

    Tend = Tend.T
    Tend['lon'] = [ice_mask_prior.longitude.mean().data, ice_mask_prior.longitude.std().data, 'lontigude']
    Tend['lat'] = [ice_mask_prior.latitude[ice_mask_prior.sum('longitude') == 0].mean().data,
                   ice_mask_prior.latitude[ice_mask_prior.sum('longitude') == 0].std().data, 'latitude']
    Tend = Tend.T

    Prior = dict()
    Prior['incident_angle'] = {'value': Tend['mean']['dp'].astype('float'), 'name': Tend['name']['dp']}
    Prior['spread'] = {'value': Tend['mean']['spr'].astype('float'), 'name': Tend['name']['spr']}
    Prior['Hs'] = {'value': Tend['mean']['hs'].astype('float'), 'name': Tend['name']['hs']}
    Prior['peak_period'] = {'value': 1 / Tend['mean']['fp'].astype('float'), 'name': '1/' + Tend['name']['fp']}
    Prior['center_lon'] = {'value': Tend['mean']['lon'].astype('float'), 'name': Tend['name']['lon']}
    Prior['center_lat'] = {'value': Tend['mean']['lat'].astype('float'), 'name': Tend['name']['lat']}

    MT.save_pandas_table({'priors_hindcast': Tend}, save_name, save_path)

    run.info(timestamp=str(timestamp),
             lon_range=[float(lon_range[0]), float(lon_range[1])],
             lat_range=[float(lat_range[0]), float(lat_range[1])],
             lat_range_prior=[float(lat_range_prior[0]), float(lat_range_prior[1])],
             n_lat_shifts=n_shifts, ice_fraction_prior=ice_fraction,
             prior_hs_m=float(Prior['Hs']['value']),
             prior_fp_hz=float(Tend['mean']['fp']), prior_tp_s=float(Prior['peak_period']['value']),
             prior_dp_deg=float(Prior['incident_angle']['value']), prior_dir_deg=float(Tend['mean']['dir']),
             prior_spread_deg=float(Prior['spread']['value']),
             prior_center=[float(Prior['center_lon']['value']), float(Prior['center_lat']['value'])])

    # %% plot 2nd figure: prior (optional, must not fail the stage)
    try:
        font_for_print()
        F = M.figure_axis_xy(2, 4.5, view_scale=0.9, container=False)
        ax1 = F.ax
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        ax1.tick_params(labelbottom=True, bottom=True)

        plot_prior(Prior, ax1, lon_range)

        str_list = list()
        for i in np.arange(0, 6):
            str_list.append(' ' + str(np.round(Tend.loc['ptp' + str(i)]['mean'], 1)) + 'sec\n ' + str(np.round(Tend.loc['phs' + str(i)]['mean'], 1))
                            + 'm ' + str(np.round(Tend.loc['pdp' + str(i)]['mean'], 1)) + 'd')
        plt.text(lon_range[1], lat_range[0], '\n '.join(str_list))

        for vv in zip(['pdp0', 'pdp1', 'pdp2', 'pdp3', 'pdp4', 'pdp5'], ['phs0', 'phs1', 'phs3', 'phs4', 'phs5']):
            angle_plot = - Tend.loc[vv[0]]['mean'] - 90
            vsize = (1 / Tend.loc[vv[1]]['mean']) ** (1 / 2) * 5
            ax1.quiver(Prior['center_lon']['value'], Prior['center_lat']['value'],
                       - np.cos(angle_plot * np.pi / 180), - np.sin(angle_plot * np.pi / 180),
                       scale=vsize, zorder=5, width=0.1, headlength=4.5, minshaft=4, alpha=0.6, color='green')

        plt.plot(G1['lons'], G1['lats'], '.r', markersize=5)
        draw_range(lon_range, lat_range_prior, c='red', linewidth=1, zorder=11)
        draw_range(lon_range, lat_range, c='blue', linewidth=0.7, zorder=10)
        plt.pcolor(lon, lat, G_beam['ice'], cmap=fcmap[-1])

        plt.title('Prior\n' + DATA_LABEL + '\n' + track_name + '\nIncident angle', loc='left')
        ax1.axis('equal')

        save_fig(F, plot_path, 'A02_hindcast_prior')
    except Exception as e:
        print('2nd figure (A02_hindcast_prior) failed:', repr(e))
        plt.close('all')

    print('done')


def main(ID, batch_key):
    prm = load_params(batch_key)
    section = dict(prm[STAGE], version=prm['version'])
    with StageRun(STAGE, ID, batch_key, params=section, script=__file__) as run:
        run_stage(ID, batch_key, prm[STAGE], run)
    return run.status


if __name__ == '__main__':
    ID, batch_key = cli_args(sys.argv, default=('SH_20190502_05180312', 'SH_testSLsinglefile2'))
    sys.exit(0 if main(ID, batch_key) != 'fail' else 1)
