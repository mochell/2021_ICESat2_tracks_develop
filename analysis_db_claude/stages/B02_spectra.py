# %%
"""
B02: generalized-Fourier (gFT) and FFT slope spectra per beam.

    python stages/B02_spectra.py <ID> <batch_key>

Reads  work/<batch>/B01_regrid/<ID>_B01_binned.h5                      (B01)
Writes work/<batch>/B02_spectra/B02_<ID>_gFT_k.nc, _gFT_x.nc, _FFT.nc, B02_<ID>_params.h5
       status/B02/<ID>.json

Port of analysis_db/B02_make_spectra_gFT.py: same algorithm, parameters from params/<v>.toml [B02],
exit() -> SkipTrack, the per-beam 'no data' branch is a real skip of that beam.
"""
import sys
import copy
import time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline_config import mconfig, np, pd, xr, plt, M, MT, col, paths_for, cli_args
from pipeline_status import StageRun, SkipTrack, require_upstream, read_status
from pipeline_params import load_params

import h5py
from threadpoolctl import threadpool_limits
from scipy.ndimage import label
import ICEsat2_SI_tools.io as io
import ICEsat2_SI_tools.spectral_estimates as spec
import spicke_remover
import generalized_FT as gFT

STAGE = 'B02'


def linear_gap_fill(F, key_lead, key_int):
    """linear interpolation over nans of column key_int of DataFrame F"""
    y_g = np.array(F[key_int])
    nans, x2 = np.isnan(y_g), lambda z: z.nonzero()[0]
    y_g[nans] = np.interp(x2(nans), x2(~nans), y_g[~nans])
    return y_g


def repack_attributes(DD):
    for k in list(DD.keys()):
        for ka in list(DD[k].attrs.keys()):
            I = DD[k]
            I.coords[ka] = ('beam', np.expand_dims(I.attrs[ka], 0))
    return DD


def dense_segments(x, xlims, Lmeters, dx, fill_min, window_m):
    """
    Largest number of 'dense' L-segments inside any window of window_m along the track.
    The track [xlims] is cut into consecutive segments of Lmeters; a segment is dense when at least
    fill_min of its dx slots hold a point (the gFT itself drops stancils below 40 %).
    """
    edges = np.arange(xlims[0], xlims[1] + Lmeters, Lmeters)
    if edges.size < 2 or x.size == 0:
        return 0
    counts, _ = np.histogram(x, edges)
    dense = (counts / (Lmeters / dx) >= fill_min).astype(int)
    n_win = max(int(window_m // Lmeters), 1)
    return int(np.convolve(dense, np.ones(n_win, dtype=int), mode='valid').max()) if dense.size >= n_win else int(dense.sum())


def make_dummy_beam(GG, beam):
    dummy = GG.copy(deep=True)
    for var in list(dummy.var()):
        dummy[var] = dummy[var] * np.nan
    dummy['beam'] = [beam]
    return dummy


def run_stage(ID, batch_key, prm, run):
    P = paths_for(batch_key, ID)
    require_upstream(run, ['B01'])
    all_beams = mconfig['beams']['all_beams']
    # beams B01 found usable (thin weak beams are left out; older B01 records have no list -> all)
    b01_info = (read_status(batch_key, 'B01', ID) or {}).get('info') or {}
    beams_ok = [b for b in all_beams if b in b01_info.get('beams_ok', all_beams)]
    run.info(beams_ok=beams_ok, beam_mode=b01_info.get('beam_mode', 'all'))

    load_path = P.stage_dir('B01_regrid')
    save_path = P.stage_dir('B02_spectra')
    save_name = 'B02_' + ID

    # %% load and quality tests
    Gd = h5py.File(load_path + ID + '_B01_binned.h5', 'r')

    nan_fraction = list()
    for k in beams_ok:
        xk = io.get_beam_var_hdf_store(Gd[k], 'x')
        nan_fraction.append(np.sum(np.isnan(xk)) / max(xk.shape[0], 1))
    nan_fraction = float(np.array(nan_fraction).mean())

    # strong/weak point ratio per complete pair; outside [1/r, r] the weak beam is too sparse -> drop it
    # (the track continues on the remaining beams, like a thin weak beam in B01)
    ratios, beams_dropped_ratio = {}, []
    for group in mconfig['beams']['groups']:
        if not set(group) <= set(beams_ok):       # incomplete pair (thin weak beam): nothing to compare
            continue
        na, nb = Gd[group[0]]['x'][:].size, Gd[group[1]]['x'][:].size
        ratio = na / nb if nb else np.inf
        ratios['/'.join(group)] = ratio
        if (ratio > prm['beam_ratio_max']) | (ratio < 1 / prm['beam_ratio_max']):
            weak = group[1] if na >= nb else group[0]
            print('bad data ratio ', group, ratio, '-> drop', weak)
            beams_dropped_ratio.append(weak)
    beams_ok = [b for b in beams_ok if b not in beams_dropped_ratio]
    run.info(nan_fraction=nan_fraction, beam_ratios=ratios, beams_dropped_ratio=beams_dropped_ratio, beams_used=beams_ok)
    if nan_fraction > prm['max_nan_fraction']:
        raise SkipTrack(f'nan fraction {nan_fraction:.2f} > {prm["max_nan_fraction"]}', nan_fraction=nan_fraction)

    # %% spectral limits
    dist = io.get_beam_var_hdf_store(Gd[beams_ok[0]], 'x')
    T_max = prm['T_max']
    k_0 = (2 * np.pi / T_max) ** 2 / 9.81
    x = np.array(dist).squeeze()
    dx = np.round(np.median(np.diff(x)), 1)
    min_datapoint = 2 * np.pi / k_0 / dx
    Lpoints = int(np.round(min_datapoint) * prm['L_factor'])
    Lmeters = Lpoints * dx
    print('L number of gridpoint:', Lpoints, ' L length in km:', Lmeters / 1e3)

    T_min = prm['T_min']
    lambda_min = 9.81 * T_min ** 2 / (2 * np.pi)
    dlambda = Lmeters * prm['oversample']
    kk = np.arange(0, 1 / lambda_min, 1 / dlambda) * 2 * np.pi
    kk = kk[k_0 <= kk]
    kk = kk[::prm['k_stride']]
    print('2 M = ', kk.size * 2)

    # %% global xlims over beams
    dist_list = np.array([np.nan, np.nan])
    for k in beams_ok:
        xk = Gd[k + '/x'][:]
        dist_list = np.vstack([dist_list, [xk[0], xk[-1]]])
    xlims = np.nanmin(dist_list[:, 0]) - dx, np.nanmin(dist_list[:, 1])
    print('xlims: ', xlims)
    run.info(dx=float(dx), Lmeters=float(Lmeters), Lpoints=int(Lpoints), n_k=int(kk.size),
             x_range_km=[float(xlims[0] / 1e3), float(xlims[1] / 1e3)])
    if xlims[1] - xlims[0] < Lmeters:
        raise SkipTrack(f'common x range {(xlims[1] - xlims[0]) / 1e3:.1f} km shorter than one stancil ({Lmeters / 1e3:.0f} km)')

    # %% per-beam gFT and FFT
    G_gFT, G_gFT_x, G_rar_fft, Pars_optm = dict(), dict(), dict(), dict()
    beams_skipped, spike_remover_failed = [b for b in all_beams if b not in beams_ok], []
    dense_count = {}
    hkey, hkey_sigma = 'h_mean', 'h_sigma'

    for k in beams_ok:
        Gi = io.get_beam_hdf_store(Gd[k])
        x_mask = (Gi['x'] > xlims[0]) & (Gi['x'] < xlims[1])
        n_dense = dense_segments(np.asarray(Gi['x'][x_mask]), xlims, Lmeters, dx, prm['segment_fill_min'], prm['dense_window_km'] * 1e3)
        dense_count[k] = n_dense
        print(k, 'dense segments in best window', n_dense)
        if n_dense < prm['min_dense_segments']:
            print('------------------- too few dense segments in beam; skip beam', k)
            beams_skipped.append(k)
            continue

        Gd_cut = Gi[x_mask]
        x = Gd_cut['x']
        del Gi
        x_mask = (x >= xlims[0]) & (x <= xlims[1])
        x = x[x_mask]
        dd = np.copy(Gd_cut[hkey])
        dd_error = np.copy(Gd_cut[hkey_sigma])
        dd_error[np.isnan(dd_error)] = prm['dd_error_fill']

        # slope spectra
        if dd.size < prm['min_points_beam']:
            print('------------------- too few points in beam', k, dd.size, '; skip beam')
            beams_skipped.append(k)
            continue
        dd = np.gradient(dd)
        try:
            dd, _ = spicke_remover.spicke_remover(dd, spreed=prm['spike_spreed'], verbose=False)
        except (ValueError, IndexError) as e:
            # spicke_remover cannot handle spikes at the very end / very short series; keep the raw slopes
            print('spike remover failed for beam', k, repr(e), '; using unfiltered slopes')
            spike_remover_failed.append(k)
        dd_nans = (np.isnan(dd)) + (Gd_cut['N_photos'] <= prm['N_photos_min'])

        dd_no_nans = dd[~dd_nans]
        x_no_nans = x[~dd_nans]
        dd_error_no_nans = dd_error[~dd_nans]
        if dd_no_nans.size < prm['min_points_beam']:
            print('------------------- too few valid points in beam', k, dd_no_nans.size, '; skip beam')
            beams_skipped.append(k)
            continue

        print('gFT', k)
        with threadpool_limits(limits=prm['n_threads'], user_api='blas'):
            S = gFT.wavenumber_spectrogram_gFT(np.array(x_no_nans), np.array(dd_no_nans), Lmeters, dx, kk,
                                               data_error=dd_error_no_nans, ov=None)
            try:
                GG, GG_x, Params = S.cal_spectrogram(xlims=xlims, max_nfev=prm['max_nfev'], plot_flag=False)
            except (StopIteration, ValueError) as e:
                # generalized_FT: not a single stancil of the beam produced a fit (StopIteration, or with
                # newer xarray a ValueError 'must supply at least one object to concatenate')
                if isinstance(e, ValueError) and 'at least one object' not in str(e):
                    raise
                print('------------------- no stancil converged in beam', k, '; skip beam')
                beams_skipped.append(k)
                continue

        S.parceval(add_attrs=True, weight_data=False)

        GG.coords['beam'] = GG_x.coords['beam'] = str(k)
        GG, GG_x = GG.expand_dims(dim='beam', axis=1), GG_x.expand_dims(dim='beam', axis=1)
        GG.coords['N_per_stancil'] = (('x', 'beam'), np.expand_dims(GG['N_per_stancil'], 1))
        GG.coords['spec_adjust'] = (('x', 'beam'), np.expand_dims(GG['spec_adjust'], 1))

        x_coord_no_gaps = linear_gap_fill(Gd_cut, 'x', 'x')
        y_coord_no_gaps = linear_gap_fill(Gd_cut, 'x', 'y')
        mapped_coords = np.atleast_2d(spec.sub_sample_coords(Gd_cut['x'], x_coord_no_gaps, y_coord_no_gaps, S.stancil_iter, map_func=None))
        if mapped_coords.ndim != 2 or mapped_coords.shape[1] != 3 or mapped_coords.shape[0] != GG.x.size:
            print('------------------- stancil/coordinate mismatch in beam', k, mapped_coords.shape, '; skip beam')
            beams_skipped.append(k)
            continue
        GG.coords['x_coord'] = GG_x.coords['x_coord'] = (('x', 'beam'), np.expand_dims(mapped_coords[:, 1], 1))
        GG.coords['y_coord'] = GG_x.coords['y_coord'] = (('x', 'beam'), np.expand_dims(mapped_coords[:, 2], 1))

        if (GG.coords['N_per_stancil'] == 0).squeeze()[0].data:
            nlabel = label((GG.coords['N_per_stancil'] == 0).squeeze())[0]
            nan_mask = nlabel == nlabel[0]
            GG.coords['x_coord'][nan_mask] = np.nan
            GG.coords['y_coord'][nan_mask] = np.nan

        lons_no_gaps = linear_gap_fill(Gd_cut, 'x', 'lons')
        lats_no_gaps = linear_gap_fill(Gd_cut, 'x', 'lats')
        mapped_coords = spec.sub_sample_coords(Gd_cut['x'], lons_no_gaps, lats_no_gaps, S.stancil_iter, map_func=None)
        GG.coords['lon'] = GG_x.coords['lon'] = (('x', 'beam'), np.expand_dims(mapped_coords[:, 1], 1))
        GG.coords['lat'] = GG_x.coords['lat'] = (('x', 'beam'), np.expand_dims(mapped_coords[:, 2], 1))

        def get_stancil_nans(stancil):
            m = (stancil[0] < x) & (x <= stancil[-1])
            return stancil[1], Gd_cut['N_photos'][m].sum()

        photon_list = np.array(list(dict(map(get_stancil_nans, copy.copy(S.stancil_iter))).values()))
        GG.coords['N_photons'] = (('x', 'beam'), np.expand_dims(photon_list, 1))

        G_gFT[k], G_gFT_x[k], Pars_optm[k] = GG, GG_x, Params

        # standard FFT
        print('FFT', k)
        dd[dd_nans] = 0
        S = spec.wavenumber_spectrogram(x, dd, Lpoints)
        try:
            G = S.cal_spectrogram()
        except ValueError as e:               # 'must supply at least one object to concatenate': beam shorter than one FFT chunk
            print('FFT failed for beam', k, repr(e), '; skip beam')
            G_gFT.pop(k), G_gFT_x.pop(k), Pars_optm.pop(k)
            beams_skipped.append(k)
            continue
        S.mean_spectral_error()
        S.parceval(add_attrs=True)
        G.coords['beam'] = str(k)
        G = G.expand_dims(dim='beam', axis=2)
        G.coords['mean_El'] = (('k', 'beam'), np.expand_dims(G['mean_El'], 1))
        G.coords['mean_Eu'] = (('k', 'beam'), np.expand_dims(G['mean_Eu'], 1))
        # FFT x is index*dx from the first data point; shift to the absolute track coordinate so it
        # can be cut to the gFT range below (the old script cut a relative x with absolute limits)
        G.coords['x'] = G.coords['x'] * dx + float(np.asarray(x)[0])
        G.attrs['x_absolute'] = 1

        stancil_iter = spec.create_chunk_boundaries(int(Lpoints), dd_nans.size)

        def get_stancil_N(stancil):
            idata = dd_nans[stancil[0]:stancil[-1]]
            return stancil[1], idata.size - idata.sum()

        N_list = np.array(list(dict(map(get_stancil_N, stancil_iter)).values()))
        G.coords['N_per_stancil'] = (('x', 'beam'), np.expand_dims(N_list, 1))
        try:
            G_rar_fft[k] = G.sel(x=slice(GG.x[0], GG.x[-1].data))
        except Exception:
            G_rar_fft[k] = G.isel(x=(GG.x[0].data < G.x.data) & (G.x.data < GG.x[-1].data))
        plt.close('all')

    Gd.close()
    run.info(dense_segments=dense_count)
    if not G_gFT:
        dense_ok = [k for k, n in dense_count.items() if n >= prm['min_dense_segments']]
        reason = (f'no beam with {prm["min_dense_segments"]} dense segments in {prm["dense_window_km"]} km' if not dense_ok
                  else 'no beam produced a gFT spectrum (dense beams failed later: too few valid points / no stancil fit)')
        raise SkipTrack(reason, beams_skipped=beams_skipped, dense_segments=dense_count)
    run.info(beams_skipped=beams_skipped, spike_remover_failed=spike_remover_failed, n_x=int(list(G_gFT.values())[0].x.size))

    # %% fill missing beams with nan dummies, save
    MT.save_pandas_table(Pars_optm, save_name + '_params', save_path)

    for beam in set(all_beams) - set(G_gFT.keys()):
        GG = list(G_gFT.values())[0]
        dummy = make_dummy_beam(GG, beam)
        dummy['N_photons'] = dummy['N_photons'] * 0
        dummy['N_per_stancil'] = dummy['N_per_stancil'] * 0
        G_gFT[beam] = dummy
        G_gFT_x[beam] = make_dummy_beam(list(G_gFT_x.values())[0], beam)
        GG = list(G_rar_fft.values())[0].copy(deep=True)
        GG.data = GG.data * np.nan
        GG['beam'] = [beam]
        G_rar_fft[beam] = GG

    G_gFT, G_gFT_x, G_rar_fft = repack_attributes(G_gFT), repack_attributes(G_gFT_x), repack_attributes(G_rar_fft)

    G_gFT_DS = xr.merge(G_gFT.values())
    G_gFT_DS['Z_hat_imag'] = G_gFT_DS.Z_hat.imag
    G_gFT_DS['Z_hat_real'] = G_gFT_DS.Z_hat.real
    G_gFT_DS = G_gFT_DS.drop_vars('Z_hat')
    G_gFT_DS.attrs['name'] = 'gFT_estimates'
    G_gFT_DS.to_netcdf(save_path + save_name + '_gFT_k.nc')
    # diagnostics: stancils whose PSD blew up (seen on the last stancil of the test track, ~1e19)
    psd_max_x = G_gFT_DS.gFT_PSD_data.max(('k', 'beam'), skipna=True)
    run.info(n_psd_blowup=int((psd_max_x > prm['psd_blowup']).sum()), psd_max=float(psd_max_x.max()))

    G_gFT_x_DS = xr.merge(G_gFT_x.values())
    G_gFT_x_DS.attrs['name'] = 'gFT_estimates_real_space'
    G_gFT_x_DS.to_netcdf(save_path + save_name + '_gFT_x.nc')

    G_fft_DS = xr.merge(G_rar_fft.values())
    G_fft_DS.attrs['name'] = 'FFT_power_spectra'
    G_fft_DS.to_netcdf(save_path + save_name + '_FFT.nc')
    print('saved and done')


def main(ID, batch_key):
    prm = load_params(batch_key)
    section = dict(prm[STAGE], version=prm['version'])
    with StageRun(STAGE, ID, batch_key, params=section, script=__file__) as run:
        run_stage(ID, batch_key, prm[STAGE], run)
    return run.status


if __name__ == '__main__':
    ID, batch_key = cli_args(sys.argv, default=('SH_20190502_05180312', 'SH_testSLsinglefile2'))
    sys.exit(0 if main(ID, batch_key) != 'fail' else 1)
