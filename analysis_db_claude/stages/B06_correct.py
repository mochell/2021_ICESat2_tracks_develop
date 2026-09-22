# %%
"""
B06: noise cut-off wavenumber per stancil, low-pass reconstruction of the gFT height model,
residual heights and incident-angle correction of the wavenumber / distance axes.

    python stages/B06_correct.py <ID> <batch_key>

Reads  work/<batch>/B01_regrid/<ID>_B01_binned.h5                                  (B01)
       work/<batch>/B02_spectra/B02_<ID>_gFT_k.nc, _gFT_x.nc, _FFT.nc                (B02)
       work/<batch>/B04_angle/B05_<ID>_angle_pdf.nc                                  (B05)
Writes work/<batch>/B06_corrected_separated/B06_<ID>_gFT_k_corrected.nc, B06_<ID>_gFT_x_corrected.nc,
       B06_<ID>_binned_resid.h5
       status/B06/<ID>.json
       plots/<hemis>/<batch>/<ID>/B06_correction/B06_k_cutoff, B06_atten_ov_simple, B06_atten_ov, B06_angle_def

Port of analysis_db/B06_correct_separate_var.py: same algorithm, parameters from params/<v>.toml [B06],
missing B05 angle file -> SkipTrack (was a bare try/except that silently set theta = 0),
B06_<ID>_B06_corrected_resid.h5 dropped (was always empty), B06_success.json dropped.
"""
import os
import sys
import copy
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline_config import (mconfig, np, pd, xr, plt, GridSpec, M, MT, col, lstrings, fig_sizes,
                             font_for_print, font_for_pres, paths_for, save_fig, cli_args)
from pipeline_status import StageRun, SkipTrack, require_upstream
from pipeline_params import load_params

import h5py
import piecewise_regression
import ICEsat2_SI_tools.io as io
import generalized_FT as gFT

STAGE = 'B06'


# ----------------------------------------------------------------------------- cut-off wavenumber
def get_correct_breakpoint(pw_results):
    """pick the breakpoint at which the steepest negative slope starts ('start' if no negative slope)"""
    br_points = list()
    for i in pw_results.keys():
        [br_points.append(i) if 'breakpoint' in i else None]
    br_points_df = pw_results[br_points]
    br_points_sorted = br_points_df.sort_values()

    alphas_sorted = [i.replace('breakpoint', 'alpha') for i in br_points_df.sort_values().index]
    alphas_sorted.append('alpha' + str(len(alphas_sorted) + 1))

    betas_sorted = [i.replace('breakpoint', 'beta') for i in br_points_df.sort_values().index]

    alphas_v2 = list()
    alpha_i = pw_results['alpha1']
    for i in [0] + list(pw_results[betas_sorted]):
        alpha_i += i
        alphas_v2.append(alpha_i)

    alphas_v2_sorted = pd.Series(index=alphas_sorted, data=alphas_v2)
    br_points_sorted['breakpoint' + str(br_points_sorted.size + 1)] = 'end'

    print('all alphas')
    print(alphas_v2_sorted)
    slope_mask = alphas_v2_sorted < 0

    if sum(slope_mask) == 0:
        print('no negative slope found, set to lowest')
        breakpoint = 'start'
    else:
        # take steepest slope
        alpah_v2_sub = alphas_v2_sorted[slope_mask]
        print(alpah_v2_sub)
        print(alpah_v2_sub.argmin())
        break_point_name = alpah_v2_sub.index[alpah_v2_sub.argmin()].replace('alpha', 'breakpoint')
        breakpoint = br_points_sorted[break_point_name]

    return breakpoint


def get_breakingpoints(xx, dd, n_breakpoints, n_breakpoints_max):
    """
    piecewise linear fit in log-log space; retries with one more breakpoint until it converges or
    n_breakpoints_max is reached. Returns (pw_fit, breakpoint) with breakpoint False if not converged.
    """
    x2, y2 = xx, dd
    convergence_flag = True
    while convergence_flag:
        pw_fit = piecewise_regression.Fit(x2, y2, n_breakpoints=n_breakpoints)
        print('n_breakpoints', n_breakpoints, pw_fit.get_results()['converged'])
        convergence_flag = not pw_fit.get_results()['converged']
        n_breakpoints += 1
        if n_breakpoints >= n_breakpoints_max:
            convergence_flag = False

    pw_results = pw_fit.get_results()

    if pw_results['converged']:
        pw_results_df = pd.DataFrame(pw_results['estimates']).loc['estimate']
        breakpoint = get_correct_breakpoint(pw_results_df)
        return pw_fit, breakpoint
    else:
        return pw_fit, False


def define_noise_wavenumber_piecewise(data_xr, n_breakpoints, n_breakpoints_max, plot_flag=False):
    """
    returns (breakpoint_k, pw_fit): the wavenumber at which the (log) displacement spectrum data_xr
    starts to decay steepest, from a piecewise regression in log-log space.
    """
    data_log = np.log(data_xr)

    k = data_log.k.data
    k_log = np.log(k)

    pw_fit, breakpoint_log = get_breakingpoints(k_log, data_log.data, n_breakpoints, n_breakpoints_max)

    # 'is' identity checks on strings are unreliable (the value comes back from a pandas Series)
    if isinstance(breakpoint_log, str) and breakpoint_log == 'start':
        print('no decay, set to lowerst wavenumber')
        breakpoint_log = k_log[0]
    if (isinstance(breakpoint_log, str) and breakpoint_log == 'end') or (breakpoint_log is False):
        print('higest wavenumner')
        breakpoint_log = k_log[-1]

    breakpoint_pos = abs(k_log - breakpoint_log).argmin()
    breakpoint_k = k[breakpoint_pos]

    if plot_flag:
        pw_fit.plot()
        plt.plot(k_log, data_log)

    return breakpoint_k, pw_fit


# ----------------------------------------------------------------------------- reconstruction
def tanh_filter(x, x_cutoff, sigma_g=0.01):
    """smooth low-pass: 1 below x_cutoff, 0 above, tanh transition of width sigma_g"""
    decay = 0.5 - np.tanh((x - x_cutoff) / sigma_g) / 2
    return decay


def reconstruct_displacement(Gx_1, Gk_1, k_thresh, sigma_g, dx):
    """
    reconstructs photon displacement heights for one stancil given the model parameters in Gk_1,
    low-pass filtered with a tanh filter at k_thresh.

    inputs:
    Gk_1     model data per stencil from _gFT_k file with sin and cos coefficients
    Gx_1     real data per stencil from _gFT_x file with mean photon heights and coordinate systems
    k_thresh threshold wavenumber for the low-pass filter
    sigma_g  width of the tanh transition
    dx       eta grid spacing of Gx

    returns:
    height_model  reconstructed displacement heights of the stancil
    dist_nanmask  mask where no observed data is
    """
    gFT_cos_coeff_sel = np.copy(Gk_1.gFT_cos_coeff)
    gFT_sin_coeff_sel = np.copy(Gk_1.gFT_sin_coeff)

    gFT_cos_coeff_sel = gFT_cos_coeff_sel * tanh_filter(Gk_1.k, k_thresh, sigma_g=sigma_g)
    gFT_sin_coeff_sel = gFT_sin_coeff_sel * tanh_filter(Gk_1.k, k_thresh, sigma_g=sigma_g)

    FT_int = gFT.generalized_Fourier(Gx_1.eta + Gx_1.x, None, Gk_1.k)
    _ = FT_int.get_H()
    FT_int.p_hat = np.concatenate([-gFT_sin_coeff_sel / Gk_1.k, gFT_cos_coeff_sel / Gk_1.k])

    height_model = FT_int.model() / dx

    dist_nanmask = np.isnan(Gx_1.y_data)

    return height_model, dist_nanmask


def save_pandas_table_overwrite(table_dict, name, save_path):
    """io.save_pandas_table appends to an existing HDFStore; remove the old file first"""
    f = save_path + name + '.h5'
    if os.path.exists(f):
        os.remove(f)
    io.save_pandas_table(table_dict, name, save_path)


# ----------------------------------------------------------------------------- stage
def run_stage(ID, batch_key, prm, run):
    P = paths_for(batch_key, ID)
    require_upstream(run, ['B02', 'B05'])
    ID_name = ID

    all_beams = mconfig['beams']['all_beams']
    high_beams = mconfig['beams']['high_beams']
    low_beams = mconfig['beams']['low_beams']

    load_path_regrid = P.stage_dir('B01_regrid')
    load_path_spectra = P.stage_dir('B02_spectra')
    load_path_angle = P.stage_dir('B04_angle')
    save_path = P.stage_dir('B06_corrected_separated')
    plot_path = P.track_plot_dir('B06_correction')

    # %% load data
    B3_hdf5 = h5py.File(load_path_regrid + ID_name + '_B01_binned.h5', 'r')
    B3 = dict()
    for b in all_beams:
        B3[b] = io.get_beam_hdf_store(B3_hdf5[b])
    B3_hdf5.close()

    load_file = load_path_spectra + 'B02_' + ID_name
    Gk = xr.open_dataset(load_file + '_gFT_k.nc')
    Gx = xr.open_dataset(load_file + '_gFT_x.nc')
    Gfft = xr.open_dataset(load_file + '_FFT.nc')

    angle_file = load_path_angle + 'B05_' + ID_name + '_angle_pdf.nc'
    if not os.path.exists(angle_file):
        raise SkipTrack('no B05 angle pdf', angle_file=angle_file)

    col.colormaps2(31, gamma=1)
    col_dict = col.rels

    # %% weighted means over beams
    G_gFT_wmean = (Gk.where(~np.isnan(Gk['gFT_PSD_data']), 0) * Gk['N_per_stancil']).sum('beam') / Gk['N_per_stancil'].sum('beam')
    G_gFT_wmean['N_photons'] = Gk['N_photons'].sum('beam')

    G_fft_wmean = (Gfft.where(~np.isnan(Gfft), 0) * Gfft['N_per_stancil']).sum('beam') / Gfft['N_per_stancil'].sum('beam')
    G_fft_wmean['N_per_stancil'] = Gfft['N_per_stancil'].sum('beam')

    # %% derive spectral errors
    Lpoints = Gk.Lpoints.mean('beam').data
    N_per_stancil = Gk.N_per_stancil.mean('beam').data

    G_error_model = dict()
    G_error_data = dict()

    for bb in Gk.beam.data:
        I = Gk.sel(beam=bb)
        b_bat_error = np.concatenate([I.model_error_k_cos.data, I.model_error_k_sin.data])
        Z_error = gFT.complex_represenation(b_bat_error, Gk.k.size, Lpoints)
        PSD_error_data, PSD_error_model = gFT.Z_to_power_gFT(Z_error, np.diff(Gk.k)[0], N_per_stancil, Lpoints)

        G_error_model[bb] = xr.DataArray(data=PSD_error_model, coords=I.drop_vars('N_per_stancil').coords, name='gFT_PSD_data_error').expand_dims('beam')
        G_error_data[bb] = xr.DataArray(data=PSD_error_data, coords=I.drop_vars('N_per_stancil').coords, name='gFT_PSD_data_error').expand_dims('beam')

    # NOTE: the original assigned the model error first and then overwrote it with the data error;
    # only the data error ends up in the files. G_error_model is kept for reference but not written.
    # join/coords/compat are passed explicitly (= the current xarray defaults) to silence the FutureWarnings
    # about changing defaults while keeping the result identical.
    gFT_PSD_data_error_mean = xr.concat(G_error_data.values(), dim='beam', join='outer', coords='different', compat='equals')
    gFT_PSD_data_error_mean = (gFT_PSD_data_error_mean.where(~np.isnan(gFT_PSD_data_error_mean), 0) * Gk['N_per_stancil']).sum('beam') / Gk['N_per_stancil'].sum('beam')

    G_gFT_wmean['gFT_PSD_data_err'] = gFT_PSD_data_error_mean
    Gk['gFT_PSD_data_err'] = xr.concat(G_error_data.values(), dim='beam', join='outer', coords='different', compat='equals')

    # %% smoothed weighted-mean spectrum
    G_gFT_smth = G_gFT_wmean['gFT_PSD_data'].rolling(k=prm['smooth_k'], center=True, min_periods=1).mean()
    G_gFT_smth['N_photons'] = G_gFT_wmean.N_photons
    G_gFT_smth["N_per_stancil_fraction"] = Gk['N_per_stancil'].T.mean('beam') / Gk.Lpoints.mean('beam')

    k = G_gFT_smth.k

    F = M.figure_axis_xy()
    plt.loglog(k, G_gFT_smth / k)
    plt.title('displacement power Spectra', loc='left')

    # %% cut-off wavenumber per stancil (displacement power spectrum)
    k_lim_list = list()
    k_end_previous = np.nan
    k = G_gFT_smth.k.data
    n_not_converged = 0

    for x in G_gFT_smth.x.data:
        print(x)
        k_end, pw_fit = define_noise_wavenumber_piecewise(G_gFT_smth.sel(x=x) / k, prm['n_breakpoints'],
                                                          prm['n_breakpoints_max'], plot_flag=False)
        if not pw_fit.get_results()['converged']:
            n_not_converged += 1

        k_save = k_end_previous if k_end == k[0] else k_end
        k_end_previous = k_save
        k_lim_list.append(k_save)
        print('--------------------------')

    # %% write k limits to datasets
    font_for_pres()
    G_gFT_smth.coords['k_lim'] = ('x', k_lim_list)
    G_gFT_smth.k_lim.plot()
    k_lim_smth = G_gFT_smth.k_lim.rolling(x=prm['smooth_k_lim_x'], center=True, min_periods=1).mean()
    k_lim_smth.plot(c='r')

    plt.title('k_c filter', loc='left')
    save_fig(F, plot_path, 'B06_k_cutoff', pdf=False)

    G_gFT_smth['k_lim'] = k_lim_smth
    G_gFT_wmean.coords['k_lim'] = k_lim_smth

    coverage = G_gFT_smth["N_per_stancil_fraction"]
    run.info(n_x=int(G_gFT_smth.x.size),
             n_pw_not_converged=int(n_not_converged),
             k_lim_min=float(np.nanmin(k_lim_smth.data)), k_lim_max=float(np.nanmax(k_lim_smth.data)),
             k_lim_n_nan=int(np.isnan(k_lim_smth.data).sum()),
             coverage_mean=float(coverage.mean()),
             coverage_frac_above_min=float((coverage >= prm['coverage_min']).mean()))

    # %% overview figure
    font_for_print()

    fn = copy.copy(lstrings)
    F = M.figure_axis_xy(fig_sizes['two_column'][0], fig_sizes['two_column'][0] * 0.9, container=True, view_scale=1)

    plt.suptitle('Cut-off Frequency for Displacement Spectral\n' + io.ID_to_str(ID_name), y=0.97)
    gs = GridSpec(8, 3, wspace=0.1, hspace=1.5)

    k_lims = G_gFT_wmean.k_lim
    xlims = G_gFT_wmean.k[0], G_gFT_wmean.k[-1]

    for pos, k, pflag in zip([gs[0:2, 0], gs[0:2, 1], gs[0:2, 2]], high_beams, [True, False, False]):
        ax0 = F.fig.add_subplot(pos)
        Gplot = Gk.sel(beam=k).isel(x=slice(0, -1)).gFT_PSD_data.squeeze().rolling(k=prm['plot_smooth_k'], x=prm['plot_smooth_x'], min_periods=1, center=True).mean()

        Gplot = Gplot.where(Gplot["N_per_stancil"] / Gplot["Lpoints"] >= prm['coverage_min'])

        alpha_range = iter(np.linspace(1, 0, Gplot.x.data.size))
        for x in Gplot.x.data:
            ialpha = next(alpha_range)
            plt.loglog(Gplot.k, Gplot.sel(x=x) / Gplot.k, linewidth=0.5, color=col.rels[k], alpha=ialpha)
            ax0.axvline(k_lims.sel(x=x), linewidth=0.4, color='black', zorder=0, alpha=ialpha)

        plt.title(next(fn) + k, color=col_dict[k], loc='left')
        plt.xlim(xlims)

        if pflag:
            ax0.tick_params(labelbottom=False, bottom=True)
            plt.ylabel("Power (m$^2$/k')")
            plt.legend()
        else:
            ax0.tick_params(labelbottom=False, bottom=True, labelleft=False)

    for pos, k, pflag in zip([gs[2:4, 0], gs[2:4, 1], gs[2:4, 2]], low_beams, [True, False, False]):
        ax0 = F.fig.add_subplot(pos)
        Gplot = Gk.sel(beam=k).isel(x=slice(0, -1)).gFT_PSD_data.squeeze().rolling(k=prm['plot_smooth_k'], x=prm['plot_smooth_x'], min_periods=1, center=True).mean()

        Gplot = Gplot.where(Gplot["N_per_stancil"] / Gplot["Lpoints"] >= prm['coverage_min'])

        alpha_range = iter(np.linspace(1, 0, Gplot.x.data.size))
        for x in Gplot.x.data:
            ialpha = next(alpha_range)
            plt.loglog(Gplot.k, Gplot.sel(x=x) / Gplot.k, linewidth=0.5, color=col.rels[k], alpha=ialpha)
            ax0.axvline(k_lims.sel(x=x), linewidth=0.4, color='black', zorder=0, alpha=ialpha)

        plt.title(next(fn) + k, color=col_dict[k], loc='left')
        plt.xlim(xlims)
        plt.xlabel("observed wavenumber k' ")

        if pflag:
            ax0.tick_params(bottom=True)
            plt.ylabel("Power (m$^2$/k')")
            plt.legend()
        else:
            ax0.tick_params(bottom=True, labelleft=False)

    save_fig(F, plot_path, 'B06_atten_ov_simple', close=False)

    # %% add mean displacement spectrogram and coverage panels
    pos = gs[5:, 0:2]
    ax0 = F.fig.add_subplot(pos)

    lat_str = str(np.round(Gx.isel(x=0).lat.mean().data, 2)) + ' to ' + str(np.round(Gx.isel(x=-1).lat.mean().data, 2))
    plt.title(next(fn) + 'Mean Displacement Spectra\n(lat=' + lat_str + ')', loc='left')

    dd = (10 * np.log((G_gFT_smth / G_gFT_smth.k).isel(x=slice(0, -1))))
    dd = dd.where(~np.isinf(dd), np.nan)

    # filter out segments with less then coverage_min of data points
    dd = dd.where(G_gFT_smth["N_per_stancil_fraction"] >= prm['coverage_min'])

    dd_lims = np.round(dd.quantile(0.01).data * 0.95, 0), np.round(dd.quantile(0.95).data * 1.05, 0)
    plt.pcolor(dd.x / 1e3, dd.k, dd, vmin=dd_lims[0], vmax=dd_lims[-1], cmap=col.white_base_blgror)
    cb = plt.colorbar(orientation='vertical')

    cb.set_label('Power (m$^2$/k)')
    plt.plot(G_gFT_smth.isel(x=slice(0, -1)).x / 1e3, G_gFT_smth.isel(x=slice(0, -1)).k_lim, color=col.black, linewidth=1)
    plt.ylabel('wavenumber k')
    plt.xlabel('X (km)')

    pos = gs[6:, -1]
    ax9 = F.fig.add_subplot(pos)

    plt.title('Data Coverage (%)', loc='left')
    plt.plot(G_gFT_smth.x / 1e3, G_gFT_smth["N_per_stancil_fraction"] * 100, linewidth=0.8, color='black')
    ax9.spines['left'].set_visible(False)
    ax9.spines['right'].set_visible(True)
    ax9.tick_params(labelright=True, right=True, labelleft=False, left=False)
    ax9.axhline(prm['coverage_min'] * 100, linewidth=0.8, linestyle='--', color='black')
    plt.xlabel('X (km)')

    save_fig(F, plot_path, 'B06_atten_ov')

    # %% reconstruct low-passed height model per stancil
    G_height_model = dict()
    dx_eta = Gx.eta.diff('eta').mean().data
    # NOTE: the original selects the binned table of 'gt2l' for every beam bb (k was hard-coded);
    # it is only used to decide whether the stancil has data. Kept as in the original.
    k = 'gt2l'
    for bb in Gx.beam.data:
        G_height_model_temp = dict()
        for i in np.arange(Gx.x.size):

            Gx_1 = Gx.isel(x=i).sel(beam=bb)
            Gk_1 = Gk.isel(x=i).sel(beam=bb)
            # NOTE: uses k_lim of x=0 for all stancils (as in the original)
            k_thresh = G_gFT_smth.k_lim.isel(x=0).data

            dist_stencil = Gx_1.eta + Gx_1.x
            dist_stencil_lims = dist_stencil[0].data, dist_stencil[-1].data

            T3_sel = B3[k].loc[((B3[k]['dist'] >= dist_stencil_lims[0]) & (B3[k]['dist'] <= dist_stencil_lims[1]))]

            if T3_sel.shape[0] != 0:
                height_model, dist_nanmask = reconstruct_displacement(Gx_1, Gk_1, k_thresh=k_thresh,
                                                                      sigma_g=prm['tanh_sigma_g'], dx=dx_eta)
                G_height_model_temp[str(i) + bb] = xr.DataArray(height_model, coords=Gx_1.coords, dims=Gx_1.dims, name='height_model')
            else:
                G_height_model_temp[str(i) + bb] = xr.DataArray(Gx_1.y_model.data, coords=Gx_1.coords, dims=Gx_1.dims, name='height_model')

        G_height_model[bb] = xr.concat(G_height_model_temp.values(), dim='x', join='outer', coords='different', compat='equals').T

    Gx['height_model'] = xr.concat(G_height_model.values(), dim='beam', join='outer', coords='different', compat='equals').transpose('eta', 'beam', 'x')

    # %% residual heights on the binned table
    B3_v2 = dict()
    for bb in Gx.beam.data:
        print(bb)
        Gx_k = Gx.sel(beam=bb)
        Gh = Gx['height_model'].sel(beam=bb).T
        Gh_err = Gx_k['model_error_x'].T
        Gnans = np.isnan(Gx_k.y_model)

        concented_heights = Gh.data.reshape(Gh.data.size)
        concented_err = Gh_err.data.reshape(Gh.data.size)
        concented_nans = Gnans.data.reshape(Gnans.data.size)
        concented_x = (Gh.x + Gh.eta).data.reshape(Gh.data.size)

        dx = Gh.eta.diff('eta')[0].data
        continous_x_grid = np.arange(concented_x.min(), concented_x.max(), dx)
        continous_height_model = np.interp(continous_x_grid, concented_x, concented_heights)
        concented_err = np.interp(continous_x_grid, concented_x, concented_err)
        continous_nans = np.interp(continous_x_grid, concented_x, concented_nans) == 1

        T3 = B3[bb]
        T3 = T3.sort_values('x')
        T3 = T3.sort_values('dist')

        T3['heights_c_model'] = np.interp(T3['dist'], continous_x_grid, continous_height_model)
        T3['heights_c_model_err'] = np.interp(T3['dist'], continous_x_grid, concented_err)
        T3['heights_c_residual'] = T3['heights_c_weighted_mean'] - T3['heights_c_model']

        B3_v2[bb] = T3

    # %% wave incident direction from the B05 angle pdf
    G_angle = xr.open_dataset(angle_file)

    Ga_abs = (G_angle.weighted_angle_PDF_smth.isel(angle=G_angle.angle > 0).data + G_angle.weighted_angle_PDF_smth.isel(angle=G_angle.angle < 0).data[:, ::-1]) / 2
    # Ga_abs is (x, angle); use explicit dims, Dataset.dims ordering is not guaranteed in recent xarray
    Ga_abs = xr.DataArray(data=Ga_abs, dims=('x', 'angle'), coords=G_angle.isel(angle=G_angle.angle > 0).coords)

    Ga_abs_front = Ga_abs.isel(x=slice(0, prm['angle_front_x']))
    Ga_best = ((Ga_abs_front * Ga_abs_front.N_data).sum('x') / Ga_abs_front.N_data.sum('x'))

    theta = Ga_best.angle[np.argmax(Ga_best.data)].data     # flat index; DataArray.argmax() without dim is deprecated
    theta_flag = True
    run.info(theta_deg=float(theta * 180 / np.pi), theta_applied=bool(theta_flag))

    font_for_print()
    F = M.figure_axis_xy(3, 5, view_scale=0.7)

    plt.subplot(2, 1, 1)
    plt.pcolor(Ga_abs)
    plt.xlabel('abs angle')
    plt.ylabel('x')

    ax = plt.subplot(2, 1, 2)
    Ga_best.plot()
    plt.title('angle front ' + str(theta * 180 / np.pi), loc='left')
    ax.axvline(theta, color='red')
    save_fig(F, plot_path, 'B06_angle_def', pdf=False)

    # %% corrected wavenumber and distance axes
    lam_p = 2 * np.pi / Gk.k
    lam = lam_p * np.cos(theta)

    if theta_flag:
        k_corrected = 2 * np.pi / lam
        x_corrected = Gk.x * np.cos(theta)
    else:
        k_corrected = 2 * np.pi / lam * np.nan
        x_corrected = Gk.x * np.cos(theta) * np.nan

    # %% spectral save
    G5 = G_gFT_wmean.expand_dims(dim='beam', axis=1)
    G5.coords['beam'] = ['weighted_mean']
    G5 = G5.assign_coords(N_photons=G5.N_photons)
    G5['N_photons'] = G5['N_photons'].expand_dims('beam')
    G5['N_per_stancil_fraction'] = G5['N_per_stancil_fraction'].expand_dims('beam')

    Gk_v2 = xr.merge([Gk, G5], join='outer', compat='no_conflicts')   # explicit current defaults, see above

    Gk_v2 = Gk_v2.assign_coords(x_corrected=("x", x_corrected.data)).assign_coords(k_corrected=("k", k_corrected.data))

    Gk_v2.attrs['best_guess_incident_angle'] = theta

    Gk_v2.to_netcdf(save_path + 'B06_' + ID_name + '_gFT_k_corrected.nc')

    # %% save real space data
    Gx.to_netcdf(save_path + 'B06_' + ID_name + '_gFT_x_corrected.nc')
    save_pandas_table_overwrite(B3_v2, 'B06_' + ID_name + '_binned_resid', save_path)   # regridded heights + model + residual
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
