# %%
"""
B04: MCMC inversion of the wave propagation angle per beam group and x position.

    python stages/B04_angle.py <ID> <batch_key>

Reads  work/<batch>/B02_spectra/B02_<ID>_gFT_x.nc, _gFT_k.nc                 (B02)
       work/<batch>/A02_prior/A02_<ID>.h5   table 'priors_hindcast'          (A02)
Writes work/<batch>/B04_angle/B04_<ID>_marginals.nc, B04_<ID>_res_table.h5
       status/B04/<ID>.json
       figures B04_prior_angle, B04_data_avail, B04_marginal_distributions

Port of analysis_db/B04_angle.py: same algorithm, parameters from params/<v>.toml [B04],
exit() -> SkipTrack, serial loop over wavenumbers (the ProcessPoolExecutor variant was dead code).
The B01 binned table was loaded but never used and is not read any more.
"""
import sys
import time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline_config import (mconfig, np, pd, xr, plt, GridSpec, M, MT, col, paths_for, save_fig,
                             font_for_print, font_for_pres, cli_args)
from pipeline_status import StageRun, SkipTrack, require_upstream
from pipeline_params import load_params

from ICEsat2_SI_tools import angle_optimizer

STAGE = 'B04'

col.colormaps2(21)
col_dict = col.rels


def define_wavenumber_weights_tot_var(dd, m=3, variance_frac=0.33, k_upper_lim=None, verbose=False):
    """
    return peaks of a power spectrum dd that in the format such that they can be used as weights for the frequencies based fitting

    inputs:
    dd             xarray with PSD as data amd coordindate wavenumber k
    m               running mean half-width in gridpoints
    variance_frac  (0 to 1) How much variance should be explained by the returned peaks
    verbose        if true it plots some stuff

    return:
    mask           size of dd. where True the data is identified as having significant amplitude
    k              wanumbers where mask is true
    dd_rm          smoothed version of dd
    positions      postions where of significant data in array
    """
    if len(dd.shape) == 2:
        dd_use = dd.mean('beam')

    if m is None:
        dd_rm = dd_use.data
    else:
        dd_rm = M.runningmean(dd_use, m, tailcopy=True)

    k = dd_use.k[~np.isnan(dd_rm)].data
    dd_rm = dd_rm[~np.isnan(dd_rm)]

    orders = dd_rm.argsort()[::-1]
    var_mask = dd_rm[orders].cumsum() / dd_rm.sum() < variance_frac
    pos_cumsum = orders[var_mask]
    mask = var_mask[orders.argsort()]
    if k_upper_lim is not None:
        mask = (k < k_upper_lim) & mask

    if verbose:
        plt.plot(dd.k, dd, '-', color=col_dict[str(dd.beam[0].data)], markersize=20, alpha=0.6)
        plt.plot(k, dd_rm, '-k', markersize=20)
        plt.plot(k[mask], dd_rm[mask], '.r', markersize=10, zorder=12)
        if k_upper_lim is not None:
            plt.gca().axvline(k_upper_lim, color='black')

    return mask, k, dd_rm, pos_cumsum


def define_wavenumber_weights_threshold(dd, m=3, Nstd=2, verbose=False):

    if m is None:
        dd_rm = dd
    else:
        dd_rm = M.runningmean(dd, m, tailcopy=True)

    k = dd.k[~np.isnan(dd_rm)]
    dd_rm = dd_rm[~np.isnan(dd_rm)]

    treshold = np.nanmean(dd_rm) + np.nanstd(dd_rm) * Nstd
    mask = dd_rm > treshold

    if verbose:
        plt.plot(dd.k, dd, '-k', markersize=20)
        plt.plot(k, dd_rm, '-b', markersize=20)
        plt.plot(k[mask], dd_rm[mask], '.r', markersize=10, zorder=12)

    return mask, k, dd_rm, np.arange(0, mask.size)[mask]


def plot_instance(z_model, fargs, key, SM, non_dim=False, title_str=None, brute=False, optimze=False, sample=False, view_scale=0.3):

    x_concat, y_concat, z_concat = fargs

    F = M.figure_axis_xy(5, 6, view_scale=view_scale, container=True)
    plt.suptitle(title_str)
    gs = GridSpec(4, 3, wspace=0.4, hspace=1.2)
    F.gs = gs

    import itertools
    col_list = itertools.cycle([col.cascade2, col.rascade2, col.cascade1, col.rascade1, col.cascade3, col.rascade3])

    beam_list = list(set(y_concat))
    for y_pos, pos in zip(beam_list, [gs[0, :], gs[1, :]]):

        F.ax2 = F.fig.add_subplot(pos)
        plt.title(str(y_pos))
        plt.plot(x_concat[y_concat == y_pos], z_concat[y_concat == y_pos], c=col.gray, linewidth=1)
        plt.plot(x_concat[y_concat == y_pos], z_model[y_concat == y_pos], '-', c=next(col_list))
        plt.xlim(x_concat[y_concat == y_pos][0], x_concat[y_concat == y_pos][-1])

    plt.xlabel('meter')
    F.ax3 = F.fig.add_subplot(gs[2:, 0:-1])
    if brute is True:
        plt.title('Brute-force costs', loc='left')
        SM.plot_brute(marker='.', color='blue', markersize=15, label='Brute', zorder=10)
    if optimze is True:
        SM.plot_optimze(color='r', markersize=10, zorder=12, label='Dual Annealing')
    if sample is True:
        SM.plot_sample(markersize=2, linewidth=0.8, alpha=0.2, color='black', zorder=8)

    F.ax4 = F.fig.add_subplot(gs[2:, -1])
    return F


def run_stage(ID, batch_key, prm, run):
    track_name = ID
    P = paths_for(batch_key, ID)
    require_upstream(run, ['B02', 'A02'])

    beam_groups = mconfig['beams']['groups']

    save_path = P.stage_dir('B04_angle')
    save_name = 'B04_' + track_name
    plot_path = P.plot_track

    # %% load B02 spectra
    load_path = P.stage_dir('B02_spectra')
    Gx = xr.load_dataset(load_path + 'B02_' + track_name + '_gFT_x.nc')
    Gk = xr.load_dataset(load_path + 'B02_' + track_name + '_gFT_k.nc')

    # %% load prior information
    load_path = P.stage_dir('A02_prior')
    try:
        Prior = MT.load_pandas_table_dict('/A02_' + track_name, load_path)['priors_hindcast']
    except Exception as e:
        raise SkipTrack('prior table not found', error=str(e)[:200])

    if np.isnan(Prior['mean']['dir']):
        raise SkipTrack('prior direction is nan')

    # cast to float numpy arrays: the 'mean' column is object dtype and pandas >= 3 no longer
    # falls back to positional indexing when a Series is indexed with argsort() positions
    Pperiod = Prior.loc[['ptp0', 'ptp1', 'ptp2', 'ptp3', 'ptp4', 'ptp5']]['mean'].astype('float').to_numpy()
    Pdir = Prior.loc[['pdp0', 'pdp1', 'pdp2', 'pdp3', 'pdp4', 'pdp5']]['mean'].astype('float').to_numpy()
    Pspread = Prior.loc[['pspr0', 'pspr1', 'pspr2', 'pspr3', 'pspr4', 'pspr5']]['mean'].astype('float').to_numpy()

    # partitions that are absent in WW3 have nan period/direction but spread 0.0 (not nan); the old
    # code masked on the spread only and the nan then propagated through np.interp into the prior
    valid = ~(np.isnan(Pperiod) | np.isnan(Pdir) | np.isnan(Pspread))
    Pperiod, Pdir, Pspread = Pperiod[valid], Pdir[valid], Pspread[valid]
    run.info(n_prior_partitions=int(valid.sum()))

    # this is a hack since the current data does not have a spread
    Pspread[Pspread == 0] = prm['prior_spread_fill_deg']

    # reset dirs:
    Pdir[Pdir > 180] = Pdir[Pdir > 180] - 360
    Pdir[Pdir < -180] = Pdir[Pdir < -180] + 360

    # reorder dirs
    dir_best = [0]
    for dir in Pdir:
        ip = np.argmin([abs(dir_best[-1] - dir), abs(dir_best[-1] - (dir - 360)), abs(dir_best[-1] - (dir + 360))])
        new_dir = np.array([dir, (dir - 360), (dir + 360)])[ip]
        dir_best.append(new_dir)
    dir_best = np.array(dir_best[1:])

    # %% interpolate prior on the wavenumber axis
    n_smooth = prm['prior_smooth_points']
    if len(Pperiod) == 0:
        print('constant peak wave number')
        kk = Gk.k
        Pwavenumber = kk * 0 + (2 * np.pi / (1 / Prior.loc['fp']['mean'])) ** 2 / 9.81
        dir_best = kk * 0 + Prior.loc['dp']['mean']
        dir_interp_smth = dir_interp = kk * 0 + Prior.loc['dp']['mean']
        spread_smth = spread_interp = kk * 0 + Prior.loc['spr']['mean']

    else:
        Pwavenumber = (2 * np.pi / Pperiod) ** 2 / 9.81
        kk = Gk.k
        dir_interp = np.interp(kk, Pwavenumber[Pwavenumber.argsort()], dir_best[Pwavenumber.argsort()])
        dir_interp_smth = M.runningmean(dir_interp, n_smooth, tailcopy=True)
        dir_interp_smth[-1] = dir_interp_smth[-2]

        spread_interp = np.interp(kk, Pwavenumber[Pwavenumber.argsort()], Pspread[Pwavenumber.argsort()].astype('float'))
        spread_smth = M.runningmean(spread_interp, n_smooth, tailcopy=True)
        spread_smth[-1] = spread_smth[-2]

    font_for_pres()

    F = M.figure_axis_xy(5, 4.5, view_scale=0.5)
    plt.subplot(2, 1, 1)
    plt.title('Prior angle smoothed\n' + track_name, loc='left')

    plt.plot(Pwavenumber, dir_best, '.r', markersize=8)
    plt.plot(kk, dir_interp, '-', color='red', linewidth=0.8, zorder=11)
    plt.plot(kk, dir_interp_smth, color=col.green1)

    plt.fill_between(kk, dir_interp_smth - spread_smth, dir_interp_smth + spread_smth, zorder=1, color=col.green1, alpha=0.2)
    plt.ylabel('Angle (deg)')

    ax2 = plt.subplot(2, 1, 2)
    plt.title('Prior angle adjusted ', loc='left')

    # adjust angle def:
    dir_interp_smth[dir_interp_smth > 180] = dir_interp_smth[dir_interp_smth > 180] - 360
    dir_interp_smth[dir_interp_smth < -180] = dir_interp_smth[dir_interp_smth < -180] + 360

    plt.fill_between(kk, dir_interp_smth - spread_smth, dir_interp_smth + spread_smth, zorder=1, color=col.green1, alpha=0.2)
    plt.plot(kk, dir_interp_smth, '.', markersize=1, color=col.green1)

    ax2.axhline(85, color='gray', linewidth=2)
    ax2.axhline(-85, color='gray', linewidth=2)

    plt.ylabel('Angle (deg)')
    plt.xlabel(r'wavenumber ($2 \pi/\lambda$)')

    save_fig(F, plot_path, 'B04_prior_angle')

    # save
    dir_interp_smth = xr.DataArray(data=dir_interp_smth * np.pi / 180, dims='k', coords={'k': kk}, name='Prior_direction')
    spread_smth = xr.DataArray(data=spread_smth * np.pi / 180, dims='k', coords={'k': kk}, name='Prior_spread')
    Prior_smth = xr.merge([dir_interp_smth, spread_smth])

    # %% prior angle gate
    prior_angle = Prior_smth.Prior_direction * 180 / np.pi
    run.info(prior_dir_deg=float(prior_angle.median()),
             prior_spread_deg=float((Prior_smth.Prior_spread * 180 / np.pi).median()))
    if (abs(prior_angle) > prm['angle_gate_deg']).all():
        print('Prior angle is ', prior_angle.min().data, prior_angle.max().data, '. quit.')
        raise SkipTrack('prior angle outside gate',
                        angle_min=float(prior_angle.min().data), angle_max=float(prior_angle.max().data),
                        angle_median=float(prior_angle.median()))

    # %% define paramater range
    params_dict = {'alpha': [prm['alpha_range'][0] * np.pi / 2, prm['alpha_range'][1] * np.pi / 2, prm['alpha_walkers']],
                   'phase': [0, 2 * np.pi, prm['phase_walkers']]}

    alpha_dx = prm['alpha_dx']
    max_wavenumbers = prm['max_wavenumbers']

    sample_flag = True
    optimize_flag = False
    brute_flag = False

    plot_flag = False

    N_sample_chain = prm['N_sample_chain']
    N_sample_chain_burn = prm['N_sample_chain_burn']

    max_x_pos = prm['max_x_pos']
    x_pos_jump = prm['x_pos_jump']

    def make_fake_data(xi, group):
        ki = Gk.k[0:2]

        bins = np.arange(params_dict['alpha'][0], params_dict['alpha'][1] + alpha_dx, alpha_dx)
        bins_pos = (bins[0:-1] + np.diff(bins) / 2)
        marginal_stack = xr.DataArray(np.nan * np.vstack([bins_pos, bins_pos]).T, dims=('angle', 'k'), coords={'angle': bins_pos, 'k': ki.data})

        group_name = str('group' + group[0].split('gt')[1].split('l')[0])
        marginal_stack.coords['beam_group'] = group_name
        marginal_stack.coords['x'] = xi
        marginal_stack.name = 'marginals'
        marginal_stack.expand_dims(dim='x', axis=2).expand_dims(dim='beam_group', axis=3)
        return marginal_stack

    # %% isolate x positions with data
    data_mask = Gk.gFT_PSD_data.mean('k')
    data_mask.coords['beam_group'] = ('beam', ['beam_group' + g[2] for g in data_mask.beam.data])
    data_mask_group = data_mask.groupby('beam_group').mean(skipna=False)
    # stancils that are actually usable: at least min_groups_with_data beam pairs have a finite
    # spectrum there (the old code took any stancil where *some* pair had data and then found no
    # data for the other pairs -> all-dummy tracks on sparse / ice-covered tracks)
    n_groups_with_data = (~np.isnan(data_mask_group)).sum('beam_group')
    data_sel_mask = n_groups_with_data >= prm.get('min_groups_with_data', 1)
    run.info(n_x_any_group=int((n_groups_with_data > 0).sum()), n_x_selected=int(data_sel_mask.sum()))

    x_list = data_sel_mask.x[data_sel_mask]  # iterate over these x posistions
    x_list_flag = ~np.isnan(data_mask_group.sel(x=x_list))  # flag that is False if there is no data

    # limit number of x coordinates
    x_list = x_list[::x_pos_jump]
    if len(x_list) > max_x_pos:
        x_list = x_list[0:max_x_pos]
    x_list_flag = x_list_flag.sel(x=x_list)
    run.info(n_x_pos=int(x_list.size), n_groups=len(beam_groups))

    # plot
    font_for_print()
    F = M.figure_axis_xy(5.5, 3, view_scale=0.8)
    plt.suptitle(track_name)
    ax1 = plt.subplot(2, 1, 1)
    plt.title('Data in Beam', loc='left')
    plt.pcolormesh(data_mask.x / 1e3, data_mask.beam, data_mask, cmap=plt.cm.OrRd)
    for i in np.arange(1.5, 6, 2):
        ax1.axhline(i, color='black', linewidth=0.5)
    plt.xlabel('Distance from Ice Edge')

    ax2 = plt.subplot(2, 1, 2)
    plt.title('Data in Group', loc='left')
    plt.pcolormesh(data_mask.x / 1e3, data_mask_group.beam_group, data_mask_group, cmap=plt.cm.OrRd)

    for i in np.arange(0.5, 3, 1):
        ax2.axhline(i, color='black', linewidth=0.5)

    plt.plot(x_list / 1e3, x_list * 0 + 0, '.', markersize=2, color=col.cascade1)
    plt.plot(x_list / 1e3, x_list * 0 + 1, '.', markersize=2, color=col.cascade1)
    plt.plot(x_list / 1e3, x_list * 0 + 2, '.', markersize=2, color=col.cascade1)

    plt.xlabel('Distance from Ice Edge')

    save_fig(F, plot_path, 'B04_data_avail')

    # %% MCMC loop over beam groups and x positions
    Marginals = dict()
    L_collect = dict()
    n_wavenumbers_used = dict()
    runtime_s = dict()
    n_dummy = 0

    group_number = np.arange(len(beam_groups))
    ggg, xxx = np.meshgrid(group_number, x_list.data)

    for gi in zip(ggg.flatten(), xxx.flatten()):
        print(gi)

        group, xi = beam_groups[gi[0]], gi[1]
        ikey = str(xi) + '_' + '_'.join(group)

        if bool(x_list_flag.sel(x=xi).isel(beam_group=gi[0]).data) is False:
            print('no data, fill with dummy')
            Marginals[ikey] = make_fake_data(xi, group)
            n_dummy += 1
            continue

        t0 = time.time()
        GGx = Gx.sel(beam=group).sel(x=xi)
        GGk = Gk.sel(beam=group).sel(x=xi)

        # define data
        # normalize data
        key = 'y_data'
        amp_Z = (GGx[key] - GGx[key].mean(['eta'])) / GGx[key].std(['eta'])
        N = amp_Z.shape[0]

        # define x,y positions
        eta_2d = GGx.eta + GGx.x_coord - GGx.x_coord.mean()
        nu_2d = GGx.eta * 0 + GGx.y_coord - GGx.y_coord.mean()

        # repack as np arrays
        x_concat = eta_2d.data.T.flatten()
        y_concat = nu_2d.data.T.flatten()
        z_concat = amp_Z.data.flatten()

        # x_coord/y_coord can be nan for a beam without data at this stancil (B02 sets them nan);
        # the old code only masked on z and then fed nan positions to the likelihood
        finite = ~np.isnan(z_concat) & ~np.isnan(x_concat) & ~np.isnan(y_concat)
        x_concat, y_concat, z_concat = x_concat[finite], y_concat[finite], z_concat[finite]
        N_data = x_concat.size

        mean_dist = (nu_2d.isel(beam=0) - nu_2d.isel(beam=1)).mean().data
        k_upper_lim = 2 * np.pi / (mean_dist * 1)
        if N_data == 0 or not np.isfinite(k_upper_lim):
            print('no finite data/positions in this beam pair, fill with dummy')
            Marginals[ikey] = make_fake_data(xi, group)
            n_dummy += 1
            continue

        print('k_upper_lim ', k_upper_lim)

        # variance method
        amp_data = np.sqrt(GGk.gFT_cos_coeff**2 + GGk.gFT_sin_coeff**2)
        mask, k, weights, positions = define_wavenumber_weights_tot_var(amp_data, m=1, k_upper_lim=k_upper_lim,
                                                                        variance_frac=prm['variance_frac'], verbose=False)

        if (len(k[mask]) == 0):
            print('no good k found, fill with dummy')
            Marginals[ikey] = make_fake_data(xi, group)
            n_dummy += 1
            continue

        # prepare loop: init object and test
        SM = angle_optimizer.sample_with_mcmc(params_dict)
        SM.set_objective_func(angle_optimizer.objective_func)

        SM.fitting_args = fitting_args = (x_concat, y_concat, z_concat)

        # test:
        k_prime_max = prm['k_prime_max_test']  # chose a test wavenumber
        amp_Z = 1
        prior_sel = {'alpha': (Prior_smth.sel(k=k_prime_max, method='nearest').Prior_direction.data,
                               Prior_smth.sel(k=k_prime_max, method='nearest').Prior_spread.data)}
        SM.fitting_kargs = fitting_kargs = {'prior': prior_sel, 'prior_weight': 3}
        # test if it works
        SM.params.add('K_prime', k_prime_max, vary=False, min=k_prime_max * 0.5, max=k_prime_max * 1.5)
        SM.params.add('K_amp', amp_Z, vary=False, min=amp_Z * .0, max=amp_Z * 5)
        try:
            SM.test_objective_func()
        except Exception as e:
            raise ValueError('Objective function test fails') from e

        def get_instance(k_pair):

            k_prime_max, Z_max = k_pair

            prior_sel = {'alpha': (Prior_smth.sel(k=k_prime_max, method='nearest').Prior_direction.data,
                                   Prior_smth.sel(k=k_prime_max, method='nearest').Prior_spread.data)}

            SM.fitting_kargs = fitting_kargs = {'prior': prior_sel, 'prior_weight': prm['prior_weight']}

            amp_Z = 1

            SM.params.add('K_prime', k_prime_max, vary=False, min=k_prime_max * 0.5, max=k_prime_max * 1.5)
            SM.params.add('K_amp', amp_Z, vary=False, min=amp_Z * .0, max=amp_Z * 5)
            L_sample_i = None
            L_optimize_i = None
            L_brute_i = None
            if sample_flag:
                SM.sample(verbose=False, steps=N_sample_chain, progress=False, workers=None)
                L_sample_i = list(SM.fitter.params.valuesdict().values())  # mcmc results

            elif optimize_flag:
                SM.optimize(verbose=False)
                L_optimize_i = list(SM.fitter_optimize.params.valuesdict().values())

            elif brute_flag:
                SM.brute(verbose=False)
                L_brute_i = list(SM.fitter_brute.params.valuesdict().values())
            else:
                raise ValueError('non of sample_flag,optimize_flag, or brute_flag  are True')

            y_hist, bins, bins_pos = SM.get_marginal_dist('alpha', alpha_dx, burn=N_sample_chain_burn, plot_flag=False)
            fitter = SM.fitter  # MCMC results
            z_model = SM.objective_func(fitter.params, *fitting_args, test_flag=True)
            cost = (fitter.residual**2).sum() / (z_concat**2).sum()

            if plot_flag:

                F = plot_instance(z_model, fitting_args, 'y_data_normed', SM, brute=brute_flag, optimze=optimize_flag,
                                  sample=sample_flag, title_str='k=' + str(np.round(k_prime_max, 4)), view_scale=0.6)

                if (fitting_kargs['prior'] is not None):
                    F.ax3.axhline(prior_sel['alpha'][0], color='green', linewidth=2, label='Prior')
                    F.ax3.axhline(prior_sel['alpha'][0] - prior_sel['alpha'][1], color='green', linewidth=0.7)
                    F.ax3.axhline(prior_sel['alpha'][0] + prior_sel['alpha'][1], color='green', linewidth=0.7)

                F.ax3.axhline(fitter.params['alpha'].min, color='gray', linewidth=2, alpha=0.6)
                F.ax3.axhline(fitter.params['alpha'].max, color='gray', linewidth=2, alpha=0.6)

                plt.sca(F.ax3)
                plt.legend()
                plt.xlabel('Phase')
                plt.ylabel('Angle')
                plt.xlim(0, np.pi * 2)

                plt.sca(F.ax4)
                plt.xlabel('Density')
                plt.stairs(y_hist, bins, orientation='horizontal', color='k')

                F.ax4.axhline(fitter.params['alpha'].min, color='gray', linewidth=2, alpha=0.6)
                F.ax4.axhline(fitter.params['alpha'].max, color='gray', linewidth=2, alpha=0.6)

                F.ax3.set_ylim(min(-np.pi / 2, prior_sel['alpha'][0] - 0.2), max(np.pi / 2, prior_sel['alpha'][0] + 0.2))
                F.ax4.set_ylim(min(-np.pi / 2, prior_sel['alpha'][0] - 0.2), max(np.pi / 2, prior_sel['alpha'][0] + 0.2))

                save_fig(F, plot_path, track_name + '_fit_k' + str(k_prime_max))

            marginal_stack_i = xr.DataArray(y_hist, dims=('angle'), coords={'angle': bins_pos})
            marginal_stack_i.coords['k'] = np.array(k_prime_max)

            rdict = {
                'marginal_stack_i': marginal_stack_i,
                'L_sample_i': L_sample_i,
                'L_optimize_i': L_optimize_i,
                'L_brute_i': L_brute_i,
                'cost': cost
            }
            return k_prime_max, rdict

        k_list, weight_list = k[mask], weights[mask]
        print('# of wavenumber: ', len(k_list))
        if len(k_list) > max_wavenumbers:
            print('cut wavenumber list to', max_wavenumbers)
            k_list = k_list[0:max_wavenumbers]
            weight_list = weight_list[0:max_wavenumbers]

        # drop wavenumbers without a finite prior (the likelihood cannot handle nan); keep k and
        # weight lists consistent because 'weight' is stored per k below
        prior_ok = np.array([bool(np.isfinite(Prior_smth.sel(k=kk_, method='nearest').Prior_direction.data)
                                  and np.isfinite(Prior_smth.sel(k=kk_, method='nearest').Prior_spread.data))
                             for kk_ in k_list])
        if (~prior_ok).sum():
            print(f'{int((~prior_ok).sum())} wavenumbers skipped: prior is nan')
        k_list, weight_list = k_list[prior_ok], weight_list[prior_ok]
        if len(k_list) == 0:
            print('no wavenumber with a finite prior, fill with dummy')
            Marginals[ikey] = make_fake_data(xi, group)
            n_dummy += 1
            continue

        A = dict()
        for k_pair in zip(k_list, weight_list):
            kk, I = get_instance(k_pair)
            A[kk] = I

        cost_stack = dict()
        marginal_stack = dict()
        L_sample = pd.DataFrame(index=['alpha', 'group_phase', 'K_prime', 'K_amp'])
        L_optimize = pd.DataFrame(index=['alpha', 'group_phase', 'K_prime', 'K_amp'])
        L_brute = pd.DataFrame(index=['alpha', 'group_phase', 'K_prime', 'K_amp'])

        for kk, I in A.items():
            L_sample[kk] = I['L_sample_i']
            L_optimize[kk] = I['L_optimize_i']
            L_brute[kk] = I['L_brute_i']

            marginal_stack[kk] = I['marginal_stack_i']
            cost_stack[kk] = I['cost']

        # add beam_group dimension
        marginal_stack = xr.concat(marginal_stack.values(), dim='k').sortby('k')
        L_sample = L_sample.T.sort_values('K_prime')
        L_optimize = L_optimize.T.sort_values('K_prime')
        L_brute = L_brute.T.sort_values('K_prime')

        print('done with ', group, xi / 1e3)

        # collect
        marginal_stack.name = 'marginals'
        marginal_stack = marginal_stack.to_dataset()
        marginal_stack['cost'] = (('k'), list(cost_stack.values()))
        marginal_stack['weight'] = (('k'), weight_list)

        group_name = str('group' + group[0].split('gt')[1].split('l')[0])
        marginal_stack.coords['beam_group'] = group_name
        marginal_stack.coords['x'] = xi

        Marginals[ikey] = marginal_stack.expand_dims(dim='x', axis=0).expand_dims(dim='beam_group', axis=1)
        Marginals[ikey].coords['N_data'] = (('x', 'beam_group'), np.expand_dims(np.expand_dims(N_data, 0), 1))

        L_sample['cost'] = cost_stack
        L_sample['weight'] = weight_list
        L_collect[group_name, str(int(xi))] = L_sample

        n_wavenumbers_used[group_name + '_' + str(int(xi))] = int(len(k_list))
        runtime_s[group_name + '_' + str(int(xi))] = round(time.time() - t0, 1)

    run.info(n_dummy=n_dummy, n_wavenumbers_used=n_wavenumbers_used, runtime_per_instance_s=runtime_s)

    # %% save
    if not L_collect:                      # every (group, x) instance was a dummy; the merge below would fail
        raise SkipTrack('no data in any beam group / x position', n_dummy=n_dummy)
    # explicit join/compat = the current xarray defaults (they change in a future release)
    MM = xr.merge(Marginals.values(), join='outer', compat='no_conflicts')
    MM = xr.merge([MM, Prior_smth], join='outer', compat='no_conflicts')
    MM.to_netcdf(save_path + save_name + '_marginals.nc')

    if not L_collect:
        raise SkipTrack('no data in any beam group / x position', n_dummy=n_dummy)

    LL = pd.concat(L_collect)
    MT.save_pandas_table({'L_sample': LL}, save_name + '_res_table', save_path)

    # %% plot
    font_for_print()
    F = M.figure_axis_xy(6, 5.5, view_scale=0.7, container=True)

    gs = GridSpec(4, 6, wspace=0.2, hspace=.8)

    ax0 = F.fig.add_subplot(gs[0:2, -1])
    ax0.tick_params(labelleft=False)

    klims = 0, LL['K_prime'].max() * 1.2

    for g in MM.beam_group:
        MMi = MM.sel(beam_group=g)
        plt.plot(MMi.weight.T, MMi.k, '.', color=col_dict[str(g.data)], markersize=3, linewidth=0.8)

    plt.xlabel('Power')
    plt.ylim(klims)

    ax1 = F.fig.add_subplot(gs[0:2, 0:-1])

    for g in MM.beam_group:
        if str(g.data) not in LL.index.get_level_values(0):
            continue        # group had only dummy instances
        Li = LL.loc[str(g.data)]

        angle_list = np.array(Li['alpha']) * 180 / np.pi
        kk_list = np.array(Li['K_prime'])
        weight_list_i = np.array(Li['weight'])

        plt.scatter(angle_list, kk_list, s=(weight_list_i * 8e1)**2, color=col_dict[str(g.data)], label='mode ' + str(g.data))

    dir_best[dir_best > 180] = dir_best[dir_best > 180] - 360
    plt.plot(dir_best, Pwavenumber, '.r', markersize=6)

    dir_interp[dir_interp > 180] = dir_interp[dir_interp > 180] - 360
    plt.plot(dir_interp, Gk.k, '-', color='red', linewidth=0.3, zorder=11)

    plt.fill_betweenx(Gk.k, (dir_interp_smth - spread_smth) * 180 / np.pi, (dir_interp_smth + spread_smth) * 180 / np.pi, zorder=1, color=col.green1, alpha=0.2)
    plt.plot(dir_interp_smth * 180 / np.pi, Gk.k, '.', markersize=1, color=col.green1)

    ax1.axvline(85, color='gray', linewidth=2)
    ax1.axvline(-85, color='gray', linewidth=2)

    plt.legend()
    plt.ylabel('wavenumber (deg)')
    plt.xlabel('Angle (deg)')

    plt.ylim(klims)

    # prior_sel is the prior of the last fitted (group, x) instance at the test wavenumber, as in the old script
    prior_angle_str = str(np.round((prior_sel['alpha'][0]) * 180 / np.pi))
    plt.title(track_name + '\nprior=' + prior_angle_str + 'deg', loc='left')

    plt.xlim(min([-90, np.nanmin(dir_best)]), max([np.nanmax(dir_best), 90]))

    ax3 = F.fig.add_subplot(gs[2, 0:-1])

    for g in MM.beam_group:
        MMi = MM.sel(beam_group=g)
        wegihted_margins = ((MMi.marginals * MMi.weight).sum(['x', 'k']) / MMi.weight.sum(['x', 'k']))
        plt.plot(MMi.angle * 180 / np.pi, wegihted_margins, '.', color=col_dict[str(g.data)], markersize=2, linewidth=0.8)

    plt.ylabel('Density')
    plt.title('weight margins', loc='left')
    plt.xlim(-90, 90)

    ax3 = F.fig.add_subplot(gs[-1, 0:-1])

    for g in MM.beam_group:
        MMi = MM.sel(beam_group=g)
        wegihted_margins = MMi.marginals.mean(['x', 'k'])
        plt.plot(MMi.angle * 180 / np.pi, wegihted_margins, '.', color=col_dict[str(g.data)], markersize=2, linewidth=0.8)

    plt.ylabel('Density')
    plt.xlabel('Angle (deg)')
    plt.title('unweighted margins', loc='left')
    plt.xlim(-90, 90)

    save_fig(F, plot_path, 'B04_marginal_distributions')
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
