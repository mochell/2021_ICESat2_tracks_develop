# %%
"""
B03: overview figures of the gFT / FFT slope spectra (figures only, no data products).

    python stages/B03_plot_spectra.py <ID> <batch_key>

Reads  work/<batch>/B02_spectra/B02_<ID>_gFT_k.nc, _gFT_x.nc, _FFT.nc      (B02)
Writes plots/<hemis>/<batch>/<ID>/B03_specs_coord_check, B03_specs_L<Lmeters>
       plots/<hemis>/<batch>/<ID>/B03_spectra/B03_freq_reconst_x<i>
       status/B03/<ID>.json

Port of analysis_db/B03_plot_spectra_ov.py: same figures, plotting constants from params/<v>.toml
[B03], the 'no y_data' exit -> SkipTrack, no B03_success/B03_fail.json markers.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline_config import mconfig, np, xr, plt, GridSpec, M, MT, col, paths_for, save_fig, font_for_print, cli_args
from pipeline_status import StageRun, SkipTrack, require_upstream
from pipeline_params import load_params

import generalized_FT as gFT

STAGE = 'B03'


def plot_wavenumber_spectrogram(ax, Gi, clev, title=None, plot_photon_density=True):
    if Gi.k[0] == 0:
        Gi = Gi.sel(k=Gi.k[1:])
    x_lambda = 2 * np.pi / Gi.k
    plt.pcolormesh(Gi.x / 1e3, x_lambda, Gi, cmap=plt.cm.ocean_r, vmin=clev[0], vmax=clev[-1])
    ax.set_yscale('log')

    if plot_photon_density:
        plt.plot(Gi.x / 1e3, x_lambda[-1] + (Gi.N_per_stancil / Gi.N_per_stancil.max()) * 10, c='black', linewidth=0.8, label='NAN-density')
        plt.fill_between(Gi.x / 1e3, x_lambda[-1] + (Gi.N_per_stancil / Gi.N_per_stancil.max()) * 10, 0, color='gray', alpha=0.3)
        ax.axhline(30, color='black', linewidth=0.3)

    plt.ylim(x_lambda[-1], x_lambda[0])
    plt.title(title, loc='left')


def plot_data_eta(D, offset=0, **kargs):
    eta_1 = D.eta
    y_data = D.y_model + offset
    plt.plot(eta_1, y_data, **kargs)
    return eta_1


def plot_model_eta(D, ax, offset=0, **kargs):
    eta = D.eta
    y_data = D.y_model + offset
    plt.plot(eta, y_data, **kargs)
    ax.axvline(eta[0].data, linewidth=0.1, color=kargs['color'], alpha=0.5)
    ax.axvline(eta[-1].data, linewidth=0.1, color=kargs['color'], alpha=0.5)


def run_stage(ID, batch_key, prm, run):
    P = paths_for(batch_key, ID)
    require_upstream(run, ['B02'])
    track_name = ID
    all_beams = mconfig['beams']['all_beams']
    high_beams = mconfig['beams']['high_beams']
    low_beams = mconfig['beams']['low_beams']
    col.colormaps2(21)
    col_dict = col.rels
    fltostr = MT.float_to_str

    load_path = P.stage_dir('B02_spectra', mkdir=False)
    load_file = load_path + 'B02_' + track_name
    plot_path = P.track_plot_dir()

    Gk = xr.open_dataset(load_file + '_gFT_k.nc')
    Gx = xr.open_dataset(load_file + '_gFT_x.nc')
    Gfft = xr.open_dataset(load_file + '_FFT.nc')

    # %% check paths (again)
    F = M.figure_axis_xy(9, 3, view_scale=0.5)

    plt.subplot(1, 3, 1)
    plt.title(track_name, loc='left')
    for k in all_beams:
        I = Gk.sel(beam=k)
        I2 = Gx.sel(beam=k)
        plt.plot(I['lon'], I['lat'], '.', c=col_dict[k], markersize=0.7, linewidth=0.3)
        plt.plot(I2['lon'], I2['lat'], '|', c=col_dict[k], markersize=0.7)
    plt.xlabel('lon')
    plt.ylabel('lat')

    plt.subplot(1, 3, 2)
    xscale = 1e3
    for k in all_beams:
        I = Gk.sel(beam=k)
        plt.plot(I['x_coord'] / xscale, I['y_coord'] / xscale, '.', c=col_dict[k], linewidth=0.8, markersize=0.8)
    plt.xlabel('x_coord (km)')
    plt.ylabel('y_coord (km)')

    plt.subplot(1, 3, 3)
    for k in all_beams:
        I = Gk.sel(beam=k)
        plt.plot(I['x_coord'] / xscale, (I['y_coord'] - I['y_coord'][0]), '.', c=col_dict[k], linewidth=0.8, markersize=0.8)
    plt.xlabel('x_coord (km)')
    plt.ylabel('y_coord deviation (m)')

    save_fig(F, plot_path, 'B03_specs_coord_check')

    # %% weighted means over beams
    G_gFT_wmean = (Gk['gFT_PSD_data'].where(~np.isnan(Gk['gFT_PSD_data']), 0) * Gk['N_per_stancil']).sum('beam') / Gk['N_per_stancil'].sum('beam')
    G_gFT_wmean['N_per_stancil'] = Gk['N_per_stancil'].sum('beam')

    G_fft_wmean = (Gfft.where(~np.isnan(Gfft), 0) * Gfft['N_per_stancil']).sum('beam') / Gfft['N_per_stancil'].sum('beam')
    G_fft_wmean['N_per_stancil'] = Gfft['N_per_stancil'].sum('beam')

    # %% peak wavenumber from the first x stancils
    Gmean = G_gFT_wmean.rolling(k=prm['rolling_k'], center=True).mean()
    band = prm['k_max_band']
    try:
        k_max = Gmean.k[Gmean.isel(x=slice(0, prm['k_max_first_x'][0])).mean('x').argmax().data].data
    except Exception as e:
        print('k_max from first', prm['k_max_first_x'][0], 'stancils failed, using', prm['k_max_first_x'][1], ':', repr(e))
        k_max = Gmean.k[Gmean.isel(x=slice(0, prm['k_max_first_x'][1])).mean('x').argmax().data].data
    k_max_range = k_max * band[0], k_max * 1, k_max * band[1]

    # %% spectrogram overview figure
    font_for_print()
    F = M.figure_axis_xy(6.5, 5.6, container=True, view_scale=1)
    Lmeters = Gk.L.data[0]

    plt.suptitle('gFT Slope Spectrograms\n' + track_name, y=0.98)
    gs = GridSpec(3, 3, wspace=0.2, hspace=.5)

    # define mean first for colorbar
    rk, rx = prm['rolling_k_x']
    Gplot = G_gFT_wmean.squeeze().rolling(k=rk, min_periods=1, center=True).median().rolling(x=rx, min_periods=1, center=True).median()
    dd = 10 * np.log10(Gplot)
    dd = dd.where(~np.isinf(dd), np.nan)
    clev_log = M.clevels([dd.quantile(0.01).data, dd.quantile(0.98).data * 1.2], 31) * 1

    xlims = Gmean.x[0] / 1e3, Gmean.x[-1] / 1e3

    for pos, k, pflag in zip([gs[0, 0], gs[0, 1], gs[0, 2]], high_beams, [True, False, False]):
        ax0 = F.fig.add_subplot(pos)
        Gplot = Gk.sel(beam=k).gFT_PSD_data.squeeze()
        dd2 = 10 * np.log10(Gplot)
        dd2 = dd2.where(~np.isinf(dd2), np.nan)
        plot_wavenumber_spectrogram(ax0, dd2, clev_log, title=k + ' unsmoothed', plot_photon_density=True)
        plt.xlim(xlims)
        if pflag:
            plt.ylabel('Wave length\n(meters)')
            plt.legend()

    for pos, k, pflag in zip([gs[1, 0], gs[1, 1], gs[1, 2]], low_beams, [True, False, False]):
        ax0 = F.fig.add_subplot(pos)
        Gplot = Gk.sel(beam=k).gFT_PSD_data.squeeze()
        dd2 = 10 * np.log10(Gplot)
        dd2 = dd2.where(~np.isinf(dd2), np.nan)
        plot_wavenumber_spectrogram(ax0, dd2, clev_log, title=k + ' unsmoothed', plot_photon_density=True)
        plt.xlim(xlims)
        if pflag:
            plt.ylabel('Wave length\n(meters)')
            plt.legend()

    ax0 = F.fig.add_subplot(gs[2, 0])
    plot_wavenumber_spectrogram(ax0, dd, clev_log, title='smoothed weighted mean \n10 $\\log_{10}( (m/m)^2 m )$', plot_photon_density=True)
    plt.xlim(xlims)

    ax0.axhline(2 * np.pi / k_max_range[0], color='red', linestyle='--', linewidth=0.5)
    ax0.axhline(2 * np.pi / k_max_range[1], color='red', linestyle='-', linewidth=0.5)
    ax0.axhline(2 * np.pi / k_max_range[2], color='red', linestyle='--', linewidth=0.5)

    if pflag:
        plt.ylabel('Wave length\n(meters)')
        plt.legend()

    ax0 = F.fig.add_subplot(gs[2, 1])
    plt.title('Photons density ($m^{-1}$)', loc='left')
    for k in all_beams:
        I = Gk.sel(beam=k)['gFT_PSD_data']
        plt.plot(Gplot.x / 1e3, I.N_photons / I.L.data, label=k, linewidth=0.8)
    plt.plot(Gplot.x / 1e3, G_gFT_wmean.N_per_stancil / 3 / I.L.data, c='black', label='ave Photons', linewidth=0.8)
    plt.xlim(xlims)
    plt.xlabel('Distance from the Ice Edge (km)')

    ax0 = F.fig.add_subplot(gs[2, 2])
    ax0.set_yscale('log')
    plt.title('Peak Spectal Power', loc='left')

    # new B02 files carry an absolute FFT x; old files had x relative to the first stancil
    x0 = 0.0 if ('x_absolute' in Gfft.coords or 'x_absolute' in Gfft.attrs) else Gk.x[0].data
    for k in all_beams:
        I = Gk.sel(beam=k)['gFT_PSD_data']
        plt.scatter(I.x.data / 1e3, I.sel(k=slice(k_max_range[0], k_max_range[2])).integrate('k').data, s=0.5, marker='.', color='red', alpha=0.3)

        I = Gfft.sel(beam=k)
        plt.scatter((x0 + I.x.data) / 1e3, I.power_spec.sel(k=slice(k_max_range[0], k_max_range[2])).integrate('k'), s=0.5, marker='.', c='blue', alpha=0.3)

    Gplot = G_fft_wmean.squeeze()
    Gplot = Gplot.power_spec[:, Gplot.N_per_stancil >= Gplot.N_per_stancil.max().data * 0.9]
    plt.plot((x0 + Gplot.x) / 1e3, Gplot.sel(k=slice(k_max_range[0], k_max_range[2])).integrate('k'), '.', markersize=1.5, c='blue', label='FFT')

    Gplot = G_gFT_wmean.squeeze()
    plt.plot(Gplot.x / 1e3, Gplot.sel(k=slice(k_max_range[0], k_max_range[2])).integrate('k'), '.', markersize=1.5, c='red', label='gFT')

    plt.ylabel('1e-3 $(m)^2~m$')
    plt.legend()

    save_fig(F, plot_path, 'B03_specs_L' + str(Lmeters))

    # %% per-x reconstruction figures
    if 'y_data' not in Gx.sel(beam='gt3r').keys():
        raise SkipTrack('no y_data in gFT_x')

    font_for_print()
    spectra_path = P.track_plot_dir('B03_spectra')
    k_thresh = prm['k_thresh']
    pad = prm['xlim_pad_dx']

    x_pos_sel = np.arange(Gk.x.size)[~np.isnan(Gk.mean('beam').mean('k').gFT_PSD_data.data)]
    # position of the largest mean PSD; argmax runs over the non-nan subset, so map it back to the
    # full x index (the old script used the subset index directly -> empty figure when leading x are nan)
    x_pos_max = x_pos_sel[Gk.mean('beam').mean('k').gFT_PSD_data[~np.isnan(Gk.mean('beam').mean('k').gFT_PSD_data)].argmax().data]
    xpp = x_pos_sel[[int(i) for i in np.round(np.linspace(0, x_pos_sel.size - 1, prm['n_x_examples']))]]
    xpp = np.insert(xpp, 0, x_pos_max)

    n_reconst_figs = 0
    for i in xpp:
        F = M.figure_axis_xy(6, 8, container=True, view_scale=0.8)
        plt.suptitle('gFT Model and Spectrograms | x=' + str(Gk.x[i].data) + ' \n' + track_name, y=0.95)
        gs = GridSpec(5, 6, wspace=0.2, hspace=0.7)

        # full model reconstruction
        ax0 = F.fig.add_subplot(gs[0:2, :])
        neven = True
        offs = 0
        for k in all_beams:
            Gx_1 = Gx.isel(x=i).sel(beam=k)
            Gk_1 = Gk.isel(x=i).sel(beam=k)

            plot_model_eta(Gx_1, ax0, offset=offs, linestyle='-', color=col_dict[k], linewidth=0.4, alpha=1, zorder=12)

            # original data
            eta_1 = plot_data_eta(Gx_1, offset=offs, linestyle='-', c='k', linewidth=1, alpha=0.5, zorder=11)

            # reconstruct in gaps
            FT = gFT.generalized_Fourier(Gx_1.eta + Gx_1.x, None, Gk_1.k)
            _ = FT.get_H()
            FT.p_hat = np.concatenate([Gk_1.gFT_cos_coeff, Gk_1.gFT_sin_coeff])
            plt.plot(Gx_1.eta, FT.model() + offs, '-', c='orange', linewidth=0.3, alpha=1, zorder=2)

            if neven:
                neven = False
                offs += .3
            else:
                neven = True
                offs += 0.6

        dx = eta_1.diff('eta').mean().data
        eta_ticks = np.linspace(Gx_1.eta.data[0], Gx_1.eta.data[-1], 11)
        ax0.set_xticks(eta_ticks)
        ax0.set_xticklabels(eta_ticks / 1e3)
        plt.xlim(eta_1[0].data - pad * dx, eta_1[-1].data + pad * dx)
        plt.title('Model reconst.', loc='left')
        plt.ylabel('relative slope (m/m)')
        plt.xlabel(r'segment distance $\eta$ (km) @ x=' + fltostr(Gx_1.x.data / 1e3, 2) + 'km')

        # spectra per beam pair
        ax1_list = list()
        dd_max = list()
        for pos, kgroup, lflag in zip([gs[2, 0:2], gs[2, 2:4], gs[2, 4:]], [['gt1l', 'gt1r'], ['gt2l', 'gt2r'], ['gt3l', 'gt3r']], [True, False, False]):
            ax11 = F.fig.add_subplot(pos)
            ax11.tick_params(labelleft=lflag)
            ax1_list.append(ax11)
            for k in kgroup:
                Gx_1 = Gx.isel(x=i).sel(beam=k)
                Gk_1 = Gk.isel(x=i).sel(beam=k)
                klim = Gk_1.k[0], Gk_1.k[-1]

                if 'l' in k:
                    dd = Gk_1.gFT_PSD_data
                    plt.plot(Gk_1.k, dd, color='gray', linewidth=.5, alpha=0.5)

                dd = Gk_1.gFT_PSD_data.rolling(k=prm['rolling_k_spec'], min_periods=1, center=True).mean()
                plt.plot(Gk_1.k, dd, color=col_dict[k], linewidth=.8)
                dd_max.append(np.nanmax(dd.data))
                plt.xlim(klim)

                if lflag:
                    plt.ylabel('$(m/m)^2/k$')
                    plt.title('Energy Spectra', loc='left')
                plt.xlabel(r'wavenumber k (2$\pi$ m$^{-1}$)')

            ax11.axvline(k_thresh, linewidth=1, color='gray', alpha=1)
            ax11.axvspan(k_thresh, klim[-1], color='gray', alpha=0.5, zorder=12)

        if ~np.isnan(np.nanmax(dd_max)):
            for ax in ax1_list:
                ax.set_ylim(0, np.nanmax(dd_max) * 1.1)

        # low-wavenumber reconstruction
        ax0 = F.fig.add_subplot(gs[-2:, :])
        neven = True
        offs = 0
        for k in all_beams:
            Gx_1 = Gx.isel(x=i).sel(beam=k)
            Gk_1 = Gk.isel(x=i).sel(beam=k)

            # original data
            eta_1 = plot_data_eta(Gx_1, offset=offs, linestyle='-', c='k', linewidth=1.5, alpha=0.5, zorder=11)

            # reconstruct in gaps
            FT = gFT.generalized_Fourier(Gx_1.eta + Gx_1.x, None, Gk_1.k)
            _ = FT.get_H()
            FT.p_hat = np.concatenate([Gk_1.gFT_cos_coeff, Gk_1.gFT_sin_coeff])

            p_hat_k = np.concatenate([Gk_1.k, Gk_1.k])
            k_mask = p_hat_k < k_thresh
            FT.p_hat[~k_mask] = 0

            plt.plot(Gx_1.eta, FT.model() + offs, '-', c=col_dict[k], linewidth=0.8, alpha=1, zorder=12)

            if neven:
                neven = False
                offs += .3
            else:
                neven = True
                offs += 0.6

        dx = eta_1.diff('eta').mean().data
        eta_ticks = np.linspace(Gx_1.eta.data[0], Gx_1.eta.data[-1], 11)
        ax0.set_xticks(eta_ticks)
        ax0.set_xticklabels(eta_ticks / 1e3)
        n_edge = min(prm['eta_edge_points'], eta_1.size // 4)      # guard against short eta
        plt.xlim(eta_1[n_edge].data - pad * dx, eta_1[-n_edge].data + pad * dx)
        plt.title('Low-Wavenumber Model reconst.', loc='left')
        plt.ylabel('relative slope (m/m)')
        plt.xlabel(r'segment distance $\eta$ (km) @ x=' + fltostr(Gx_1.x.data / 1e3, 2) + 'km')

        save_fig(F, spectra_path, 'B03_freq_reconst_x' + str(i))
        n_reconst_figs += 1

    run.info(n_x=int(Gk.x.size), Lmeters=float(Lmeters), k_max=float(k_max), n_reconst_figs=n_reconst_figs,
             x_positions=[int(i) for i in xpp])
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
