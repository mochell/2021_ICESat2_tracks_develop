# %%
import os, sys
#execfile(os.environ['PYTHONSTARTUP'])

"""
This file open a ICEsat2 track applied filters and corections and returns smoothed photon heights on a regular grid in an .nc file.
This is python 3
"""

exec(open(os.environ['PYTHONSTARTUP']).read())
exec(open(STARTUP_2021_IceSAT2).read())

#%matplotlib inline

import ICEsat2_SI_tools.convert_GPS_time as cGPS
import h5py
import ICEsat2_SI_tools.io as io
import ICEsat2_SI_tools.spectral_estimates as spec

import time
import imp
import copy
import spicke_remover
import datetime
import generalized_FT as gFT
from scipy.ndimage.measurements import label

#import s3fs
# %%
track_name, batch_key, test_flag = io.init_from_input(sys.argv) # loads standard experiment
#track_name, batch_key, test_flag = '20190605061807_10380310_004_01', 'SH_batch01', False
#track_name, batch_key, test_flag = '20190601094826_09790312_004_01', 'SH_batch01', False
#track_name, batch_key, test_flag = '20190207111114_06260210_004_01', 'SH_batch02', False
#track_name, batch_key, test_flag = '20190208152826_06440210_004_01', 'SH_batch01', False
#track_name, batch_key, test_flag = '20190213133330_07190212_004_01', 'SH_batch02', False
#track_name, batch_key, test_flag = '20190207002436_06190212_004_01', 'SH_batch02', False
#track_name, batch_key, test_flag = '20190206022433_06050212_004_01', 'SH_batch02', False

#track_name, batch_key, test_flag = 'SH_20190101_00570212', 'SH_batch04', True
#track_name, batch_key, test_flag = 'SH_20190219_08070210', 'SH_batchminimal', True 



#track_name, batch_key, test_flag = '20190219073735_08070210_004_01', 'SH_batch02', False
#print(track_name, batch_key, test_flag)
hemis, batch = batch_key.split('_')

load_path   = mconfig['paths']['work'] +batch_key+'/B02_spectra/'
load_file   = load_path + 'B02_' + track_name #+ '.nc'
plot_path   = mconfig['paths']['plot'] + '/'+hemis+'/'+batch_key+'/' + track_name + '/'
MT.mkdirs_r(plot_path)

Gk = xr.open_dataset(load_file+'_gFT_k.nc')
Gx = xr.open_dataset(load_file+'_gFT_x.nc')

Gfft = xr.open_dataset(load_file+'_FFT.nc')
# print(Gk)
# print(Gx)
time.sleep(2)


# %%
# for ibeam in Gk.beam:
#     print(Gk.sel(beam=ibeam).gFT_PSD_data.data)

# %%
all_beams   = mconfig['beams']['all_beams']
high_beams  = mconfig['beams']['high_beams']
low_beams   = mconfig['beams']['low_beams']
#Gfilt   = io.load_pandas_table_dict(track_name + '_B01_regridded', load_path) # rhis is the rar photon data
#Gd      = io.load_pandas_table_dict(track_name + '_B01_binned' , load_path)  #
col.colormaps2(21)

# %% check paths (again)

col_dict= col.rels
F = M.figure_axis_xy(9, 3, view_scale  =0.5)

plt.subplot(1,3, 1)
plt.title(track_name , loc ='left')
for k in all_beams:
    I = Gk.sel(beam=k)
    I2 = Gx.sel(beam=k)
    plt.plot(I['lon'], I['lat'], '.', c= col_dict[k],  markersize = 0.7, linewidth  =0.3)
    plt.plot(I2['lon'], I2['lat'], '|', c= col_dict[k],  markersize = 0.7 )


plt.xlabel('lon')
plt.ylabel('lat')

plt.subplot(1,3, 2)

xscale= 1e3
for k in all_beams:
    I = Gk.sel(beam=k)
    plt.plot( I['x_coord']/xscale  , I['y_coord']/xscale, '.'  , c= col_dict[k] , linewidth = 0.8,  markersize = 0.8 )
    # I2 = G_gFT[k]
    # plt.plot( I2.coords['x_coord']/xscale,  I2.coords['y_coord']/xscale, '*' , markersize = 0.7)

plt.xlabel('x_coord (km)')
plt.ylabel('y_coord (km)')

plt.subplot(1,3, 3)

xscale= 1e3
for k in all_beams:
    I = Gk.sel(beam=k)
    plt.plot( I['x_coord']/xscale  , (I['y_coord']-I['y_coord'][0]), '.'  , c= col_dict[k], linewidth = 0.8,  markersize = 0.8)
    # I2 = G_gFT[k]
    # plt.plot( I2.coords['x_coord']/xscale,  I2.coords['y_coord']/xscale, '*' , markersize = 0.7)

plt.xlabel('x_coord (km)')
plt.ylabel('y_coord deviation (m)')


F.save_light(path=plot_path, name = 'B03_specs_coord_check')


# %%
def dict_weighted_mean(Gdict, weight_key):
    """
    returns the weighted meean of a dict of xarray, data_arrays
    weight_key must be in the xr.DataArrays
    """
    #Gdict = G_rar_fft
    #weight_key='N_per_stancil'

    akey = list( Gdict.keys() )[0]
    GSUM = Gdict[akey].copy()
    GSUM.data     = np.zeros(GSUM.shape)
    N_per_stancil = GSUM.N_per_stancil * 0
    N_photons     = np.zeros(GSUM.N_per_stancil.size)

    counter= 0
    for k,I in Gdict.items():
        #print(k)
        I =I.squeeze()
        print(len(I.x) )
        if len(I.x) !=0:
            GSUM                += I.where( ~np.isnan(I), 0) * I[weight_key] #.sel(x=GSUM.x)
            N_per_stancil       += I[weight_key]
        if 'N_photons' in GSUM.coords:
            N_photons    += I['N_photons']
        counter+=1

    GSUM             = GSUM  / N_per_stancil

    if 'N_photons' in GSUM.coords:
        GSUM.coords['N_photons'] = (('x', 'beam'), np.expand_dims(N_photons, 1) )

    GSUM['beam'] = ['weighted_mean']
    GSUM.name='power_spec'

    return GSUM


G_gFT_wmean = (Gk['gFT_PSD_data'].where( ~np.isnan(Gk['gFT_PSD_data']), 0) * Gk['N_per_stancil']).sum('beam')/ Gk['N_per_stancil'].sum('beam')
G_gFT_wmean['N_per_stancil'] = Gk['N_per_stancil'].sum('beam')

G_fft_wmean = (Gfft.where( ~np.isnan(Gfft), 0) * Gfft['N_per_stancil']).sum('beam')/ Gfft['N_per_stancil'].sum('beam')
G_fft_wmean['N_per_stancil'] = Gfft['N_per_stancil'].sum('beam')


# %% plot
def plot_wavenumber_spectrogram(ax, Gi, clev, title= None, plot_photon_density=True ):

    if Gi.k[0] ==0:
        Gi= Gi.sel(k=Gi.k[1:])
    x_lambda= 2 * np.pi/Gi.k
    plt.pcolormesh(Gi.x/1e3, x_lambda , Gi, cmap=plt.cm.ocean_r , vmin = clev[0], vmax = clev[-1])

    ax.set_yscale('log')
    # plt.colorbar(orientation='vertical', pad=0.06, label='Spectral Power (m^2/m)')

    if plot_photon_density:

        plt.plot(Gi.x/1e3, x_lambda[-1] + (Gi.N_per_stancil/Gi.N_per_stancil.max() ) * 10 , c='black', linewidth= 0.8, label='NAN-density' )
        plt.fill_between(Gi.x/1e3, x_lambda[-1] + (Gi.N_per_stancil/Gi.N_per_stancil.max() ) * 10,  0, color='gray', alpha = 0.3)
        ax.axhline(30, color='black', linewidth=0.3)

    #plt.xlabel('Distance from the Ice Edge (km)')
    plt.ylim(x_lambda[-1], x_lambda[0])
    plt.title(title, loc='left')

#Gplot = G.rolling(x=5, min_periods= 1, center=True).mean()
#Gmean = G_gFT_wmean.rolling(x=2, min_periods= 1, center=True).mean()
Gmean = G_gFT_wmean.rolling(k=5, center=True).mean()
#Gmean = Gmean.where(~np.isnan(Gmean), 0)
try:
    k_max_range = Gmean.k[Gmean.isel(x= slice(0, 5)).mean('x').argmax().data].data* 0.75, Gmean.k[Gmean.isel(x= slice(0, 5)).mean('x').argmax().data].data* 1, Gmean.k[Gmean.isel(x= slice(0, 5)).mean('x').argmax().data].data* 1.25
except:
    k_max_range = Gmean.k[Gmean.isel(x= slice(0, 20)).mean('x').argmax().data].data* 0.75, Gmean.k[Gmean.isel(x= slice(0, 20)).mean('x').argmax().data].data* 1, Gmean.k[Gmean.isel(x= slice(0, 20)).mean('x').argmax().data].data* 1.25


# %%
font_for_print()
F = M.figure_axis_xy(6.5, 5.6, container= True, view_scale =1)
Lmeters = Gk.L.data[0]

plt.suptitle('gFT Slope Spectrograms\n' + track_name, y = 0.98)
gs = GridSpec(3,3,  wspace=0.2,  hspace=.5)#figure=fig,
#clev=np.arange(0, 6, 0.1)*3

#%matplotlib inline

# define mean first for colorbar
Gplot = G_gFT_wmean.squeeze().rolling(k=10, min_periods= 1, center=True).median().rolling(x=3, min_periods= 1, center=True).median()
dd = 10 * np.log10(Gplot)
dd= dd.where(~np.isinf(dd), np.nan )
clev_log = M.clevels( [dd.quantile(0.01).data, dd.quantile(0.98).data * 1.2], 31)* 1

#clev = M.clevels( [Gmean.quantile(0.6).data * 1e4, Gmean.quantile(0.99).data * 1e4], 31)/ 1e4

xlims= Gmean.x[0]/1e3, Gmean.x[-1]/1e3

k =high_beams[0]
for pos, k, pflag in zip([gs[0, 0],gs[0, 1],gs[0, 2] ], high_beams, [True, False, False] ):
    ax0 = F.fig.add_subplot(pos)
    Gplot = Gk.sel(beam = k).gFT_PSD_data.squeeze()#.rolling(k=10, x=2, min_periods= 1, center=True).mean()
    #Gplot.mean('x').plot()
    dd2 = 10 * np.log10(Gplot)
    dd2= dd2.where(~np.isinf(dd2), np.nan )
    plot_wavenumber_spectrogram(ax0, dd2,  clev_log, title =k + ' unsmoothed',  plot_photon_density=True )
    plt.xlim(xlims)
    #
    if pflag:
        plt.ylabel('Wave length\n(meters)')
        plt.legend()

for pos, k, pflag in zip([gs[1, 0],gs[1, 1],gs[1, 2] ], low_beams, [True, False, False] ):
    ax0 = F.fig.add_subplot(pos)
    Gplot = Gk.sel(beam = k).gFT_PSD_data.squeeze()#.rolling(k=10, x=2, min_periods= 1, center=True).mean()
    #Gplot.mean('x').plot()
    dd2 = 10 * np.log10(Gplot)
    dd2= dd2.where(~np.isinf(dd2), np.nan )
    plot_wavenumber_spectrogram(ax0, dd2,  clev_log, title =k+ ' unsmoothed',  plot_photon_density=True )
    plt.xlim(xlims)
    #
    if pflag:
        plt.ylabel('Wave length\n(meters)')
        plt.legend()

ax0 = F.fig.add_subplot(gs[2, 0])

plot_wavenumber_spectrogram(ax0, dd, clev_log  , title ='smoothed weighted mean \n10 $\log_{10}( (m/m)^2 m )$', plot_photon_density= True)
plt.xlim(xlims)

# plt.plot(Gplot.x/1e3, 10* nan_list +20 , c='black', label='NAN-density' )
# ax0.axhline(30, color='black', linewidth=0.5)

ax0.axhline(2* np.pi/k_max_range[0], color='red', linestyle= '--', linewidth= 0.5)
ax0.axhline(2* np.pi/k_max_range[1], color='red', linestyle= '-', linewidth= 0.5)
ax0.axhline(2* np.pi/k_max_range[2], color='red', linestyle= '--', linewidth= 0.5)

if pflag:
    plt.ylabel('Wave length\n(meters)')
    plt.legend()

pos = gs[2, 1]
ax0 = F.fig.add_subplot(pos)
plt.title('Photons density ($m^{-1}$)', loc='left')

for k in all_beams:
    I = Gk.sel(beam = k)['gFT_PSD_data']
    plt.plot(Gplot.x/1e3, I.N_photons/I.L.data, label=k, linewidth=0.8)
plt.plot(Gplot.x/1e3, G_gFT_wmean.N_per_stancil/3/I.L.data , c='black', label='ave Photons' , linewidth=0.8)
plt.xlim(xlims)
plt.xlabel('Distance from the Ice Edge (km)')

pos = gs[2, 2]

ax0 = F.fig.add_subplot(pos)
ax0.set_yscale('log')

plt.title('Peak Spectal Power', loc='left')

x0 = Gk.x[0].data
for k in all_beams:
    I = Gk.sel(beam = k)['gFT_PSD_data']
    plt.scatter(I.x.data/1e3, I.sel(k=slice(k_max_range[0], k_max_range[2])).integrate('k').data ,  s=0.5, marker='.', color='red', alpha= 0.3)

    I= Gfft.sel(beam = k)#.to_array()
    #I= I[:,I.N_per_stancil >=  I.N_per_stancil.max().data*0.9]
    plt.scatter( (x0 +I.x.data)/1e3, I.power_spec.sel(k=slice(k_max_range[0], k_max_range[2])).integrate('k')  ,  s=0.5, marker='.', c='blue', alpha= 0.3)


Gplot= G_fft_wmean.squeeze()
Gplot = Gplot.power_spec[:,Gplot.N_per_stancil >=  Gplot.N_per_stancil.max().data*0.9]
plt.plot( (x0 + Gplot.x)/1e3, Gplot.sel(k=slice(k_max_range[0], k_max_range[2])).integrate('k')  , '.', markersize=1.5 , c='blue', label= 'FFT')

Gplot= G_gFT_wmean.squeeze()
plt.plot(  Gplot.x/1e3, Gplot.sel(k=slice(k_max_range[0], k_max_range[2])).integrate('k')   , '.' , markersize=1.5, c='red', label= 'gFT')

plt.ylabel('1e-3 $(m)^2~m$')
plt.legend()
#plt.ylim(Gplot.min()*1.4, Gplot.max()*1.4 )
#plt.xlim(xlims)

F.save_light(path=plot_path, name = 'B03_specs_L'+str(Lmeters))

# %%
Gk.sel(beam = k).gFT_PSD_data.plot()

# %%  define simple routines
def plot_model_eta(D, ax,  offset = 0, xscale= 1e3 , **kargs ):
    eta  = D.eta + D.x
    y_data = D.y_model+offset
    plt.plot(eta/xscale,y_data , **kargs)

    ax.axvline(eta[0].data/xscale , linewidth=2,  color=kargs['color'], alpha=0.5)
    ax.axvline(eta[-1].data/xscale, linewidth=2,  color=kargs['color'], alpha=0.5)

def add_info(D, Dk, ylims):
    eta  = D.eta + D.x
    N_per_stancil, ksize = Dk.N_per_stancil.data , Dk.k.size
    plt.text(eta[0].data, ylims[-1], '  N='+numtostr(N_per_stancil)  + ' N/2M= '+ fltostr(N_per_stancil/2/ksize, 1) )

def plot_data_eta(D,  offset = 0,xscale= 1e3 ,  **kargs ):
    eta_1  = D.eta + D.x
    y_data = D.y_model +offset
    plt.plot(eta_1/xscale,y_data , **kargs)
    return eta_1


# %% phase examples
### overlapping views
#for i in np.arange(0,29,2):
# i = 4
# c1= 'blue'
# c2= 'red'
#
# Gx_1 = Gx.isel(x= i).sel(beam = k)
# Gx_2 = Gx.isel(x= i+1).sel(beam = k)
#
# Gk_1 = Gk.isel(x= i).sel(beam = k)
# Gk_2 = Gk.isel(x= i+1).sel(beam = k)
#
# fltostr = MT.float_to_str
# numtostr = MT.num_to_str
#
# #if k%2 ==0:
# font_for_print()
# F = M.figure_axis_xy(9, 5, container =True, view_scale= 0.8)
#
# plt.suptitle('gFT Slope Spectrograms\n' + track_name, y = 0.98)
# gs = GridSpec(3,4,  wspace=0.2,  hspace=.5)#figure=fig,
#
# ax0 = F.fig.add_subplot(gs[0, :])
#
#
#
# plot_model_eta(Gx_1, ax0, linestyle='-', color=c1, linewidth=0.4, alpha=1, zorder=12 )
# plot_model_eta(Gx_2, ax0, linestyle='-', color=c2, linewidth=0.4, alpha=1, zorder=12 )
#
# ylims= -np.nanstd(Gx_1.y_data)*3, np.nanstd(Gx_1.y_data)*3
#
# add_info(Gx_1, Gk_1 , ylims )
# add_info(Gx_2, Gk_1 , ylims )
#
# # oringial data
#
# eta_1= plot_data_eta(Gx_1 , offset= 0 , linestyle= '-', c='k',linewidth=1,  alpha =0.5, zorder=11)
# eta_2= plot_data_eta(Gx_2  , offset= 0 , linestyle= '-', c='k',linewidth=1,  alpha =0.5, zorder=11)
#
# dx = eta_1.diff('eta').mean()
# plt.xlim(eta_1[0].data - 40 * dx, eta_2[-1].data +  40 * dx )
# plt.ylim(ylims[0], ylims[-1])
#

# %% Single views

def plot_data_eta(D,  offset = 0  ,  **kargs ):
    eta_1  = D.eta# + D.x
    y_data = D.y_model +offset
    plt.plot(eta_1,y_data , **kargs)
    return eta_1

def plot_model_eta(D, ax,  offset = 0,  **kargs ):
    eta  = D.eta #+ D.x
    y_data = D.y_model+offset
    plt.plot(eta ,y_data , **kargs)

    ax.axvline(eta[0].data, linewidth=0.1,  color=kargs['color'], alpha=0.5)
    ax.axvline(eta[-1].data, linewidth=0.1,  color=kargs['color'], alpha=0.5)

if ('y_data' in Gx.sel(beam = 'gt3r').keys()):
    print('ydata is ', ('y_data' in Gx.sel(beam = 'gt3r').keys()) )
else:
    print('ydata is ', ('y_data' in Gx.sel(beam = 'gt3r').keys()) )
    MT.json_save('B03_fail', plot_path, {'reason':'no y_data'})
    print('failed, exit')
    exit()


# %%
fltostr = MT.float_to_str
numtostr = MT.num_to_str

font_for_print()


#for i in x_pos_sel[::2]:
#i =x_pos_sel[20]
MT.mkdirs_r(plot_path+'B03_spectra/')

x_pos_sel =  np.arange(Gk.x.size)[~np.isnan(Gk.mean('beam').mean('k').gFT_PSD_data.data)]
x_pos_max = Gk.mean('beam').mean('k').gFT_PSD_data[~np.isnan(Gk.mean('beam').mean('k').gFT_PSD_data)].argmax().data
xpp = x_pos_sel[ [int(i) for i in np.round(np.linspace(0, x_pos_sel.size-1, 4))]]
xpp = np.insert(xpp, 0, x_pos_max)

for i in xpp:

    #i = xpp[0]
    F = M.figure_axis_xy(6, 8, container =True, view_scale= 0.8)

    plt.suptitle('gFT Model and Spectrograms | x='+str(Gk.x[i].data)+' \n' + track_name, y = 0.95)
    gs = GridSpec(5,6,  wspace=0.2,  hspace=0.7)#figure=fig,

    ax0 = F.fig.add_subplot(gs[0:2, :])
    col_d = col.__dict__['rels']

    neven = True
    offs = 0
    for k in all_beams:

        Gx_1 = Gx.isel(x= i).sel(beam = k)
        Gk_1 = Gk.isel(x= i).sel(beam = k)

        plot_model_eta(Gx_1, ax0, offset= offs,  linestyle='-', color=col_d[k], linewidth=0.4, alpha=1, zorder=12 )
        ylims= -np.nanstd(Gx_1.y_data)*3, np.nanstd(Gx_1.y_data)*3
        #add_info(Gx_1, Gk_1 , ylims )

        # oringial data
        eta_1= plot_data_eta(Gx_1 , offset= offs , linestyle= '-', c='k',linewidth=1,  alpha =0.5, zorder=11)

        # reconstruct in  gaps
        FT = gFT.generalized_Fourier(Gx_1.eta + Gx_1.x, None,Gk_1.k )
        _ = FT.get_H()
        FT.p_hat=np.concatenate([ Gk_1.gFT_cos_coeff, Gk_1.gFT_sin_coeff ])
        plt.plot(Gx_1.eta, FT.model()+offs ,'-', c='orange', linewidth=0.3, alpha=1,zorder= 2)

        if neven:
            neven = False
            offs += .3
        else:
            neven = True
            offs +=0.6


    dx = eta_1.diff('eta').mean().data

    eta_ticks =  np.linspace(Gx_1.eta.data[0], Gx_1.eta.data[-1], 11)

    ax0.set_xticks(eta_ticks)
    ax0.set_xticklabels(eta_ticks/1e3)
    plt.xlim( eta_1[0].data - 40 * dx, eta_1[-1].data+  40 * dx )
    plt.title('Model reconst.', loc ='left')


    plt.ylabel('relative slope (m/m)')
    plt.xlabel('segment distance $\eta$ (km) @ x='+fltostr(Gx_1.x.data/1e3, 2)+'km')


    # spectra
    # define threshold
    k_thresh = 0.085
    ax1_list = list()
    dd_max=list()
    for pos, kgroup, lflag in zip([ gs[2, 0:2],  gs[2, 2:4], gs[2, 4:]],  [['gt1l', 'gt1r'], ['gt2l', 'gt2r'], ['gt3l', 'gt3r']], [True, False, False] ):

        ax11 = F.fig.add_subplot(pos)
        ax11.tick_params(labelleft=lflag)
        ax1_list.append(ax11)
        for k in kgroup:

            Gx_1 = Gx.isel(x= i).sel(beam = k)
            Gk_1 = Gk.isel(x= i).sel(beam = k)

            klim= Gk_1.k[0], Gk_1.k[-1]

            if 'l' in k:
                dd = Gk_1.gFT_PSD_data#.rolling(k=10, min_periods= 1, center=True).mean()
                plt.plot(Gk_1.k,  dd, color='gray', linewidth=.5 ,alpha= 0.5 )

            dd = Gk_1.gFT_PSD_data.rolling(k=10, min_periods= 1, center=True).mean()
            plt.plot(Gk_1.k,  dd, color=col_d[k], linewidth=.8 )
            dd_max.append(np.nanmax(dd.data))
            plt.xlim(klim)

            if lflag:
                plt.ylabel('$(m/m)^2/k$')
                plt.title('Energy Spectra', loc ='left')

            plt.xlabel('wavenumber k (2$\pi$ m$^{-1}$)')

        #plt.ylim(dd.min(), max(dd_max) * 1.1)

        ax11.axvline(k_thresh, linewidth=1,  color='gray', alpha=1)
        ax11.axvspan(k_thresh , klim[-1],   color='gray', alpha=0.5, zorder=12)

    if ~np.isnan(np.nanmax(dd_max)):
        for ax in ax1_list:
            ax.set_ylim(0, np.nanmax(dd_max) * 1.1)

    ax0 = F.fig.add_subplot(gs[-2:, :])

    neven = True
    offs = 0
    for k in all_beams:

        Gx_1 = Gx.isel(x= i).sel(beam = k)
        Gk_1 = Gk.isel(x= i).sel(beam = k)

        #plot_model_eta(Gx_1, ax0, offset= offs,  linestyle='-', color=col_d[k], linewidth=0.4, alpha=1, zorder=12 )
        ylims= -np.nanstd(Gx_1.y_data)*3, np.nanstd(Gx_1.y_data)*3
        #add_info(Gx_1, Gk_1 , ylims )

        # oringial data
        eta_1= plot_data_eta(Gx_1 , offset= offs , linestyle= '-', c='k',linewidth=1.5,  alpha =0.5, zorder=11)

        # reconstruct in  gaps
        FT = gFT.generalized_Fourier(Gx_1.eta + Gx_1.x, None,Gk_1.k )
        _ = FT.get_H()
        FT.p_hat=np.concatenate([ Gk_1.gFT_cos_coeff, Gk_1.gFT_sin_coeff ])

        p_hat_k = np.concatenate([ Gk_1.k, Gk_1.k ])
        k_mask = p_hat_k < k_thresh
        FT.p_hat[~k_mask] = 0

        plt.plot(Gx_1.eta, FT.model()+offs ,'-', c=col_d[k], linewidth=0.8, alpha=1,zorder= 12)

        if neven:
            neven = False
            offs += .3
        else:
            neven = True
            offs +=0.6

    dx = eta_1.diff('eta').mean().data

    eta_ticks =  np.linspace(Gx_1.eta.data[0], Gx_1.eta.data[-1], 11)

    ax0.set_xticks(eta_ticks)
    ax0.set_xticklabels(eta_ticks/1e3)
    plt.xlim( eta_1[1000].data - 40 * dx, eta_1[-1000].data+  40 * dx )
    plt.title('Low-Wavenumber Model reconst.', loc ='left')


    plt.ylabel('relative slope (m/m)')
    plt.xlabel('segment distance $\eta$ (km) @ x='+fltostr(Gx_1.x.data/1e3, 2)+'km')

    F.save_pup(path=plot_path+'B03_spectra/', name = 'B03_freq_reconst_x'+str(i))

MT.json_save('B03_success', plot_path, {'time':'time.asctime( time.localtime(time.time()) )'})
