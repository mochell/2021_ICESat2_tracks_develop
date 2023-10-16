
import os, sys
#execfile(os.environ['PYTHONSTARTUP'])

"""
This file open a ICEsat2 track applied filters and corections and returns smoothed photon heights on a regular grid in an .nc file.
This is python 3
"""

exec(open(os.environ['PYTHONSTARTUP']).read())
exec(open(STARTUP_2021_IceSAT2).read())
xr.set_options(display_style='text')
#%matplotlib inline

import ICEsat2_SI_tools.convert_GPS_time as cGPS
import h5py
import ICEsat2_SI_tools.io as io
import ICEsat2_SI_tools.spectral_estimates as spec

import imp
import copy
import spicke_remover
import datetime
import concurrent.futures as futures

from numba import jit

from ICEsat2_SI_tools import angle_optimizer
import ICEsat2_SI_tools.wave_tools as waves
import concurrent.futures as futures

import time
import ICEsat2_SI_tools.lanczos as lanczos


col.colormaps2(21)

col_dict = col.rels
#import s3fs
# %%
track_name, batch_key, test_flag = io.init_from_input(sys.argv) # loads standard experiment
#track_name, batch_key, test_flag = '20190605061807_10380310_004_01', 'SH_batch01', False
#track_name, batch_key, test_flag = '20190601094826_09790312_004_01', 'SH_batch01', False
#track_name, batch_key, test_flag = '20190207111114_06260210_004_01', 'SH_batch02', False
#track_name, batch_key, test_flag = '20190219073735_08070210_004_01', 'SH_batch02', False
#track_name, batch_key, test_flag = '20190215184558_07530210_004_01', 'SH_batch02', False

# good track
#track_name, batch_key, test_flag = '20190502021224_05160312_004_01', 'SH_batch02', False
#track_name, batch_key, test_flag = '20190502050734_05180310_004_01', 'SH_batch02', False
#track_name, batch_key, test_flag = '20190216200800_07690212_004_01', 'SH_batch02', False

#track_name, batch_key, test_flag = '20190213133330_07190212_004_01', 'SH_batch02', False
#track_name, batch_key, test_flag = '20190219073735_08070210_004_01', 'SH_batch02', False
#track_name, batch_key, test_flag =  'SH_20190224_08800212', 'SH_publish', True

#print(track_name, batch_key, test_flag)
hemis, batch = batch_key.split('_')
#track_name= '20190605061807_10380310_004_01'
ATlevel= 'ATL03'



plot_path   = mconfig['paths']['plot'] + '/'+hemis+'/'+batch_key+'/' + track_name + '/B05_angle/'
MT.mkdirs_r(plot_path)
#bad_track_path =mconfig['paths']['work'] +'bad_tracks/'+ batch_key+'/'
# %%

all_beams   = mconfig['beams']['all_beams']
high_beams  = mconfig['beams']['high_beams']
low_beams   = mconfig['beams']['low_beams']
beam_groups = mconfig['beams']['groups']
group_names = mconfig['beams']['group_names']
#Gfilt   = io.load_pandas_table_dict(track_name + '_B01_regridded', load_path) # rhis is the rar photon data

# load_path   = mconfig['paths']['work'] +'/B01_regrid_'+hemis+'/'
# G_binned    = io.load_pandas_table_dict(track_name + '_B01_binned' , load_path)  #

load_path   = mconfig['paths']['work']+batch_key+'/B02_spectra/'
Gk          = xr.load_dataset(load_path+ '/B02_'+track_name + '_gFT_k.nc' )  #

load_path   = mconfig['paths']['work'] +batch_key+'/B04_angle/'
Marginals   = xr.load_dataset(load_path+ '/B04_'+track_name + '_marginals.nc' )  #

# %% load prior information
load_path   = mconfig['paths']['work']+batch_key+'/A02_prior/'
Prior = MT.load_pandas_table_dict('/A02_'+track_name, load_path)['priors_hindcast']

save_path =  mconfig['paths']['work'] +batch_key+ '/B04_angle/'

# font_for_print()
# F = M.figure_axis_xy(5.5, 3, view_scale= 0.8)
# plt.suptitle(track_name)
# ax1 =  plt.subplot(2, 1, 1)
# plt.title('Data in Beam', loc= 'left')
#
# xi =1

#data = Marginals.isel(x=xi).sel(beam_group= 'group1').marginals
# angle_mask = Marginals.angle[2:-2]
#
#data.T.plot(cmap= plt.cm.OrRd)

# %%


def derive_weights(weights):
    weights = (weights-weights.mean())/weights.std()
    weights = weights - weights.min()
    return weights

def weighted_means(data, weights, x_angle, color='k'):
    """
    weights should have nans when there is no data
    data should have zeros where there is no data
    """
    from scipy.ndimage.measurements import label
    # make wavenumber groups
    groups, Ngroups = label(weights.where(~np.isnan(weights), 0)  )

    for ng in np.arange(1, Ngroups+1):
        wi          = weights[groups == ng]
        weight_norm = weights.sum('k')
        k           = wi.k.data
        data_k      = data.sel(k=k).squeeze()
        data_weight = (data_k * wi)
        plt.stairs(data_weight.sum('k')/ weight_norm , x_angle, linewidth=1 , color ='k')
        if data_k.k.size > 1:
            for k in data_k.k.data:
                plt.stairs(data_weight.sel(k=k) / weight_norm, x_angle, color ='gray', alpha =0.5)

    data_weighted_mean = (data.where( (~np.isnan(data)) & (data != 0), np.nan) * weights ).sum('k')/weight_norm
    return data_weighted_mean




# cut out data at the boundary and redistibute variance
angle_mask = Marginals.angle *0 ==0
angle_mask[0], angle_mask[-1] = False, False
corrected_marginals = Marginals.marginals.isel(angle=angle_mask ) + Marginals.marginals.isel(angle=~angle_mask ).sum('angle')/sum(angle_mask).data

# get groupweights
# ----------------- thius does not work jet.ckeck with data on server how to get number of data points per stancil
#Gx['x'] = Gx.x - Gx.x[0]

# makde dummy variables
M_final      = xr.full_like(corrected_marginals.isel(k=0, beam_group =0).drop('beam_group').drop('k'), np.nan)
M_final_smth = xr.full_like(corrected_marginals.isel(k=0, beam_group =0).drop('beam_group').drop('k'), np.nan)
if M_final.shape[0] > M_final.shape[1]:
    M_final= M_final.T
    M_final_smth= M_final_smth.T
    corrected_marginals=corrected_marginals.T

Gweights = corrected_marginals.N_data
Gweights = Gweights/Gweights.max()

k_mask = corrected_marginals.mean('beam_group').mean('angle')

xticks_2pi = np.arange(-np.pi, np.pi+np.pi/4, np.pi/4)
xtick_labels_2pi = ['-$\pi$', '-$3\pi/4$', '-$\pi/2$','-$\pi/4$','0','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$']

xticks_pi = np.arange(-np.pi/2, np.pi/2+np.pi/4, np.pi/4)
xtick_labels_pi = ['-$\pi/2$','-$\pi/4$','0','$\pi/4$','$\pi/2$',]


font_for_print()
x_list = corrected_marginals.x
for xi in range(x_list.size):


    F = M.figure_axis_xy(7,3.5, view_scale= 0.8, container = True)
    gs = GridSpec(3,2,  wspace=0.1,  hspace=.8)#figure=fig,
    x_str= str(int(x_list[xi]/1e3))

    plt.suptitle('Weighted marginal PDFs\nx='+ x_str +'\n'+io.ID_to_str(track_name), y= 1.05, x = 0.125, horizontalalignment= 'left')
    group_weight = Gweights.isel(x =xi)

    ax_list= dict()
    ax_sum = F.fig.add_subplot(gs[1, 1])
    #ax_sum.tick_params(labelbottom=False)

    ax_list['sum'] = ax_sum

    data_collect = dict()
    for group, gpos in zip(Marginals.beam_group.data, [ gs[0, 0], gs[0, 1], gs[1, 0]] ):
        ax0 = F.fig.add_subplot(gpos)
        ax0.tick_params(labelbottom=False)
        ax_list[group] = ax0

        data    = corrected_marginals.isel(x=xi).sel(beam_group= group)
        weights = derive_weights( Marginals.weight.isel(x=xi).sel(beam_group= group)  )
        weights = weights**2

        # derive angle axis
        x_angle = data.angle.data
        d_angle= np.diff(x_angle)[0]
        x_angle = np.insert(x_angle, x_angle.size , x_angle[-1].data +  d_angle)

        if ( (~np.isnan(data)).sum().data == 0) | (( ~np.isnan(weights)).sum().data == 0):
            data_wmean = data.mean('k')
        else:
            data_wmean = weighted_means(data, weights, x_angle, color= col_dict[group] )
            plt.stairs(data_wmean , x_angle, color =col_dict[group], alpha =1)
        # test if density is correct
        # if np.round(np.trapz(data_wmean) * d_angle, 2) < 0.90:
        #     raise ValueError('weighted mean is not a density anymore')

        plt.title('Marginal PDF '+ group, loc ='left')
        plt.sca(ax_sum)

        # if data_collect is None:
        #     data_collect = data_wmean
        # else:
        data_collect[group] = data_wmean
        #ax0.set_yscale('log')


    data_collect = xr.concat(data_collect.values(), dim='beam_group')
    final_data   = (group_weight * data_collect).sum('beam_group')/group_weight.sum('beam_group').data

    plt.sca(ax_sum)
    plt.stairs( final_data , x_angle, color = 'k', alpha =1, linewidth =0.8)
    ax_sum.set_xlabel('Angle (rad)')
    plt.title('Weighted mean over group & wavenumber', loc='left')

    # get relevant priors
    for axx in ax_list.values():
        axx.set_ylim(0, final_data.max() * 1.5)
        #figureaxx.set_yscale('log')
        axx.set_xticks(xticks_pi)
        axx.set_xticklabels(xtick_labels_pi)

    try:
        ax_list['group3'].set_ylabel('PDF')
        ax_list['group1'].set_ylabel('PDF')
        ax_list['group3'].tick_params(labelbottom=True)
        ax_list['group3'].set_xlabel('Angle (rad)')
    except:
        pass

    ax_final = F.fig.add_subplot(gs[-1, :])
    plt.title('Final angle PDF', loc='left')

    priors_k = Marginals.Prior_direction[ ~np.isnan(k_mask.isel(x= xi))]
    for pk in priors_k:
        ax_final.axvline(pk, color =col.cascade2, linewidth= 1, alpha = 0.7)

    plt.stairs( final_data , x_angle, color = 'k', alpha =0.5, linewidth =0.8)

    final_data_smth = lanczos.lanczos_filter_1d(x_angle,final_data, 0.1)

    plt.plot(x_angle[0:-1], final_data_smth, color = 'black', linewidth= 0.8)

    ax_final.axvline( x_angle[0:-1][final_data_smth.argmax()], color =col.orange, linewidth= 1.5, alpha = 1, zorder= 1)
    ax_final.axvline( x_angle[0:-1][final_data_smth.argmax()], color =col.black, linewidth= 3.2, alpha = 1, zorder= 0)


    plt.xlabel('Angle (rad)')
    plt.xlim(-np.pi*0.8, np.pi*0.8)

    ax_final.set_xticks(xticks_pi)
    ax_final.set_xticklabels(xtick_labels_pi)

    M_final[xi,:] = final_data
    M_final_smth[xi, :] = final_data_smth

    F.save_pup(path = plot_path, name = 'B05_weigthed_margnials_x' + x_str)



M_final.name='weighted_angle_PDF'
M_final_smth.name='weighted_angle_PDF_smth'
Gpdf = xr.merge([M_final,M_final_smth])

if len(Gpdf.x) < 2:
    print('not enough x data, exit')
    MT.json_save('B05_fail', plot_path+'../',  {'time':time.asctime( time.localtime(time.time()) ) , 'reason': 'not enough x segments'})
    print('exit()')
    exit()

# %%
class plot_polarspectra(object):
        def __init__(self,k, thetas, data, data_type='fraction' ,lims=None,  verbose=False):

            """
            data_type       either 'fraction' or 'energy', default (fraction)
            lims            (None) limts of k. if None set by the limits of the vector k
            """
            self.k      =k
            self.data   =data
            self.thetas =thetas

            #self.sample_unit=sample_unit if sample_unit is not None else 'df'
            # decided on freq limit
            self.lims= lims = [self.k.min(),self.k.max()] if lims is None else lims #1.0 /lims[1], 1.0/ lims[0]
            freq_sel_bool=M.cut_nparray(self.k, lims[0], lims[1] )

            self.min=np.round(np.nanmin(data[freq_sel_bool,:]), 2)#*0.5e-17
            self.max=np.round(np.nanmax(data[freq_sel_bool,:]), 2)
            if verbose:
                print(str(self.min), str(self.max) )

            self.klabels=np.linspace(self.min, self.max, 5) #np.arange(10, 100, 20)

            self.data_type=data_type
            if data_type == 'fraction':
                self.clevs=np.linspace(np.nanpercentile(dir_data.data, 1), np.ceil(self.max* 0.9), 21)
            elif data_type == 'energy':
                self.ctrs_min=self.min+self.min*.05
                #self.clevs=np.linspace(self.min, self.max, 21)
                self.clevs=np.linspace(self.min+self.min*.05, self.max*.60, 21)


        def linear(self, radial_axis='period', ax=None, cbar_flag=True):

            """
            """
            if ax is None:
                ax          = plt.subplot(111, polar=True)
                #self.title  = plt.suptitle('  Polar Spectrum', y=0.95, x=0.5 , horizontalalignment='center')
            else:
                ax=ax
            ax.set_theta_direction(-1)  #right turned postive
            ax.set_theta_zero_location("W")

            grid=ax.grid(color='k', alpha=.5, linestyle='-', linewidth=.5)

            if self.data_type == 'fraction':
                cm=plt.cm.RdYlBu_r #brewer2mpl.get_map( 'RdYlBu','Diverging', 4, reverse=True).mpl_colormap
                colorax = ax.contourf(self.thetas,self.k, self.data, self.clevs, cmap=cm, zorder=1)# ,cmap=cm)#, vmin=self.ctrs_min)
            elif self.data_type == 'energy':
                cm=plt.cm.Paired#brewer2mpl.get_map( 'Paired','Qualitative', 8).mpl_colormap
                cm.set_under='w'
                cm.set_bad='w'
                colorax = ax.contourf(self.thetas,self.k, self.data, self.clevs, cmap=cm, zorder=1)#, vmin=self.ctrs_min)
            #divider = make_axes_locatable(ax)
            #cax = divider.append_axes("right", size="5%", pad=0.05)

            if cbar_flag:
                cbar = plt.colorbar(colorax, fraction=0.046, pad=0.1, orientation="horizontal")
                # if self.data_type == 'fraction':
                #     cbar.set_label('Energy Distribution', rotation=0, fontsize=fontsize)
                # elif self.data_type == 'energy':
                #     cbar.set_label('Energy Density ('+self.unit+')', rotation=0, fontsize=fontsize)
                cbar.ax.get_yaxis().labelpad = 30
                cbar.outline.set_visible(False)
                #cbar.ticks.
                clev_tick_names, clev_ticks =MT.tick_formatter(FP.clevs, expt_flag= False, shift= 0, rounder=4, interval=1)
                cbar.set_ticks(clev_ticks[::5])
                cbar.set_ticklabels(clev_tick_names[::5])
                self.cbar  = cbar

            if (self.lims[-1]- self.lims[0]) > 500:
                radial_ticks = np.arange(100, 1600, 300)
            else:
                radial_ticks = np.arange(100, 800, 100)
            xx_tick_names, xx_ticks = MT.tick_formatter( radial_ticks , expt_flag= False, shift= 1, rounder=0, interval=1)
            #xx_tick_names, xx_ticks = MT.tick_formatter( np.arange( np.floor(self.k.min()),self.k.max(), 20) , expt_flag= False, shift= 1, rounder=0, interval=1)
            xx_tick_names = ['  '+str(d)+'m' for d in xx_tick_names]

            ax.set_yticks(xx_ticks[::1])
            ax.set_yticklabels(xx_tick_names[::1])

            degrange    = np.arange(0,360,30)
            degrange    = degrange[(degrange<=80)| (degrange>=280)]
            degrange_label = np.copy(degrange)
            degrange_label[degrange_label > 180] = degrange_label[degrange_label > 180] - 360

            degrange_label = [str(d)+'$^{\circ}$' for d in degrange_label]

            lines, labels = plt.thetagrids(degrange, labels=degrange_label)#, frac = 1.07)

            for line in lines:
                #L=line.get_xgridlines
                line.set_linewidth(5)
                #line.set_linestyle(':')

            #ax.set_yscale('log')
            ax.set_ylim(self.lims)
            ax.spines['polar'].set_color("none")
            ax.set_rlabel_position(87)
            self.ax=ax


font_for_print()
F = M.figure_axis_xy(6, 5.5, view_scale= 0.7, container = True)
gs = GridSpec(8,6,  wspace=0.1,  hspace=3.1)#figure=fig,
col.colormaps2(21)

cmap_spec= plt.cm.ocean_r
clev_spec = np.linspace(-8, -1, 21) *10

cmap_angle= col.cascade_r
clev_angle = np.linspace(0, 4, 21)


ax1 = F.fig.add_subplot(gs[0:3, :])
ax1.tick_params(labelbottom=False)

weighted_spec   = (Gk.gFT_PSD_data * Gk.N_per_stancil).sum('beam') /Gk.N_per_stancil.sum('beam')
x_spec          = weighted_spec.x/1e3
k               = weighted_spec.k

xlims = x_spec[0], x_spec[-1]
#weighted_spec.plot()
#clev_spec = np.linspace(-8, -1, 21) *10
clev_spec = np.linspace(-80, (10* np.log(weighted_spec)).max() * 0.9, 21)

plt.pcolor(x_spec, k, 10* np.log(weighted_spec),vmin= clev_spec[0], vmax= clev_spec[-1],  cmap =cmap_spec )


plt.title(track_name + '\nPower Spectra (m/m)$^2$ k$^{-1}$', loc='left')

cbar = plt.colorbar( fraction=0.018, pad=0.01, orientation="vertical", label ='Power')
cbar.outline.set_visible(False)
clev_ticks = np.round(clev_spec[::3], 0)
#clev_tick_names, clev_ticks =MT.tick_formatter(clev_spec, expt_flag= False, shift= 0, rounder=1, interval=2)
cbar.set_ticks(clev_ticks)
cbar.set_ticklabels(clev_ticks)

plt.ylabel('wavenumber $k$')

#plt.colorbar()
ax2 = F.fig.add_subplot(gs[3:5, :])
ax2.tick_params(labelleft=True)

dir_data = Gpdf.interp(x= weighted_spec.x).weighted_angle_PDF_smth.T

x = Gpdf.x/1e3
angle = Gpdf.angle
plt.pcolor(x_spec, angle, dir_data , vmin= clev_angle[0], vmax= clev_angle[-1], cmap =cmap_angle)
#plt.contourf(x_spec, angle, dir_data ,clev_angle, cmap =cmap_angle)

cbar = plt.colorbar( fraction=0.01, pad=0.01, orientation="vertical", label ='Density')
plt.title('Direction PDF', loc='left')

plt.xlabel('x (km)')
plt.ylabel('angle')

ax2.set_yticks(xticks_pi)
ax2.set_yticklabels(xtick_labels_pi)


x_ticks  = np.arange(0, xlims[-1].data, 50)
x_tick_labels, x_ticks = MT.tick_formatter(x_ticks, expt_flag= False, shift= 0, rounder=1, interval=2)

ax1.set_xticks(x_ticks)
ax2.set_xticks(x_ticks)
ax1.set_xticklabels(x_tick_labels)
ax2.set_xticklabels(x_tick_labels)
ax1.set_xlim(xlims)
ax2.set_xlim(xlims)


xx_list = np.insert(corrected_marginals.x.data, 0, 0)
x_chunks = spec.create_chunk_boundaries( int(xx_list.size/3),  xx_list.size,  iter_flag= False )
x_chunks = x_chunks[:, ::2]
x_chunks[-1, -1] = xx_list.size-1
#x_chunks#.shape

for x_pos, gs in zip( x_chunks.T , [ gs[-3:, 0:2], gs[-3:, 2:4], gs[-3:, 4:]] ):
    #print( x_pos)
    #print( xx_list[x_pos])
    x_range = xx_list[[x_pos[0], x_pos[-1]]]

    ax1.axvline(x_range[0]/1e3, linestyle= ':', color= 'white', alpha = 0.5)
    ax1.axvline(x_range[-1]/1e3, color = 'gray', alpha = 0.5)

    ax2.axvline(x_range[0]/1e3, linestyle= ':', color= 'white', alpha = 0.5)
    ax2.axvline(x_range[-1]/1e3, color = 'gray', alpha = 0.5)


    i_spec  = weighted_spec.sel(x= slice(x_range[0], x_range[-1]) )
    i_dir   = corrected_marginals.sel(x= slice(x_range[0], x_range[-1]) )

    dir_data  = (i_dir * i_dir.N_data).sum([ 'beam_group', 'x'])/ i_dir.N_data.sum([ 'beam_group', 'x'])
    lims = dir_data.k[ (dir_data.sum('angle')!=0) ][0].data, dir_data.k[ (dir_data.sum('angle')!=0)  ][-1].data

    N_angle = i_dir.angle.size
    dir_data2 =  dir_data#.where( dir_data.sum('angle') !=0, 1/N_angle/d_angle )

    plot_data  = dir_data2  * i_spec.mean('x')
    plot_data  = plot_data.rolling(angle =5, k =10).median()#.plot()

    plot_data = plot_data.sel(k=slice(lims[0],lims[-1] ) )
    xx = 2 * np.pi/plot_data.k

    #F = M.figure_axis_xy(5, 4)
    #ax = plt.subplot(1, 1, 1, polar=True)
    #
    if np.nanmax(plot_data.data) != np.nanmin(plot_data.data):

        ax3 = F.fig.add_subplot(gs, polar=True)
        FP= plot_polarspectra(xx, plot_data.angle, plot_data, lims=None , verbose= False, data_type= 'fraction')
        FP.clevs=np.linspace(np.nanpercentile(plot_data.data, 1), np.round(plot_data.max(), 4), 21)
        FP.linear(ax = ax3, cbar_flag=False)
        #FP.cbar.set_label('Energy Density ( (m/m)$^2$ k$^{-1}$ deg$^{-1}$ )', rotation=0, fontsize=10)
        #plt.show()

F.save_pup(path = plot_path + '../', name = 'B05_dir_ov')

# save data
Gpdf.to_netcdf(save_path+ '/B05_'+track_name + '_angle_pdf.nc' )

MT.json_save('B05_success', plot_path + '../', {'time':time.asctime( time.localtime(time.time()) )})
