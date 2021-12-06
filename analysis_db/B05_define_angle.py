
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

from contextlib import contextmanager
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
track_name, batch_key, test_flag = '20190502021224_05160312_004_01', 'SH_batch02', False
#track_name, batch_key, test_flag = '20190502050734_05180310_004_01', 'SH_batch02', False


#print(track_name, batch_key, test_flag)
hemis, batch = batch_key.split('_')
#track_name= '20190605061807_10380310_004_01'
ATlevel= 'ATL03'



plot_path   = mconfig['paths']['plot'] + '/'+hemis+'/'+batch_key+'/' + track_name + '/'
MT.mkdirs_r(plot_path)
bad_track_path =mconfig['paths']['work'] +'bad_tracks/'+ batch_key+'/'
# %%

all_beams   = mconfig['beams']['all_beams']
high_beams  = mconfig['beams']['high_beams']
low_beams   = mconfig['beams']['low_beams']
beam_groups = mconfig['beams']['groups']
group_names = mconfig['beams']['group_names']
#Gfilt   = io.load_pandas_table_dict(track_name + '_B01_regridded', load_path) # rhis is the rar photon data

# load_path   = mconfig['paths']['work'] +'/B01_regrid_'+hemis+'/'
# G_binned    = io.load_pandas_table_dict(track_name + '_B01_binned' , load_path)  #

load_path   = mconfig['paths']['work'] +'/B02_spectra_'+hemis+'/'
Gk          = xr.load_dataset(load_path+ '/B02_'+track_name + '_gFT_k.nc' )  #

load_path   = mconfig['paths']['work'] + '/B04_angle_'+hemis+'/'
Marginals   = xr.load_dataset(load_path+ '/B04_'+track_name + '_marginals.nc' )  #

# %% load prior information
load_path   = mconfig['paths']['work'] +'/A02_prior_'+hemis+'/'
Prior = MT.load_pandas_table_dict('/A02b_'+track_name, load_path)['priors_hindcast']


#### Define Prior
# Use partitions
Prior2              = Prior.loc[['ptp0','ptp1','ptp2','ptp3','ptp4','ptp5']]['mean']
dominat_period      = Prior2[Prior2.max() ==Prior2]
aa = Prior.loc[['pdp0','pdp1','pdp2','pdp3','pdp4','pdp5']]['mean'].astype('float')
dominant_dir        = waves.get_ave_amp_angle(aa *0+1,aa  )[1]
dominant_dir_spread = Prior.loc[['pspr0','pspr1','pspr2','pspr3','pspr4','pspr5']]['mean'].median()

prior_sel= {'alpha': ( dominant_dir *np.pi/180 , dominant_dir_spread *np.pi/180) } # to radiens
#prior_sel= {'alpha': ( -60 *np.pi/180 , dominant_dir_spread *np.pi/180) } # to radiens

# prior_angle =prior_sel['alpha'][0] * 180/np.pi
# if abs(prior_angle) > 80:
#     print('Prior angle is ', prior_angle, '. quit.')
#     prior_sel['time'] = time.asctime( time.localtime(time.time()) )
#     prior_sel['angle'] =  prior_angle
#     MT.json_save('B04_fail', plot_path, prior_sel)
#     print('exit()')
#     exit()
# %%


# font_for_print()
# F = M.figure_axis_xy(5.5, 3, view_scale= 0.8)
# plt.suptitle(track_name)
# ax1 =  plt.subplot(2, 1, 1)
# plt.title('Data in Beam', loc= 'left')
#
# xi =1

# data = Marginals.isel(x=xi).sel(beam_group= 'group1').sel(angle=angle_mask).marginals
# angle_mask = Marginals.angle[2:-2]
#
# data.T.plot(cmap= plt.cm.OrRd)

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
Gx['x'] = Gx.x - Gx.x[0]
Gweights = Gx.sel(beam =all_beams).sel(x = corrected_marginals.x.data, method='nearest')

Gweights

Gweights = Gweights/Gweights.max()

Gweights_all= dict()
for g, name in zip(beam_groups, ['group1', 'group2', 'group3']):
    Gweights2 =Gweights.sel(beam =g).sum('beam')
    Gweights2.expand_dims('beam_group')
    Gweights2['beam_group'] = name
    Gweights2.name = 'beam_weights'
    Gweights_all[name] = Gweights2
    #Gweights.expand_dims(beam_group=name)

Gweights_all = xr.concat(Gweights_all.values(), dim ='beam_group')


xi =2
group_weight = Gweights_all.sel(x =Marginals.x.isel(x= xi).data )
group_weight = group_weight/ group_weight.max()

font_for_print()
F = M.figure_axis_xy(3.5, 6, view_scale= 0.8, container = True)

gs = GridSpec(4,1,  wspace=0.2,  hspace=.2)#figure=fig,

ax_sum = F.fig.add_subplot(gs[-1, 0])

data_collect = dict()
for group, gpos in zip(Marginals.beam_group.data, [ gs[0, 0], gs[1, 0], gs[2, 0]] ):
    ax0 = F.fig.add_subplot(gpos)
    ax0.tick_params(labelbottom=False)

    data    = corrected_marginals.isel(x=xi).sel(beam_group= group)
    weights = derive_weights( Marginals.weight.isel(x=xi).sel(beam_group= group)  )
    weights = weights**2

    # derive angle axis
    x_angle = data.angle.data
    d_angle= np.diff(x_angle)[0]
    x_angle = np.insert(x_angle, x_angle.size , x_angle[-1].data +  d_angle)

    data_wmean = weighted_means(data, weights, x_angle, color= col_dict[group] )
    plt.stairs(data_wmean , x_angle, color =col_dict[group], alpha =1)
    # test if density is correct
    if np.round(np.trapz(data_wmean) * d_angle, 2) < 0.90:
        raise ValueError('weighted mean is not a density anymore')

    plt.sca(ax_sum)
    plt.stairs(data_wmean , x_angle, color =col_dict[group], alpha =0.6)
    # if data_collect is None:
    #     data_collect = data_wmean
    # else:
    data_collect[group] = data_wmean


data_collect = xr.concat(data_collect.values(), dim='beam_group')

plt.sca(ax_sum)
final_data = (group_weight * data_collect).sum('beam_group')/group_weight.sum('beam_group')
plt.stairs( final_data , x_angle, color = 'k', alpha =1)



# %%
#k_list = data.k[~np.isnan(data.mean('angle'))]
#for k in k_list:


k_list = data.k[~np.isnan(data.mean('angle'))]
for k in k_list:
    a = data.angle
    a = np.insert(a.data, a.data.size , a[-1 ].data + np.diff(a)[0] )
    plt.stairs(data.sel(k=k) * weights.sel(k=k)**2 , a )


#plt.pcolormesh(data_mask.x/1e3, data_mask.beam, data_mask, cmap= plt.cm.OrRd)
for i in np.arange(1.5, 6, 2):
    ax1.axhline(i, color= 'black', linewidth =0.5)
plt.xlabel('Distance from Ice Edge')

# ax2 = plt.subplot(2, 1, 2)
# plt.title('Data in Group', loc= 'left')
# plt.pcolormesh(data_mask.x/1e3, data_mask_group.beam_group, data_mask_group, cmap= plt.cm.OrRd)
#
# for i in np.arange(0.5, 3, 1):
#     ax2.axhline(i, color= 'black', linewidth =0.5)
#
# plt.plot( x_list/1e3, x_list*0 +0, '.', markersize= 2, color= col.cascade1 )
# plt.plot( x_list/1e3, x_list*0 +1, '.', markersize= 2, color= col.cascade1 )
# plt.plot( x_list/1e3, x_list*0 +2, '.', markersize= 2, color= col.cascade1 )
#
# plt.xlabel('Distance from Ice Edge')
#
# F.save_pup(path= plot_path, name = 'B04_data_avail')










# %% plot
font_for_print()
F = M.figure_axis_xy(6, 5.5, view_scale= 0.7, container = True)

gs = GridSpec(4,6,  wspace=0.2,  hspace=.8)#figure=fig,

ax0 = F.fig.add_subplot(gs[0:2, -1])
ax0.tick_params(labelleft=False)

#klims = k_list.min()*0.2 , k_list.max()*1.2
klims = 0, k_list.max()*1.2

for g in MM.beam_group:
    MMi = MM.sel(beam_group= g)
    plt.plot(  MMi.weight.T,MMi.k,   '.', color= col_dict[str(g.data)], markersize= 3, linewidth = 0.8)

plt.xlabel('Power')
plt.ylim(klims)

ax1 = F.fig.add_subplot(gs[0:2 , 0:-1])

for g in MM.beam_group:
    Li = LL.loc[str(g.data)]

    angle_list = np.array(Li['alpha']) *  180 /np.pi
    kk_list = np.array(Li['K_prime'])
    weight_list_i = np.array(Li['weight'])

    plt.scatter( angle_list, kk_list, s= (weight_list_i*8e1)**2 , c=col_dict[str(g.data)], label ='mode ' + str(g.data) )


lflag= 'paritions ww3'
for i in np.arange(6):
    i_dir, i_period = Prior.loc['pdp'+ str(i)]['mean'], Prior.loc['ptp'+ str(i)]['mean']
    i_k = (2 * np.pi/ i_period)**2 / 9.81
    i_dir = [i_dir -360 if i_dir > 180 else i_dir][0]
    i_dir = [i_dir +360 if i_dir < -180 else i_dir][0]

    plt.plot(i_dir, i_k,  '.', markersize = 6, color= col.black, label= lflag)
    lflag = None


#ax1.axvline(  best_guess * 180/ np.pi , color=col.blue, linewidth = 1.5, label ='best guess fitting')

ax1.axvline(  (prior_sel['alpha'][0])  *  180 /np.pi, color='k', linewidth = 1.5, label ='prior')
ax1.axvline(  (prior_sel['alpha'][0]- prior_sel['alpha'][1])  *  180 /np.pi, color='k', linewidth = 0.7, label ='prior uncertrainty')
ax1.axvline(  (prior_sel['alpha'][0]+ prior_sel['alpha'][1]) *  180 /np.pi , color='k', linewidth = 0.7)

plt.legend()
plt.ylabel('wavenumber (deg)')
#plt.xlim(- 170, 170)
plt.xlim(- 90, 90)
plt.ylim(klims)

prior_angle_str =str(np.round( (prior_sel['alpha'][0])  *  180 /np.pi))
plt.title(track_name + '\nprior=' + prior_angle_str + 'deg', loc= 'left' )

ax3 = F.fig.add_subplot(gs[2 , 0:-1])

for g in MM.beam_group:
    MMi = MM.sel(beam_group= g)
    wegihted_margins = ( (MMi.marginals * MMi.weight).sum(['x','k'] )/MMi.weight.sum(['x', 'k']) )
    plt.plot(  MMi.angle * 180/ np.pi, wegihted_margins , '.', color= col_dict[str(g.data)], markersize= 2, linewidth = 0.8)

plt.ylabel('Density')
plt.xlabel('Angle (deg)')
plt.title('weight margins', loc='left')

#plt.plot(marginal_stack.angle *  180 /np.pi, marginal_stack.T ,  c=col.gray, label ='weighted mean BF')

#plt.plot(cost_wmean.angle *  180 /np.pi, cost_wmean ,  c=col.rascade3, label ='weighted mean BF')
plt.xlim(- 90, 90)
#plt.xlim(- 125, 125)

ax3 = F.fig.add_subplot(gs[-1 , 0:-1])

for g in MM.beam_group:
    MMi = MM.sel(beam_group= g)
    wegihted_margins = MMi.marginals.mean(['x','k'] )# ( (MMi.marginals * MMi.weight).sum(['x','k'] )/MMi.weight.sum(['x', 'k']) )
    plt.plot(  MMi.angle * 180/ np.pi, wegihted_margins , '.', color= col_dict[str(g.data)], markersize= 2, linewidth = 0.8)

plt.ylabel('Density')
plt.xlabel('Angle (deg)')
plt.title('unweighted margins', loc='left')

#plt.plot(marginal_stack.angle *  180 /np.pi, marginal_stack.T ,  c=col.gray, label ='weighted mean BF')

#plt.plot(cost_wmean.angle *  180 /np.pi, cost_wmean ,  c=col.rascade3, label ='weighted mean BF')
plt.xlim(- 90, 90)
#plt.xlim(- 125, 125)

F.save_pup(path= plot_path, name = 'B04_marginal_distributions')

MT.json_save('B04_success', plot_path, {'time':'time.asctime( time.localtime(time.time()) )'})

# %%
# # %% define weights
# #weights = xr.DataArray(Z_max_list_gauss, dims ='k', coords = {'k': k_list_gauss})
# weights_costs = xr.DataArray(weight_list, dims ='k', coords = {'k': k_list})
# weights_costs = weights_costs/weights_costs.max()
# #weights_costs = weights_costs.sel(k= cost_stack_rolled.k)
#
# weight_k = (1/weights_costs.k)
# weight_k = weight_k/weight_k.max()
#
# weights_sum= weights_costs# * weight_k.data**2
# plt.plot(weights_sum)
#
#
#
# cost_wmean  = (weights_sum * marginal_stack /weights_sum.sum() ).sum('k')
#
# cost_wmean.T.plot()
# best_guess = cost_wmean.angle[cost_wmean.argmax()].data# * 180/ np.pi
# best_guess
