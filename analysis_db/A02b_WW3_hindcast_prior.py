
import os, sys
#execfile(os.environ['PYTHONSTARTUP'])

"""

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
#import s3fs

col.colormaps2(21)

# %%
track_name, batch_key, test_flag = io.init_from_input(sys.argv) # loads standard experiment
#track_name, batch_key, test_flag = '20190605061807_10380310_004_01', 'SH_batch01', False
#track_name, batch_key, test_flag = '20190601094826_09790312_004_01', 'SH_batch01', False
#track_name, batch_key, test_flag = '20190207111114_06260210_004_01', 'SH_batch02', False
#track_name, batch_key, test_flag = '20190219073735_08070210_004_01', 'SH_batch02', False
track_name_short = track_name[0:-16]



#print(track_name, batch_key, test_flag)
hemis, batch = batch_key.split('_')
#track_name= '20190605061807_10380310_004_01'
ATlevel= 'ATL03'

save_path   = mconfig['paths']['work'] + '/A02_prior_'+hemis+'/'
plot_path   = mconfig['paths']['plot'] + '/'+hemis+'/'+batch_key+'/' + track_name + '/'
save_name   = 'A02b_'+track_name
plot_name    = 'A02b_'+track_name_short
MT.mkdirs_r(plot_path)
MT.mkdirs_r(save_path)
bad_track_path =mconfig['paths']['work'] +'bad_tracks/'+ batch_key+'/'
# %%

all_beams   = mconfig['beams']['all_beams']
high_beams  = mconfig['beams']['high_beams']
low_beams   = mconfig['beams']['low_beams']
#Gfilt   = io.load_pandas_table_dict(track_name + '_B01_regridded', load_path) # rhis is the rar photon data

load_path_WAVE_GLO  = mconfig['paths']['work'] +'/GLOBMULTI_ERA5_GLOBCUR_01/'
file_name_base      = 'LOPS_WW3-GLOB-30M_'


load_path   = mconfig['paths']['work'] +'/B01_regrid_'+hemis+'/'
Gd          = io.load_pandas_table_dict(track_name + '_B01_binned' , load_path)  #

# %%
G1 = dict()
for b in all_beams:
    # find 1st point
    G1[b] = Gd[b].iloc[abs(Gd[b]['lats']).argmin()]

G1 = pd.DataFrame.from_dict(G1).T

dlon_deg = 1 # degree range aroud 1st point
dlat_deg = 30, 5 # degree range aroud 1st point
dlat_deg_prior = 2, 5 # degree range aroud 1st point

dtime = 4 # in hours

lon_range       = G1['lons'].min() - dlon_deg , G1['lons'].max() + dlon_deg
lat_range       = np.sign(G1['lats'].min())*85 , G1['lats'].max() + dlat_deg[1]
lat_range_prior = G1['lats'].min() - dlat_deg_prior[0] , G1['lats'].max() + dlat_deg_prior[1]
timestamp       = pd.to_datetime(G1[['year', 'month', 'day', 'hour', 'minute', 'second']]).mean()
time_range      = np.datetime64(timestamp) - np.timedelta64(dtime, 'h') , np.datetime64(timestamp) + np.timedelta64(dtime, 'h')
#print(time_range)

# create timestamp according to fiels on ftp server:
# time_stamps_search = np.arange(time_range[0].astype('datetime64[3h]') - np.timedelta64(12*24, 'h') , time_range[1].astype('datetime64[3h]') +  np.timedelta64(3, 'h'), np.timedelta64(3, 'h'))
# time_stamps_search_str = [str(t).replace('-', '') for t in time_stamps_search]

# delete to save memory
#del Gd

# %%
import glob

if time_range[0].astype('M8[M]') != time_range[1].astype('M8[M]'): # spanning two years
    MM_str =  str(time_range[0].astype('M8[M]')).replace('-', '_')
    f_list = glob.glob(load_path_WAVE_GLO+'/*'+MM_str+'_'+hemis+'*.nc')

    MM_str =  str(time_range[-1].astype('M8[M]')).replace('-', '_')
    f_list = glob.glob(load_path_WAVE_GLO+'/*'+MM_str+'_'+hemis+'*.nc')
    f_list = f_list + f_list
else:
    MM_str =  str(time_range[0].astype('M8[M]')).replace('-', '_')
    f_list = glob.glob(load_path_WAVE_GLO+'/*'+MM_str+'_'+hemis+'*.nc')

print(f_list)

def sel_data(I, lon_range, lat_range, timestamp = None):
    lon_flag = (lon_range[0] < I.longitude.data) & (I.longitude.data < lon_range[1])
    lat_flag = (lat_range[0] < I.latitude.data) & (I.latitude.data < lat_range[1])
    time_flag = (time_range[0] < I.time.data) & (I.time.data < time_range[1])
    if timestamp is None:
        I = I.isel(latitude = lat_flag, longitude = lon_flag)
    else:
        I = I.interp(time=np.datetime64(timestamp)).isel(latitude = lat_flag, longitude = lon_flag)

    #I = I.isel(latitude = lat_flag, longitude = lon_flag)
    return I

try:

    Gww3 = xr.open_mfdataset(f_list)
    G_beam  = sel_data(Gww3     , lon_range, lat_range      , timestamp).load()
    G_prior = sel_data(G_beam   , lon_range, lat_range_prior)


    # % create Ice mask
    ice_mask = (G_beam.ice > 0) | np.isnan(G_beam.ice)

    # mask all latitudes that are completely full with sea ice.
    lats = list(ice_mask.latitude.data)
    lats.sort(reverse= True)
    #(ice_mask.sum('longitude') == ice_mask.longitude.size).sel(latitude = lats)
    # find 1st latitude that is completely full with sea ice.
    ice_lat_pos = next((i for i, j in enumerate((ice_mask.sum('longitude') == ice_mask.longitude.size).sel(latitude = lats)) if j), None)
    # recreate lat mask based on this criteria
    lat_mask = lats < lats[ice_lat_pos]
    lat_mask = xr.DataArray( lat_mask.repeat(ice_mask.longitude.size ).reshape(ice_mask.shape), dims = ice_mask.dims, coords = ice_mask.coords )
    lat_mask['latitude'] =lats

    ice_mask = ice_mask + lat_mask


    # for vv in list(G_prior.variables.keys()):
    #     print(vv, '  ', G_prior[vv].long_name)
    # %
    #print(FWi.load().time.min().data, FWi.load().time.max().data )
    def draw_range(lon_range, lat_range, *args, **kargs):
        plt.plot( [lon_range[0], lon_range[1], lon_range[1], lon_range[0], lon_range[0]] , [lat_range[0], lat_range[0], lat_range[1], lat_range[1], lat_range[0]] , *args, **kargs)

    # G_beam.dp.plot()
    # G_beam.ice.plot()


    fvar = ['ice',          'dir',                   'dp',       'spr',      'fp',    'hs']
    fcmap =[plt.cm.Blues_r, col.circle_medium_triple, plt.cm.Blues, plt.cm.Blues, plt.cm.Blues, plt.cm.Blues  ]
    fpos = [0, 1, 2, 3, 4, 5]

    font_for_print()
    #plt.rc('pcolor', shading = 'auto')
    F = M.figure_axis_xy(4, 2.5, view_scale= 0.9, container = True)
    plt.suptitle(track_name_short + ' | ' +file_name_base[0:-1].replace('_', ' '), y = 1.4)
    lon, lat= G_beam.longitude, G_beam.latitude

    gs = GridSpec(1,6,  wspace=0.1,  hspace=0.4)#figure=fig,
    #pos0,pos1,pos2 = gs[0:3, 0],gs[3, 0],gs[4, 0]#,gs[3, 0]

    for fv, fp, fc in zip(fvar, fpos, fcmap):

        ax1 = F.fig.add_subplot(gs[0, fp]) #plt.subplot(1, 6, fp)
        if fp ==0:
            ax1.spines['bottom'].set_visible(False)
            ax1.spines['left'].set_visible(False)
            ax1.tick_params(labelbottom=True, bottom=True)

        else:
            ax1.axis('off')

        plt.plot(G1['lons'], G1['lats'], '.r', markersize=5)
        draw_range(lon_range, lat_range_prior, c='red', linewidth = 1, zorder=12)
        draw_range(lon_range, lat_range, c='blue', linewidth = 0.7, zorder=10)
        #G_beam.ice.plot(cmap=plt.cm.Blues_r, )
        if fv != 'ice':
            plt.pcolor(lon, lat,G_beam[fv].where(~ice_mask, np.nan), cmap=fc)
            plt.contour(lon, lat,G_beam.ice, colors= 'black', linewidths = 0.6)
        else:
            plt.pcolor(lon, lat,G_beam[fv], cmap=fc)

        #plt.title(fv, loc='center')
        plt.title(G_beam[fv].long_name.replace(' ', '\n') +'\n'+ fv, loc='left')
        ax1.axis('equal')

    F.save_pup(path= plot_path, name =plot_name+'_hindcast_data')


    # % derive prior:
    #G_beam_masked['dir']
    G_beam_masked = G_beam.where(~ice_mask, np.nan)

    # %
    # make pandas table with obs track end postitions
    key_list = list(G_beam_masked.keys())
    Tend = pd.DataFrame(index = key_list)

    dlist = list()
    for kk in key_list:
        dlist.append( G_beam_masked[kk].mean().data)
    Tend['mean']  = dlist

    dlist = list()
    for kk in key_list:
        dlist.append( G_beam_masked[kk].std().data)
    Tend['std']  = dlist

    dlist = list()
    for kk in key_list:
        dlist.append( G_beam_masked[kk].long_name)
    Tend['name']  = dlist


    Tend = Tend.T
    Tend['lon'] = [ice_mask.longitude.mean().data ,ice_mask.longitude.std().data , 'lontigude']
    Tend['lat'] = [ice_mask.latitude[ice_mask.sum('longitude') ==0].mean().data ,ice_mask.latitude[ice_mask.sum('longitude') ==0].std().data , 'latitude']
    Tend = Tend.T

    # %

    Prior = dict()

    Prior['incident_angle'] = {'value': Tend['mean']['dp'].astype('float') , 'name': Tend['name']['dp']}
    Prior['spread'] = {'value': Tend['mean']['spr'].astype('float') , 'name': Tend['name']['spr']}
    Prior['Hs'] = {'value': Tend['mean']['hs'].astype('float') , 'name': Tend['name']['hs']}
    Prior['peak_period'] = {'value': 1/Tend['mean']['fp'].astype('float') , 'name': '1/' +Tend['name']['fp']}

    Prior['center_lon'] = {'value': Tend['mean']['lon'].astype('float') , 'name': Tend['name']['lon']}
    Prior['center_lat'] = {'value': Tend['mean']['lat'].astype('float') , 'name': Tend['name']['lat']}


    def plot_prior(Prior, axx):
        angle = Prior['incident_angle']['value']
        axx.quiver(Prior['center_lon']['value'], Prior['center_lat']['value'], -np.sin( angle *np.pi/180), - np.cos( angle *np.pi/180) , scale=4.5, zorder =12, width=0.1 ,headlength = 4.5, minshaft=2, alpha = 0.6, color = 'black' )
        axx.plot(Prior['center_lon']['value'], Prior['center_lat']['value'] , '.', markersize= 6, zorder =12, alpha = 1, color = 'black' )
        tstring=  ' ' +str(np.round(  Prior['peak_period']['value'], 1) )+'sec \n   ' +  str( np.round(Prior['Hs']['value'], 1))+'m'
        plt.text(Prior['center_lon']['value'], Prior['center_lat']['value'], tstring)

    font_for_print()
    F = M.figure_axis_xy(2, 4.5, view_scale= 0.9, container = False)

    ax1 = F.ax
    lon, lat= G_beam.longitude, G_beam.latitude
    #gs = GridSpec(1,6,  wspace=0.1,  hspace=0.4)#figure=fig,
    #pos0,pos1,pos2 = gs[0:3, 0],gs[3, 0],gs[4, 0]#,gs[3, 0]
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.tick_params(labelbottom=True, bottom=True)

    plot_prior(Prior, ax1)
    plt.plot(G1['lons'], G1['lats'], '.r', markersize=5)
    draw_range(lon_range, lat_range_prior, c='red', linewidth = 1, zorder=11)
    draw_range(lon_range, lat_range, c='blue', linewidth = 0.7, zorder=10)

    #G_beam.ice.plot(cmap=plt.cm.Blues_r, )
    fv ='ice'
    if fv != 'ice':
        plt.pcolor(lon, lat,G_beam[fv].where(~ice_mask, np.nan), cmap=fc)
        plt.contour(lon, lat,G_beam.ice, colors= 'black', linewidths = 0.6)
    else:
        plt.pcolor(lon, lat,G_beam[fv], cmap=fc)

    #plt.title(fv, loc='center')
    plt.title('Prior\n' + file_name_base[0:-1].replace('_', ' ') +'\n'+ track_name_short, loc='left')
    ax1.axis('equal')
    #F.save_pup(path= plot_path, name =plot_name+'_hindcast_prior')


    MT.save_pandas_table({'priors_hindcast':Tend}, save_name, save_path)

    target_name = 'A02b_'+track_name+'_hindcast_success'

except:
    target_name = 'A02b_'+track_name+'hindcast_fail'

MT.json_save(target_name, save_path, str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")) )
