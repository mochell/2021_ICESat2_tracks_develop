import os, sys
#execfile(os.environ['PYTHONSTARTUP'])

"""
This file open a ICEsat2 track applied filters and corections and returns smoothed photon heights on a regular grid in an .nc file.
This is python 3
"""
# exec(open(os.environ['PYTHONSTARTUP']).read())
# exec(open(STARTUP_2019_DP).read())

exec(open(os.environ['PYTHONSTARTUP']).read())
exec(open(STARTUP_2021_IceSAT2).read())

import datetime
import h5py
from random import sample
import imp
import ICEsat2_SI_tools.convert_GPS_time as cGPS
import ICEsat2_SI_tools.io as io
from spectral_estimates import create_chunk_boundaries_unit_lengths, create_chunk_boundaries
import spectral_estimates as spec
import m_tools_ph3 as MT
import filter_regrid as regrid

from numba import jit

import concurrent.futures as futures

# memory test
#from guppy import hpy


def get_size(x):
    from pympler import asizeof
    ss = asizeof.asizeof(x)/1e6
    return str(ss)


#get_size(hemis)

#import s3fs
#processed_ATL03_20190605061807_10380310_004_01.h5

#imp.reload(io)
track_name, batch_key, test_flag = io.init_from_input(sys.argv) # loads standard experiment
#track_name, batch_key, test_flag = '20190605061807_10380310_004_01', 'SH_batch01', False
#track_name, batch_key, test_flag = '20190207234532_06340210_004_01', 'SH_batch02', False
#track_name, batch_key, test_flag = '20190215184558_07530210_004_01', 'SH_batch02', False
#track_name, batch_key, test_flag = '20190219073735_08070210_004_01', 'SH_batch02', False

track_name, batch_key, test_flag = '20190207235856_06340212_004_01', 'SH_batch02', False
#track_name, batch_key, test_flag = '20190502033317_05170310_004_01', 'SH_batch02', False


# equatorward track
#track_name, batch_key, test_flag = '20190208154150_06440212_004_01', 'SH_batch02', False

# poleward track
#track_name, batch_key, test_flag = '20190209150245_06590210_004_01', 'SH_batch02', False
#conner

#print(track_name, batch_key, test_flag)
hemis, batch = batch_key.split('_')
#track_name= '20190605061807_10380310_004_01'
ATlevel= 'ATL03'

load_path   = mconfig['paths']['scratch'] +'/'+ batch_key +'/'
load_file   = load_path + 'processed_'+ATlevel+'_'+track_name+'.h5'

save_path  = mconfig['paths']['work'] +'/B01_regrid_'+hemis+'/'

plot_path = mconfig['paths']['plot']+ '/'+hemis+'/'+batch_key+'/'+track_name +'/B01/'
bad_track_path =mconfig['paths']['work'] +'bad_tracks/'+ batch_key+'/'
MT.mkdirs_r(save_path)

# set pars

# define parameters:
Lmeter      = 20 # stencil length in units of 'dist'; likely in meters the resulting resolution is L/2
Nphoton_min = 5 # mininum need photons per stancil to return results

Lmeter_large= 100e3 # stancil width for testing photon density. stancils do not overlab.
minium_photon_density = 0.02 # minimum photon density per meter in Lmeter_large chunk to be counted as real signal

plot_flag   = True
Nworkers    = 1        # number of threads for parallel processing # inner loop
Nworkers_process = 6  # number of threads for parallel processing  # outer loop
# %%
# test which beams exist:
all_beams   = mconfig['beams']['all_beams']
high_beams  = mconfig['beams']['high_beams']
# low_beams   = mconfig['beams']['low_beams']


all_beams = high_beams


f         = h5py.File(load_file, 'r')
beams     = [b if b in f.keys() else None for b in all_beams]
imp.reload(regrid)
track_poleward    = regrid.track_pole_ward_file(f)
print('poleward track is ' , track_poleward)
# Load fata and apply height corrections
# This needs version 2 of the ALT 03 dataset

# ATL03       =   h5py.File(load_file, 'r')

#ATL03[k+'/heights/ph_id_channel'][0:100]
#accent = regrid.track_type( B[beams_list[0]] )

# %%

#for k in beams:
def load_data_and_cut(k):
    print(k)

    T, seg = io.getATL03_beam(load_file, beam= k)
    print('loaded')
    T = T[T['mask_seaice']] # only take sea ice points, no ocean points
    T = T.drop(labels=[ 'year', 'month', 'day', 'hour', 'minute', 'second', 'ph_id_count', 'mask_seaice'], axis= 1)
    print( 'T MB '  + get_size(T) )

    ho = k
    ho = MT.add_line_var(ho, 'size', str(T.shape[0]))
    ho = MT.add_line_var(ho, 'by confidence levels:' + str(np.arange(0, 5)), [ (T['signal_confidence'] == i).sum() for i in np.arange(0, 5) ])

    # filter:
    Tsel    = T[ (T['heights']<100) & (T['heights'] > -100) ]# & (T['delta_time']>5) & (T['delta_time']<24) ]

    # if len(Tsel) == 0:
    #     ho  = MT.add_line_var(ho, 'no photons found', '')
    #     #Tsel= T[(T['signal_confidence'] ==-1 ) & (T['heights']<100)  & (T['heights'] > -100) ]# & (T['delta_time']>5) & (T['delta_time']<24) ]

    Tsel_c  = io.getATL03_height_correction(load_file)
    Tsel_c  = Tsel_c[Tsel_c['dem_h'] < 1e5] # cute out weird references
    # needs only dem_h and heihgts
    Tsel = regrid.correct_heights(Tsel, Tsel_c).reset_index(drop=True)# correct height
    print('height corrected')


    ### cut data at the rear that has too much variance
    # cut last segments of data until variance is similar
    rear_mask = np.array(Tsel.index) > -1 # True
    nsize0 = Tsel.shape[0]
    N_seg= 20
    cut_flag = True
    dd_old = -1
    dd0 = np.array(Tsel['heights_c'])
    print('inital length' , nsize0)

    @jit(nopython=True, parallel= False)
    def adjust_length(var_list, rear_mask, cut_flag, track_poleward):

        var_list = var_list if track_poleward else var_list[::-1]

        if var_list[0:3].mean()*10 < var_list[-1]:
            #print('cut last '+ str(100/N_seg) +'% of data')
            rear_mask[int(nsize* (N_seg-1) / N_seg):] = False
        else:
            cut_flag =  False

        rear_mask = rear_mask if track_poleward else rear_mask[::-1]

        return rear_mask, cut_flag

    #@jit(nopython=True, parallel= True)
    def get_var(sti):
        return dd[sti[0]: sti[1]].var()


    while cut_flag:
        dd= dd0[rear_mask]
        nsize = dd.size
        print('new length', nsize)

        stencil_iter = create_chunk_boundaries( int(nsize/N_seg), nsize,ov =0, iter_flag=True )
        var_list = np.array(list(map(get_var, stencil_iter)))
        #print(k, var_list)
        rear_mask, cut_flag = adjust_length(var_list, rear_mask, cut_flag, track_poleward)

        if nsize == dd_old:
            print('--- lengthen segments')
            N_seg -=1
            #cut_flag = False

        dd_old = nsize

    print( 'Tsel MB '  + get_size(Tsel) )

    return k, rear_mask, Tsel, seg, ho

hist    = 'Beam stats'
B       = dict()
B1save  = dict()
SEG     = dict()
#k = beams[0]

# A = list()
# for k in all_beams:
#     A.append(load_data_and_cut(k))

with futures.ProcessPoolExecutor(max_workers=Nworkers_process) as executor:
    A = list( executor.map(load_data_and_cut, all_beams)  )

print( 'A MB '  + get_size(A) )

for I in A: # collect returns from from mapping
    k, rear_mask, Tsel, seg, ho  = I

    Tsel['process_mask']  = rear_mask
    B1save[k]             = Tsel
    B[k]                  = Tsel[rear_mask].drop(columns='process_mask')
    SEG[k]                = seg

    ho      = MT.add_line_var(ho, 'cutted ', str( np.round( 100 *(rear_mask.size -rear_mask.sum() ) / rear_mask.size, 0)) + '% in the back of the track' )
    ho      = MT.add_line_var(ho, 'selected size', str(Tsel.shape[0]))
    #ho      = MT.add_line_var(ho, 'final size ', str(Tsel_c.shape[0]))
    print(ho)
    hist    = MT.write_log(hist, ho)

print('done with 1st loop')

#del T
#del Tsel_c
del Tsel
del A

print( 'B_save MB '  + get_size(B1save) )
print( 'B MB '  +  get_size(B) )
print( 'SEG MB '  +  get_size(SEG) )



# get_size(Tsel)
# Tsel
# Tsel.dtypes
#Tsel.memory_usage()

# %% define x- coodindate

# find earliest segment length that is used.
# this is on the equatorward side for poleward track, or
# on the poleward side for equatorward tracks.

total_segment_dist_x_min = list()
for k,I in SEG.items():
    total_segment_dist_x_min.append( I['segment_dist_x'].min() )
total_segment_dist_x_min = min(total_segment_dist_x_min)

# %%
def make_x_coorindate(k):

    """
    Returns the "true" along track coordindate but finding the correpsonding segment length
    also adds the segment_ID to the main table T
    """
    print(k, ' make coodindate')
    T, seg  = B[k], SEG[k]

    # make sure data is strictly ordered by delta_time
    T   = T.sort_values('delta_time').reset_index(drop=True)

    # find positions where segmetn length is reset
    # shifts segment length postions in time
    delta_onehalf           = seg['delta_time'].diff()/2
    delta_onehalf.iloc[0]   = delta_onehalf.iloc[1]
    seg['delta_half_time']  = seg['delta_time']  - delta_onehalf - 1e-5

    # cur phontos that are not in segmentns
    T           = T[ (T['delta_time'] > seg['delta_half_time'].iloc[0]) & (T['delta_time'] < seg['delta_half_time'].iloc[-1])]
    bin_labels  = np.digitize( T['delta_time'], seg['delta_half_time'], right = True )

    # select relevant data
    SS      = seg['segment_dist_x']
    SS_sid  = seg['segment_id']

    repeats = np.bincount(bin_labels, minlength =SS.shape[0])

    # check if repeats sumup
    if repeats.sum() != T.shape[0]:
        print('repeats do not sum up')

    # repeat  segment dat accoridng to photons
    SS       = SS.repeat(repeats)
    SS.index = T.index
    SS_sid       = SS_sid.repeat(repeats)
    SS_sid.index = T.index

    # define new coordinate
    T['x']          =  SS + T['along_track_distance']
    T['segment_id'] =  SS_sid

    # find bad photons
    def find_anomalie_photons(Ti2, segi):
        x_interp = np.interp(Ti2['delta_time'],  segi['delta_time'], segi['segment_dist_x'] )

        diff_x  = x_interp -  Ti2['x']
        diff_x = abs(diff_x-diff_x.mean())
        return diff_x > 3 *diff_x.std() , x_interp

    fail_mask, x_interp = find_anomalie_photons(T, seg)

    print('weird photon fraction:' ,  sum(fail_mask)/ fail_mask.size)
    # aply fail mask
    #T= T[~fail_mask]
    T['x'][fail_mask] = x_interp[fail_mask]

    return k, T



with futures.ProcessPoolExecutor(max_workers=Nworkers_process) as executor:
    A = list( executor.map(make_x_coorindate, all_beams)  )


# %%
B= dict()
for I in A: # collect returns from from mapping
    k               = I[0]
    B[ k ]          = I[1][::-1]

    if ~track_poleward: # invert x- coordinate if there is an equatorward track
        B[k]            = B[k].reset_index(drop=True)
        B[k]['x_true']  = B[k]['x']
        B[k]['x']       = abs(B[k]['x'] - B[k]['x'].iloc[0])
    else:
        B[k]['x_true']  = B[k]['x']


dist_list   = np.array([np.nan, np.nan])
for k in B.keys():
    dist_list = np.vstack([ dist_list, [  B[k]['x'].iloc[0] , B[k]['x'].iloc[-1] ]  ])

print( 'B MB '  + get_size(B) )

del A
del SEG
#del T


# define latitude limits
# lat_lims, lon_lims, accent = lat_min_max(B, all_beams, accent = None)


# %% test ave photon density abnd quit if necessary
track_dist_bounds   = [ np.nanmin(dist_list[:, 0], 0) , np.nanmax(dist_list[:, 1], 0) ]
length_meter        = abs(track_dist_bounds[1] - track_dist_bounds[0])
#length_meter       = (abs(lat_lims[1])  - abs(lat_lims[0])) * 110e3
p_densities_r       = list()
p_densities_l       = list()

for k, I in B.items():
    if 'r' in k:
        p_densities_r.append( I.shape[0] /length_meter)
    else:
        p_densities_l.append( I.shape[0] /length_meter)

if (np.array(p_densities_l).mean() < 0.5) & (np.array(p_densities_r).mean() < 0.5): # in photons per meter
    print('photon density too low, this track is classified as bad track and pushed to bad list')
    MT.json_save(track_name, bad_track_path, {'densities': [ np.array(p_densities_r).mean(), np.array(p_densities_l).mean()] , 'date': str(datetime.date.today()) })
    print('exit.')
    exit()


# %% save corrected data and delete from cash
io.save_pandas_table(B1save, track_name + '_B01_corrected'  , save_path) # all photos but heights adjusted and with distance coordinate
del B1save

# for testing
# T2 = B['gt1r']
# plt.plot( T2['delta_time'] , T2['along_track_distance'] , '.')
# plt.plot(T2['x'])
# for k,I in B.items():
#     plt.plot( I['x']  , I['across_track_distance'] - I[I['seg_ID_local'] == 0]['across_track_distance'].mean(), '.' , markersize = 0.3)
#     #plt.xlim(3e6, 3.25e6)
#
# for k,I in B.items():
#     plt.plot( I['x']  , I['across_track_distance'], '.' , markersize = 0.3)
#     #plt.xlim(3e6, 3.25e6)


# %%
F = M.figure_axis_xy(4, 3, view_scale = 0.7)

for k,I in B.items():
    plt.plot( I['lats'] ,  I['x']  , '.' , markersize = 0.2)
    #plt.xlim(3e6, 3.25e6)
plt.xlabel('lats')
plt.ylabel('x')
F.save_light(path= plot_path, name='B01_ALT03_'+track_name+'_tracks_check_lat_x')

# %%
F = M.figure_axis_xy(4, 3, view_scale = 0.7)
for k,I in B.items():
    plt.plot( I['delta_time']  , I['lats'], '.' , markersize = 0.3)


plt.xlabel('delta time')
plt.ylabel('lat')
F.save_light(path= plot_path, name='B01_ALT03_'+track_name+'_tracks_check_time_lat')
# %%
F = M.figure_axis_xy(4, 3, view_scale = 0.7)

for k,I in B.items():
    plt.plot( I['delta_time']  , I['x'], '.' , markersize = 0.3)

plt.xlabel('delta time')
plt.ylabel('x')

F.save_light(path= plot_path, name='B01_ALT03_'+track_name+'_tracks_check_time_x')

# %%
##### 1.) derive common axis for beams and filter out low density area at the beginning
print('filter out low density area at the beginning')

#@profile
def derive_axis_and_boundaries(key):
    #key = 'gt3r'
    print(key)
    #T2   #      = regrid.derive_axis(B[key], lat_lims).sort_values('dist')
    T2 = B[key]#['x']

    x0         = get_better_lower_boundary(Lmeter_large, np.array(T2['x']))

    print( 'cut off ' , 100 * (1 - T2[T2['x'] > x0].shape[0]/T2.shape[0])  , '% off all data points at the beginning' )
    T2         = T2[T2['x'] >x0] # cut off low density start

    return key, T2, [T2['x'].min(), T2['x'].max()]

def get_better_lower_boundary(Lmeter_large, dd):

    # T2         = regrid.derive_axis(B[key], lat_lims)
    #dd = np.array(T2['x'])
    stencil_iter = spec.create_chunk_boundaries_unit_lengths( Lmeter_large, [ dd.min(), dd.max()],ov =0, iter_flag= True)

    def get_density(sti):
        return sti[0], sum( (sti[0] <= dd) & (dd < sti[-1]) ) / Lmeter_large

    with futures.ThreadPoolExecutor(max_workers= Nworkers_process) as executor_sub:
        #var_list = np.array(list(map(get_density, stencil_iter)))
        var_list = np.array(list(executor_sub.map(get_density, stencil_iter)))

    var_list = np.array(var_list)
    #var_list[:,0] = np.random.rand(10)
    #sort var_list
    var_list = var_list[var_list[:,0].argsort(), :]
    #print(var_list)
    if sum(var_list[:,1] > minium_photon_density) > 1:
        first_stencil = next((i for i, j in enumerate(var_list[:,1] > minium_photon_density) if j), None) #- 1
        stencil_iter  = spec.create_chunk_boundaries_unit_lengths( Lmeter_large, [ dd.min(), dd.max()],ov =0, iter_flag= False)
        return stencil_iter[0, first_stencil]

    else:
        #first_stencil = next((i for i, j in enumerate(var_list[:,1] > 0) if j), None)# - 1
        print('no sufficient photon density found. return short stencil')
        # first_stencil= len(var_list[:,1]) -1
        # stencil_iter  = spec.create_chunk_boundaries_unit_lengths( Lmeter_large, [ dd.min(), dd.max()],ov =0, iter_flag= False)
        # return stencil_iter[0, first_stencil]

        # #print(first_stencil)
        return var_list[-1,0]#[-1,0]


with futures.ProcessPoolExecutor(max_workers=Nworkers_process) as executor:
    A = list( executor.map(derive_axis_and_boundaries, all_beams)  )


# %%
B2          = dict()
dist_list   = np.array([np.nan, np.nan])
for I in A:
    k         = I[0]
    B2[k]     = I[1]
    #B2[k]['dist'] = B2[k]['x']
    dist_list = np.vstack([dist_list,I[2] ])

del A
del B
track_dist_bounds     = [ np.round( np.nanmin(dist_list[:, 0], 0),1) , np.round(np.nanmax(dist_list[:, 1], 0), 1) ]

print( 'B2 MB '  + get_size(B2) )

# for testing:
#ts_s = np.copy(track_dist_bounds)

# %%
xscale= 1e3
F= M.figure_axis_xy(5, 3, view_scale= 0.6)
for k,I in B2.items():
    plt.plot( I['x']/xscale  , I['across_track_distance']/xscale , '.' , markersize = 0.3)
    #plt.xlim(3e6, 3.25e6)


for k in high_beams:

    Ii = B2[k].iloc[0]
    plt.text(Ii.x/xscale+ 5, Ii.across_track_distance/xscale , k + '\n'+ str(Ii[[ 'lats', 'lons'] ]).split('Name')[0] )

    Ii = B2[k].iloc[-1]
    plt.text(Ii.x/xscale+ 5, Ii.across_track_distance/xscale , str(Ii[[ 'lats', 'lons'] ]).split('Name')[0], ha ='right' )

F.ax.axvline(track_dist_bounds[0]/xscale, color='gray', zorder= 2)
F.ax.axvline(track_dist_bounds[1]/xscale, color='gray', zorder= 2)
F.ax.axhline(0, color='gray', zorder= 2)

plt.title('B01 filter and regrid | ' + track_name +'\npoleward '+str(track_poleward)+' \n \n', loc='left')
plt.xlabel('along track distance (km)')
plt.ylabel('across track distance (km)')

F.save_light(path= plot_path +'../', name='B01_ALT03_'+track_name+'_tracks_all')


# %%
# for testing
#track_dist_bounds[1]  = track_dist_bounds[0] + (track_dist_bounds[1] - track_dist_bounds[0])/20
#track_dist_bounds = ts_s[0] ,ts_s[0] + (ts_s[1] - ts_s[0]) /6




# print(key, Ti.shape,2* Ti.shape[0]/Lmeter)
#
# stencil_iter = create_chunk_boundaries_unit_lengths( Lmeter, track_dist_bounds, iter_flag=False )
# stencil_iter.shape
# print(str(stencil_iter.shape[1]/60/60)+'h')
#
# stencil_iter = create_chunk_boundaries_unit_lengths( Lmeter, track_dist_bounds, iter_flag=True )
#Ti.dtypes

#type = np.float64

#Ti = Ti.astype({ 'heights_c':ftype, 'x':ftype, 'lons':ftype, 'lats':ftype, 'x_true':ftype})
#Ti.dtypes
#Bi = regrid.get_stencil_stats( Ti, stencil_iter, 'heights_c', 'x' , stancil_width= Lmeter/2, Nphoton_min=Nphoton_min, map_func = map)

# with futures.ThreadPoolExecutor(max_workers= 4) as executor_sub:
#     Bi = regrid.get_stencil_stats( Ti, stencil_iter, 'heights_c', 'x' , stancil_width= Lmeter/2, Nphoton_min=Nphoton_min, map_func = executor_sub.map)

# %%
bin_labels= np.searchsorted(stancil_set[0,:][0:50], Ti_sel['x'][0:30])
bin_labels
stencil_center = stancil_set[1,bin_labels-1]


plt.plot(  Ti_sel['x'][0:30], Ti_sel['x'][0:30] *0 , '.')
plt.plot(  stencil_center, stencil_center *0 , 'y.')
plt.plot(stancil_set[0,:][0:20] , stancil_set[0,:][0:20] *0, 'r.')


# %%

T2 = Ti
key_var , key_x_coord ='heights_c', 'x'
def get_stencil_stats_shift( T2, stencil_iter,  key_var , key_x_coord, stancil_width ,  Nphoton_min = 5, plot_flag= False):

    """
    T2              pd.Dataframe with beam data needs at least 'dist' and key
    stencil_iter    np.array that constains the stancil boundaries and center [left boundary, center, right boundary]
    key_var         coloumn index used in T2
    key_x_coord     coloumn index of x coordinate
    stancil_width   width of stencil. is used to normalize photon positions within each stancil.
    Nphoton_min     minimum required photots needed to return meaning full averages

    returns:
    pandas DataFrame with the same as T2 but not taken the median of each column
    the following columns are also added:
    key+ '_weighted_mean'   x-weighted gaussian mean of key for each stencil
    key+ '_mode'            mode of key for each stencil
    'N_photos'              Number of Photons for each stencil
    key+ '_std'             standard deviation for each stencil

    the column 'key' is rename to key+'_median'

    """

    stencil_1       = stencil_iter[:, ::2]
    stencil_1half   = stencil_iter[:, 1::2]


    @jit(nopython=False, parallel= False)
    def weighted_mean(x_rel, y):
        "returns the gaussian weighted mean for stencil"

        #@jit(nopython=True, parallel= False)
        def weight_fnk(x):
            "returns gaussian weight given the distance to the center x"
            return np.exp(- (x/.5)**2 )

        w = weight_fnk(x_rel)
        return np.sum(w*y)/np.sum(w)


    def calc_stencil_stats(group, key,  key_x_coord, stancil_width, stancils):

        "returns stats per stencil"
        #import time
        #tstart = time.time()
        Nphoton     = group.shape[0]
        istancil = group['x_bins'].iloc[int(Nphoton/2)]
        stencil_center = stancils[1, istancil-1]


        if Nphoton > Nphoton_min:

            x_rel   = (group[key_x_coord] - stencil_center)/ stancil_width
            y   = group[key]

            #Tmedian[key+ '_weighted_mean']
            key_weighted_mean = weighted_mean(np.array(x_rel), np.array(y))
            key_std           = y.std()

        else:

            #Nphoton           = 0
            key_weighted_mean = np.nan
            #Tmedian[key+ '_mode']           = np.nan
            key_std            = np.nan

        #Tweight = pd.DataFrame([key_weighted_mean, key_std, Nphoton], index= [key+ '_weighted_mean', key+ '_std', 'N_photos' ])
        Tweight = pd.Series([key_weighted_mean, key_std, Nphoton], index= [key+ '_weighted_mean', key+ '_std', 'N_photos' ])


        #print ( str( istancil) + ' s' + str(time.time() - tstart))
        return Tweight.T

    T_sets = list()
    stancil_set = stencil_1
    for stancil_set in [stencil_1, stencil_1half]:

        # select photons that are in bins
        Ti_sel = T2[  (stancil_set[0,0] < T2['x']) &  (T2['x'] < stancil_set[2,-1]) ]

        # put each photon in a bin
        bin_labels  = np.searchsorted(stancil_set[0,:], Ti_sel['x'])
        #bin_labels2 = np.digitize( Ti_sel['x'], stancil_set[0,:], right = True )

        Ti_sel['x_bins'] =bin_labels
        # group data by this bin
        Ti_g = Ti_sel.groupby(Ti_sel['x_bins'], dropna= False , as_index = True )#.median()

        # take median of the data
        Ti_median = Ti_g.median()

        # apply weighted mean and count photons
        args = [ key_var, key_x_coord, Lmeter/2, stancil_set]

        #%timeit -r 1 -n 1 Ti_weight  = Ti_g.apply(calc_stencil_stats, *args)
        Ti_weight  = Ti_g.apply(calc_stencil_stats, *args)

        #merge both datasets
        T_merged = pd.concat( [Ti_median, Ti_weight], axis= 1)

        # rename columns
        T_merged             =  T_merged.rename(columns={key_var: key_var+'_median', key_x_coord: key_x_coord+ '_median'})
        T_merged[ key_var+  '_median'][ np.isnan(T_merged[key_var+ '_std']) ] = np.nan # replace median calculation with nans

        # set stancil center an new x-coodinate
        T_merged['x'] =  stancil_set[1, T_merged.index-1]

        T_sets.append(T_merged)

    # mergeboth stancils
    T3 = pd.concat(T_sets ).sort_values(by= 'x').reset_index()

    if plot_flag:
        Ti_1, Ti_1half =  T_sets

        plt.plot( Ti_1half.iloc[0:60].x, Ti_1half.iloc[0:60]['heights_c_median'], '.' )
        plt.plot( Ti_1.iloc[0:60].x, Ti_1.iloc[0:60]['heights_c_median'], '.' )
        plt.plot( T3.iloc[0:120].x, T3.iloc[0:120]['heights_c_median'], '-' )


    return T3


Lmeter= 20
k = 'gt1r'
key, Ti = k, B2[k].copy().sort_values('x')
stencil_iter = create_chunk_boundaries_unit_lengths( Lmeter, track_dist_bounds, iter_flag=False )

Bi = get_stencil_stats_shift( Ti, stencil_iter, 'heights_c', 'x' , stancil_width= Lmeter/2, Nphoton_min=Nphoton_min, plot_flag = False)


# Bii = Bi[~np.isnan(Bi['heights_c_std'])]
# plt.hist(Bii['x'] - Bii['x_median'], bins=120 )

# %%
#imp.reload(regrid)
##### 2.) regridding and averaging
print('regrid')
def regridding_wrapper(I):
    key, Ti = I
    print(key, Ti.shape,2* Ti.shape[0]/Lmeter)
    stencil_iter = create_chunk_boundaries_unit_lengths( Lmeter, track_dist_bounds, iter_flag=False )
    #print(str(stencil_iter)+'sec')
    Bi = get_stencil_stats_shift( Ti, stencil_iter, 'heights_c', 'x' , stancil_width= Lmeter/2, Nphoton_min=Nphoton_min)

    #print( 'Bi MB '  + get_size(Bi) )
    print(key, 'done')
    return key, Bi

with futures.ProcessPoolExecutor(max_workers=Nworkers_process) as executor:
    B3 = dict( executor.map(regridding_wrapper, B2.items() )  )

print( 'B3 MB '  + get_size(B3) )
print( 'I MB '  + get_size(I) )
# %% ---define start and end position and same in Json file

I = B3['gt2r'].copy()
D_info = dict()
for k,I in B3.items():

    # reset x coordinate
    I['median_dist']   = I['x_median'] - track_dist_bounds[0] #- Lmeter/2
    I['dist']          = I['x']        - track_dist_bounds[0] #- Lmeter/2
    #I['index']      = I['x']
    # rename y coordinate
    I = I.rename(columns={'across_track_distance': 'y'})

    # find starting and end position
    Di_s  = dict(I[I['segment_id'] == I['segment_id'].iloc[0] ].mean()[['lons', 'lats', 'segment_id', 'delta_time']])
    Di_s['across_track_distance_0'] =track_dist_bounds[0]

    Di_e  = dict(I[I['segment_id'] == I['segment_id'].iloc[-1] ].mean()[['lons', 'lats', 'segment_id', 'delta_time']])
    Di_e['across_track_distance_0'] =track_dist_bounds[0]

    D_info[k] = {'start':Di_s,  'end':Di_e , 'poleward': str(track_poleward) }

    # reorder indexes
    column_names = ['x', 'y', 'x_median', 'median_dist', 'lons', 'lats' ,'heights_c_weighted_mean', 'heights_c_median', 'heights_c_std',  'N_photos', ]
    vars_ad = set(list(I[I['segment_id'] == I['segment_id'].iloc[0] ].mean().index)) - set(column_names)
    I = I.reindex(columns=column_names  + list(vars_ad))

    B3[k] = I

# save Json
MT.json_save(track_name + '_B01_stats',save_path, D_info, verbose= True )

# %% saving data
io.save_pandas_table(B2, track_name + '_B01_regridded'  , save_path) # all photos but heights adjusted and with distance coordinate
io.save_pandas_table(B3, track_name + '_B01_binned'     , save_path) # regridding heights


# %% plotting just for checking
key         = 'gt1r'
if plot_flag:
    MT.mkdirs_r(plot_path)
    T2  = B2[key]
    Ti2 = B3[key]

    dl = 4000

    x_key= 'x'
    latlims = (Ti2[x_key].iloc[0] , Ti2[x_key].iloc[-1] )
    #chunk_list = np.arange(latlims[0],latlims[1], dl )
    #chunk_list = sample(  list(np.arange(latlims[0],latlims[1],dl )[0:80])  ,10)
    chunk_list = np.arange(latlims[0],latlims[1], dl )[::10]
    chunk_list = np.append( chunk_list, latlims[1]- dl-1)
    for ll in chunk_list:
        F = M.figure_axis_xy(7, 3, view_scale=0.8)

        plt.plot( T2[x_key], T2['heights_c'],   'k.',  markersize= 0.5, alpha =0.8 )
        #plt.plot( ALT07['ref']['latitude'] , ALT07['heights']['height_segment_height'] , 'r.', markersize=0.8, alpha = 1, label ='ALT07 seg. heights')

        plt.plot(Ti2[x_key], Ti2['heights_c_weighted_mean'] -0.5, '.-', color='blue', linewidth=0.5, markersize=2,alpha=0.9, label='x-gauss weighted mean +1')
        plt.plot(Ti2[x_key], Ti2['heights_c_median'] +0.5, 'r.-',  linewidth=0.5, markersize=2,alpha=0.9, label='median')

        #plt.plot(Ti2['x'], Ti2['heights_c_mode']-1, 'g.-',  linewidth=0.5, markersize=2,alpha=0.9, label='mode - 1')

        #plt.plot(Ti2['x'], Ti2['heights_c_std'] - 1.8, 'k-', linewidth=0.5,alpha=1)
        #plt.fill_between(  Ti2['x'], Ti2['heights_c_std'] -1.8 , y2=-1.8, color='gray',alpha=1)

        # plt.plot( ALT03['delta_time'], ALT03['heights_c'],   'k.',  markersize= 0.3, alpha =0.2 )
        # plt.plot(ALT07['time']['delta_time'] , ALT07['heights']['height_segment_height'] , 'r.', markersize=0.8, alpha = 1, label ='ALT07 seg. heights')
        plt.legend(loc=1)
        plt.xlim(ll, ll+dl)
        plt.ylim(-4, 4)

        plt.xlabel('Meters from the Sea Ice Edge')
        plt.ylabel('Height Anomalie (meters)')
        F.ax.axhline(y =-1.8, color='black', linewidth=0.5)
        F.save_light(path= plot_path, name='ALT03_filt_compare'+ str(ll))
