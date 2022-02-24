import os, sys
#execfile(os.environ['PYTHONSTARTUP'])

"""
This script opens an ATL07 track and tests if there is sufficient data and maybe waves:
1) disects each beam such that they start in the MIZ and end at the pole.
2) checks if the variance decays from the MIZ poleward if data density is high enough.
3) generates plausible ATL03 track_names
4) Saves dummy files tracknames and reporting A01b_success_'track_name'.json file

"""

# exec(open(os.environ['PYTHONSTARTUP']).read())
# exec(open(STARTUP_2019_DP).read())
sys.path
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

import concurrent.futures as futures

import piecewise_regression

#import s3fs
#processed_ATL03_20190605061807_10380310_004_01.h5

#imp.reload(io)
track_name, batch_key, test_flag = io.init_from_input(sys.argv) # loads standard experiment

# track NH
#track_name, batch_key, test_flag = '20190301004639_09560201_005_01', 'NH_batch05', False
track_name, batch_key, test_flag = '20190101180843_00660201_005_01', 'NH_batch05', False

# track SH
#track_name, batch_key, test_flag = '20190101084259_00600201_005_01', 'SH_batch04', False
#track_name, batch_key, test_flag = '20190102130012_00780201_005_01', 'SH_batch04', False

#print(track_name, batch_key, test_flag)
hemis, batch = batch_key.split('_')
#track_name= '20190605061807_10380310_004_01'

ATlevel= 'ATL07-02' if hemis == 'SH' else 'ATL07-01'

load_path   = mconfig['paths']['scratch'] +'/'+ batch_key +'/'
load_file   = load_path + ATlevel+'_'+track_name+'.h5'

save_path  = mconfig['paths']['work'] +'/'+ batch_key +'/'+'/A01b_regrid_'+hemis+'/'
plot_path = mconfig['paths']['plot']+ '/'+hemis+'/'+batch_key+'/'+track_name +'/A01b/'
#bad_track_path =mconfig['paths']['work'] +'bad_tracks/'+ batch_key+'/'
MT.mkdirs_r(save_path)

plot_flag   = False
# %%
# test which beams exist:
all_beams = mconfig['beams']['all_beams']


try:
    f     = h5py.File(load_file, 'r')
except:
    print('file not found, exit')
    MT.json_save(name='A01b_success_'+track_name, path=save_path, data= {'reason':'ATL07 file not found, exit'})
    exit()

beams     = [b if b in f.keys() else None for b in all_beams]
imp.reload(regrid)
#track_poleward    = regrid.track_pole_ward_file(f, product='ATL10')
# print('poleward track is' , track_poleward)
# ATL03       =   h5py.File(load_file, 'r')
# ATL03['orbit_info'].keys()
# ATL03['orbit_info/lan'][:]
#
# ATL03['gt1l/freeboard_beam_segment/height_segments'].keys()

# %%
def cut_rear_data(xx0, dd0, N_seg= 20):
    """
    returns masks that cuts large variance in the back of the data
    """
    rear_mask = xx0*0 > -1 # True
    nsize0 = rear_mask.size

    cut_flag = True
    dd_old = -1

    print('inital length' , nsize0)

    #@jit(nopython=True, parallel= False)
    def adjust_length(var_list, rear_mask, cut_flag):

        #var_list = var_list if track_poleward else var_list[::-1]

        if var_list[0:3].mean()*2 < var_list[-1]:
            #print('cut last '+ str(100/N_seg) +'% of data')
            rear_mask[int(nsize* (N_seg-1) / N_seg):] = False
        else:
            cut_flag =  False

        #rear_mask = rear_mask if track_poleward else rear_mask[::-1]

        return rear_mask, cut_flag

    def get_var(sti):
        return dd[sti[0]: sti[1]].var()

    while cut_flag:
        dd= dd0[rear_mask]
        nsize = dd.size
        print('new length', nsize)
        if (nsize/N_seg) < 1:
            break
        stencil_iter = create_chunk_boundaries( int(nsize/N_seg), nsize,ov =0, iter_flag=True )
        var_list = np.array(list(map(get_var, stencil_iter)))
        #print(k, var_list)
        rear_mask, cut_flag = adjust_length(var_list, rear_mask, cut_flag)

        if nsize == dd_old:
            print('--- lengthen segments')
            N_seg -=1
            #cut_flag = False

        dd_old = nsize

    return rear_mask

def get_breakingpoints(xx, dd ,Lmeter= 3000):

    nsize = dd.size
    stencil_iter = spec.create_chunk_boundaries_unit_lengths( Lmeter, [ xx.min(), xx.max()],ov =Lmeter*3/4, iter_flag= True)
    iter_x = spec.create_chunk_boundaries_unit_lengths( Lmeter, [ xx.min(), xx.max()],ov =Lmeter*3/4, iter_flag= False)[1,:]

    def get_var(sti):
        mask = (sti[0] < xx) & (xx <= sti[1])
        return np.nanvar(dd[mask])

    var_list = np.array(list(map(get_var, stencil_iter)))

    x2, y2 =  iter_x/1e3, var_list

    x2= x2[~np.isnan(y2)]
    y2= y2[~np.isnan(y2)]

    convergence_flag =True
    n_breakpoints= 1
    while convergence_flag:
        pw_fit = piecewise_regression.Fit(x2, y2, n_breakpoints=1)
        print('n_breakpoints', n_breakpoints, pw_fit.get_results()['converged'])
        convergence_flag = not pw_fit.get_results()['converged']
        n_breakpoints += 1
        if n_breakpoints == 4:
            convergence_flag = False

    pw_results = pw_fit.get_results()
    if pw_results['converged']:
        if pw_results['estimates']['alpha1']['estimate'] < 0:
            print('decay at the front')
            print('n_breakpoints',pw_fit.n_breakpoints )

        breakpoint = pw_results['estimates']['breakpoint1']['estimate']
        return pw_results['estimates']['alpha1']['estimate'], pw_fit, breakpoint

    else:
        return np.nan, pw_fit, False

DD_slope  = pd.DataFrame(index =beams, columns= ['TF1', 'TF2'])
DD_data   = pd.DataFrame(index =beams, columns= ['TF1', 'TF2'])
DD_region = pd.DataFrame(index =beams, columns= ['TF1', 'TF2'])
DD_region[:] = (np.nan)
DD_pos_start = pd.DataFrame(index =beams, columns= ['TF1_lon', 'TF1_lat', 'TF2_lon', 'TF2_lat'])
DD_pos_end   = pd.DataFrame(index =beams, columns= ['TF1_lon', 'TF1_lat', 'TF2_lon', 'TF2_lat'])

plot_flag = True
for k in beams:

    #k = beams[0]
    #imp.reload(io)
    print(k)
    try:
        T_freeboard = io.getATL07_beam(load_file, beam= k)
    except:
        print('failed to load beam')
        slope_test = False
        data_density  = False
        #return data_density, slope_test
        print('break -------', k, TF,  data_density, slope_test)
        continue



    # find devide such that each hemisphere is split into two parts, if data is there
    if hemis == 'SH':
        ###### for SH tracks
        lon_list = T_freeboard['ref']['longitude']
        mask1 = (lon_list[0]-5 < lon_list) & (lon_list < lon_list[0]+5)
        mask2 = ~mask1
        tot_size =T_freeboard['ref']['latitude'].shape[0]

    else:
        ###### for NH tracks
        from scipy.ndimage.measurements import label
        mask1 = label(T_freeboard['ref']['latitude'] < 88)[0] ==1
        mask2 = ~mask1
        tot_size =T_freeboard['ref']['latitude'].shape[0]

    # cut data accordingly
    if (sum(mask1)/tot_size < 0.05) or (sum(mask1) < 1000):
        TF1 = None
        TF2 = T_freeboard[mask2]
    elif (sum(mask2)/tot_size < 0.05) or (sum(mask2) < 1000):
        TF1 = T_freeboard[mask1]
        TF2 = None
    else:
        TF1 = T_freeboard[mask1]
        TF2 = T_freeboard[mask2]

    # plot splits
    plt.figure()
    plt.plot(TF1['ref']['latitude'], TF1['ref']['longitude'], 'r.', label ='TF1')
    if TF2 is not None:
        plt.plot(TF2['ref']['latitude'], TF2['ref']['longitude'], 'b.', label ='TF2')
    plt.title(track_name + ' '+hemis, loc= 'left')
    plt.legend()
    plt.xlabel('Latitude')
    plt.ylabel('longitude')
    M.save_anyfig(plt.gcf(), name='A01b_'+track_name+'_'+ hemis+'_'+k  , path=plot_path)
    plt.close()
    #plt.show()

    # check if sub-taable goes equatorward or not, then sort accordingly and define along-track axis
    def pole_ward_table(T):
        """
        Returns true if table goes poleward
        hdf5_file is a an HFD5 object in read mode
        """

        time = T['time']['delta_time']
        lat = T['ref']['latitude']
        print('1st lat =' + str(abs(lat.iloc[time.argmin()])) , ';last lat =' + str(abs(lat.iloc[time.argmax()])) )

        return abs(lat.iloc[time.argmax()]) > abs(lat.iloc[time.argmin()])


    TF1_poleward = pole_ward_table(TF1)
    if TF2 is not None:
        TF2_poleward =  pole_ward_table(TF2)

        if TF1_poleward & TF2_poleward:
            raise ValueError('both parts are acending or decending')
    else:
        TF2_poleward = not TF1_poleward

    # flip the beam section that is not poleward
    if TF1_poleward & (TF2 is not None):
        print('TF2 poleward is ', TF2_poleward)
        TF2 = TF2.sort_values(('ref','seg_dist_x'), ascending=False).reset_index()
    else:
        print('TF1 polewards is ', TF2_poleward)
        TF1 = TF1.sort_values(('ref','seg_dist_x'), ascending=False).reset_index()

    # create local x axis
    TF1['x'] = abs(TF1['ref']['seg_dist_x'] -TF1['ref']['seg_dist_x'].iloc[0])
    if TF2 is not None:
        TF2['x'] = abs(TF2['ref']['seg_dist_x'] -TF2['ref']['seg_dist_x'].iloc[0])

    # assign Region to each subset, hemisphere dependent
    for TF,Tsel,TF_polward in zip(['TF1', 'TF2'], [TF1, TF2], [TF1_poleward, TF2_poleward]):

        print(TF,TF_polward)
        if (hemis == 'SH') & TF_polward:
            region = ('10') # SO region
        elif (hemis == 'SH') & (not TF_polward):
            region = ('12') # SO region
        elif (hemis == 'NH') & (TF_polward):
            region = ('03','04') # assign subarctic and high-arctic region
        elif (hemis == 'NH') & (not TF_polward):
            region = ('05','04') # assign subarctic and high-arctic region
        else:
            region =False

        # cut high sigma values
        if (Tsel is None):
            slope_test = False
            data_density  = False
            #return data_density, slope_test
            print('break -------', k, TF,  data_density, slope_test)
            continue

        # ignore bad segments
        Tsel = Tsel[Tsel['heights']['height_segment_surface_error_est'] < 1e2]
        if (Tsel.size <= 50):
            #print('too small table, skip ')
            Tsel = None

        # if Tsel is None skip this itteration
        if Tsel is None:
            slope_test = False
            data_density  = False
            #return data_density, slope_test
            print('break -------', k, TF,  data_density, slope_test)
            continue
        else:
            data_density =Tsel.shape[0]/abs(Tsel['x'].max() - Tsel['x'].min()) # datapoints per meters

        # plt.plot(Tsel['x']/1e3, Tsel['ref']['beam_fb_sigma'], '.k')
        # plt.plot(Tsel['x']/1e3, Tsel['ref']['beam_fb_height'], '.r')
        # plt.plot(Tsel['x']/1e3, Tsel['freeboard']['height_segment_height'], '.')
        # #plt.xlim(60,70)
        # plt.ylim(-1, 5)

        # % cut data in the back: only usefull for SH:
        xx0, dd0 = np.array(Tsel['x']), np.array(Tsel['heights']['height_segment_height'])
        if hemis is 'SH':
            # cut data around the contiental margin
            rear_mask = cut_rear_data(xx0, dd0)
        else:
            # assume all data points are valid for NH ...
            rear_mask = np.array(Tsel['x'] > -1)

        #print('density post cutting', len(xx0[rear_mask])/abs(xx0[rear_mask].max() - xx0[rear_mask].min()) )

        # if cutted data is too short, skip loop
        if len(xx0[rear_mask]) < 500:
            slope_test = False
            data_density  = False
            #return data_density, slope_test
            print('break -------', k, TF,  data_density, slope_test)
            continue

        # estmiate slope at the beginning
        slope_test, pw_fit, breakpoint = get_breakingpoints(xx0[rear_mask], dd0[rear_mask], Lmeter= 3000)

        if plot_flag:
            plt.figure()
            plt.plot(xx0[rear_mask]/1e3, dd0[rear_mask], '.k', markersize= 0.4)
            pw_fit.plot()
            plt.title(k +' '+ TF + ',  data=' +str(data_density) +  ', slope='+str(slope_test)  +  '\n' + track_name , loc= 'left')
            M.save_anyfig(plt.gcf(), name='A01b_'+track_name+'_'+k +'_'+ TF , path=plot_path)
            plt.close()
            #plt.show()


        # assign to tables
        DD_slope.loc[ k, TF] = slope_test
        DD_data.loc[  k, TF] = data_density
        DD_region.loc[k, TF] = region
        DD_pos_start.loc[k, [TF+'_lon', TF+'_lat']]  =  Tsel.iloc[0]['ref']['longitude']  , Tsel.iloc[0]['ref']['latitude']
        DD_pos_end.loc[k, [TF+'_lon', TF+'_lat']]    =  Tsel.iloc[-1]['ref']['longitude'] ,Tsel.iloc[-1]['ref']['latitude']
        print('result-------', k, TF, data_density, slope_test)


# %%
DD_pos_start

TF1_start = DD_pos_start[ ['TF1_lon', 'TF1_lat'] ].iloc[  abs(DD_pos_start['TF1_lat']).astype('float').argmin() ]
TF2_start = DD_pos_start[ ['TF2_lon', 'TF2_lat'] ].iloc[  abs(DD_pos_start['TF2_lat']).astype('float').argmin() ]

DD_pos_end

TF1_end = DD_pos_end[ ['TF1_lon', 'TF1_lat'] ].iloc[  abs(DD_pos_start['TF1_lat']).astype('float').argmax() ]
TF2_end = DD_pos_end[ ['TF2_lon', 'TF2_lat'] ].iloc[  abs(DD_pos_start['TF2_lat']).astype('float').argmax() ]

DD_pos_end

# %%
# Test if 1st slope segment is negative. There might be better way to test for waves in the data
DD_slope_mask = DD_slope <0

# if there is at leat one slope pair download data, otherwise write files and exit
if (DD_slope_mask.sum() > 1).sum() > 0:
    print('download data')

else:
    print('no suffcient data, quit()')
    ll_name   = ATlevel+'_stats_'+track_name+'_fail'
    DD_merge  = pd.concat({'density_Nperm':DD_data , 'slopes':DD_slope}, axis=1)
    DD_merge.to_html(save_path+ll_name+'.html')
    #DD_merge.columns = ['-'.join(col).strip() for col in DD_merge.columns.values]
    MT.save_pandas_table({'T':DD_merge},ll_name, save_path)
    #DD_merge.columns = ['-'.join(col).strip() for col in DD_merge.columns.values]
    #MT.json_save(name=ll_name, path=save_path, data= DD_merge.where(pd.notnull(DD_merge), 0).T.to_dict())
    MT.json_save(name='A01b_success_'+track_name, path=save_path, data= DD_slope.where(pd.notnull(DD_slope), 0).to_dict())
    exit()



# distill regions of interest
region_list   = list()
#DD_region[~DD_slope_mask] = (np.nan)
for i in DD_region.to_numpy().flatten()[DD_slope_mask.to_numpy().flatten()]:

    if hemis == 'SH':
        region_list.append(i)
    else:
        [region_list.append(ii) for ii in i]

region_list = list(set(region_list))
region_list = [str(int(i)).zfill(2) for i in region_list]
print('region(s) ', region_list)

# %%
class case_ID(object):
    """docstring for case_ID"""
    def __init__(self, track_name):
        import re
        super(case_ID, self).__init__()

        #track_name_pattern = r'(\d{4})(\d{2})(\d{2})(\d{2})(\d{2})(\d{2})_(\d{4})(\d{2})(\d{2})_(\d{3})_(\d{2})'
        track_name_pattern = r'(\D{2}|\d{2})_?(\d{4})(\d{2})(\d{2})(\d{2})?(\d{2})?(\d{2})?_(\d{4})(\d{2})(\d{2})_?(\d{3})?_?(\d{2})?'
        case_ID_pattern = r'(\d{4})(\d{2})(\d{2})_(\d{4})(\d{2})(\d{2})'

        track_name_rx = re.compile(track_name_pattern)
        self.hemis,self.YY,self.MM,self.DD,self.HH,self.MN,self.SS,self.TRK,self.CYC,self.GRN,self.RL,self.VRS = track_name_rx.findall(track_name).pop()

        if self.hemis == '01':
            self.hemis = 'NH'
        elif self.hemis == '02':
            self.hemis = 'SH'
        else:
            self.hemis = self.hemis
        #self.hemis = hemis
        self.set()
        self.track_name_init = track_name

    def set(self):
        block1 = (self.YY,self.MM,self.DD)
        block2 = (self.TRK,self.CYC,self.GRN)

        self.ID = self.hemis+'_'+''.join(block1) +'_'+ ''.join(block2)
        return self.ID

    def set_ATL03_trackname(self):

        block1 = (self.YY,self.MM,self.DD)
        block1b = (self.HH,self.MN,self.SS)
        block2 = (self.TRK,self.CYC,self.GRN)
        if self.RL is '':
            raise ValueError("RL not set")
        if self.VRS is '':
            raise ValueError("VRS not set")

        block3 = (self.RL,self.VRS)

        self.ID = ''.join(block1) +''.join(block1b) +'_'+ ''.join(block2) +'_'+ '_'.join(block3)
        return self.ID

    def set_ATL10_trackname(self):

        block1 = (self.YY,self.MM,self.DD)
        block1b = (self.HH,self.MN,self.SS)
        block2 = (self.TRK,self.CYC, '01') # granule is alwasy '01' for ATL10
        if self.RL is '':
            raise ValueError("RL not set")
        if self.VRS is '':
            raise ValueError("VRS not set")

        block3 = (self.RL,self.VRS)

        if self.hemis == 'NH':
            hemis = '01'
        elif self.hemis == 'SH':
            hemis = '02'
        else:
            hemis = self.hemis

        self.ID = hemis+'_'+''.join(block1) +''.join(block1b) +'_'+ ''.join(block2) +'_'+ '_'.join(block3)
        return self.ID


ID1= 'ATL03-NH_20190102_00780200'
CID = case_ID(ID1)
CID.track_name_init
CID.set()
CID.GRN='10'
CID.set()
CID.RL, CID.VRS='001', '01'
CID.set_ATL03_trackname()
CID.set_ATL10_trackname()

CID = case_ID('ATL10-01_'+track_name)
CID.track_name_init
CID.set()
CID.GRN='40'
CID.set()
CID.set_ATL03_trackname()
CID.set_ATL10_trackname()

# %%
# generate file names for ATL03

DD_list = list()
for TF,TF_polward in zip(['TF1', 'TF2'], [TF1_poleward, TF2_poleward]):
    iregion = DD_region[TF][DD_slope_mask[TF]]
    if len(iregion) !=0:
        iregion2 =list(iregion[0])
        print(iregion2)
        # create track dict
        CID = case_ID(hemis+'_'+track_name)
        CID.GRN = iregion2[0]
        DD= {'case_ID':  CID.set() ,  'tracks' : {} }

        ATL03_list= list()
        for i in iregion2:
            CID = case_ID(hemis+'_'+track_name)
            CID.GRN = i
            # print(CID.set() )
            # print(CID.set_ATL03_trackname() )
            ATL03_list.append(CID.set_ATL03_trackname())

        DD['tracks']['ATL03']   = ['ATL03_'+i for i in ATL03_list]
        DD['tracks']['ATL10']   = 'ATL10-' +CID.set_ATL10_trackname()

        # add other pars:
        DD['pars'] ={'poleward':TF_polward, }
        print(DD)
        DD_list.append(DD)


DD_list



DD


# print results and write files to exit
print('data density N/meter')
print(DD_data)

print('slopes')
print(DD_slope)

#DD_slope.to_html()
for ll in ATL03_list:
    ll_name = 'ATL03_stats_'+ll
    DD_merge = pd.concat({'density_Nperm':DD_data , 'slopes':DD_slope}, axis=1)
    DD_merge.to_html(save_path+ll_name+'.html')
    #DD_merge.columns = ['-'.join(col).strip() for col in DD_merge.columns.values]
    MT.save_pandas_table({'T':DD_merge},ll_name, save_path)
    #MT.json_save(name=ll_name, path=save_path, data= DD_merge.where(pd.notnull(DD_merge), 0).T.to_dict())
    #DD_merge.to_json(save_path+ll_name+'.json', orient="records", lines=True)

MT.json_save(name='A01b_success_'+track_name, path=save_path, data= DD_slope.where(pd.notnull(DD_slope), 0).to_dict())
