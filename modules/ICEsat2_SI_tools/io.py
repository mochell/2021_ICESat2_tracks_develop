

def save_pandas_table(table_dict, name , save_path):
    import os
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    import warnings
    from pandas import HDFStore
    from pandas.io.pytables import PerformanceWarning
    warnings.filterwarnings('ignore',category=PerformanceWarning)

    with HDFStore(save_path+'/'+name+'.h5') as store:
        for name,table in table_dict.items():
                store[name]=table
def load_pandas_table_dict(name , save_path):
    import warnings
    from pandas import HDFStore
    from pandas.io.pytables import PerformanceWarning
    warnings.filterwarnings('ignore',category=PerformanceWarning)

    return_dict=dict()
    with HDFStore(save_path+'/'+name+'.h5') as store:
        #print(store)
        #print(store.keys())
        for k in store.keys():
            return_dict[k[1:]]=store.get(k)

    return return_dict


def get_time_for_track(delta_time, atlas_epoch):
    "returns pandas dataframe"
    import pandas as pd
    import convert_GPS_time as cGPS
    # Conversion of delta_time to a calendar date
    temp = cGPS.convert_GPS_time(atlas_epoch[0] + delta_time, OFFSET=0.0)

    year = temp['year'][:].astype('int')
    month = temp['month'][:].astype('int')
    day = temp['day'][:].astype('int')
    hour = temp['hour'][:].astype('int')
    minute = temp['minute'][:].astype('int')
    second = temp['second'][:].astype('int')

    return pd.DataFrame({'year':year, 'month':month, 'day':day, 'hour':hour, 'second':second})

def getATL03_beam(fileT, numpy=False, beam='gt1l', maxElev=1e6):
    """
    returns 'beam' from fileT as pandas table.
    fillT   path of file
    numpy=  (False), or True. if True the method returns a list of numpy instances,
            if False it returns a pandas table
    beam    key of the iceSAT2 beam.
    """
    # Add in a proper description of the function here
    import h5py
    import convert_GPS_time as cGPS
    import pandas as pd
    # Open the file
    ATL03       =   h5py.File(fileT, 'r')
    lons        =   ATL03[beam+'/heights/lon_ph'][:]
    lats        =   ATL03[beam+'/heights/lat_ph'][:]

    # Along track distance from equator i think.
    along_track_distance=ATL03[beam+'/heights/dist_ph_along'][:]
    #dem_h = ATL03[beam+'/geophys_corr/dem_h'][:]
    #delta_time_dem_h = ATL03[beam+'/geophys_corr/delta_time'][:]
    #segment_dist_x=ATL03[beam+'/geolocation/segment_dist_x'][:]

    #  Nathan says it's the number of seconds since the GPS epoch on midnight Jan. 6, 1980
    delta_time  =   ATL03[beam+'/heights/delta_time'][:]
    #podppd_flag=ATL03[beam+'/geolocation/podppd_flag'][:]

    # #Add this value to delta time parameters to compute full gps_seconds
    atlas_epoch =   ATL03['/ancillary_data/atlas_sdp_gps_epoch'][:]

    # Conversion of delta_time to a calendar date
    temp = cGPS.convert_GPS_time(atlas_epoch[0] + delta_time, OFFSET=0.0)

    # Express delta_time relative to start time of granule
    delta_time_granule=delta_time-delta_time[0]

    year = temp['year'][:].astype('int')
    month = temp['month'][:].astype('int')
    day = temp['day'][:].astype('int')
    hour = temp['hour'][:].astype('int')
    minute = temp['minute'][:].astype('int')
    second = temp['second'][:].astype('int')

    # Primary variables of interest

    # Photon height
    heights=ATL03[beam+'/heights/h_ph'][:]
    #print(heights.shape)

    # Flag for signal confidence
    # column index:  0=Land; 1=Ocean; 2=SeaIce; 3=LandIce; 4=InlandWater
    # values:
        #-- -1: Events not associated with a specific surface type
        #--  0: noise
        #--  1: buffer but algorithm classifies as background
        #--  2: low
        #--  3: medium
        #--  4: high
    signal_confidence=ATL03[beam+'/heights/signal_conf_ph'][:,2]
    #print(signal_confidence.shape)



    # Add photon rate and background rate to the reader here
    ATL03.close()

    if numpy == True:
        # list the variables you want to output here..
        return along_track_dist, elev

    else:
        dF = pd.DataFrame({'heights':heights, 'lons':lons, 'lats':lats, 'signal_confidence':signal_confidence,
                       'delta_time_granule':delta_time_granule,'delta_time':delta_time, 'along_track_distance':along_track_distance,
                        'year':year, 'month':month, 'day':day, 'hour':hour,'minute':minute , 'second':second})

        # Filter out high elevation values
        #dF = dF[(dF['signal_confidence']>2)]
        # Reset row indexing
        #dF=dF.reset_index(drop=True)
        return dF

def getATL03_height_correction(fileT, beam='gt1r'):
    """
    This method returns relevant data for wave estimates from ALT 07 tracks.
    returns: Pandas data frame
    """
    # Add in a proper description of the function here

    import h5py
    import pandas as pd
    # Open the file
    ATL03 = h5py.File(fileT, 'r')

    ### bulk positions and statistics
    vars_bulk = [
            'delta_time', # referenc time since equator crossing
            'dem_h', # best giod approxiamtion
            ]

    D_bulk= dict()
    for var in vars_bulk:
        D_bulk[var] = ATL03[beam+'/geophys_corr/'+var][:]
    dF_bulk = pd.DataFrame(D_bulk)

    ATL03.close()

    return dF_bulk


def getATL07_beam(fileT, beam='gt1r', maxElev=1e6):
    """
    This method returns relevant data for wave estimates from ALT 07 tracks.
    returns: Pandas data frame
    """
    # Add in a proper description of the function here

    import h5py
    import pandas as pd
    # Open the file
    ATL07 = h5py.File(fileT, 'r')

    ### bulk positions and statistics
    vars_bulk = [
            'longitude',
            'latitude',
            'height_segment_id',# Height segment ID (10 km segments)
            'seg_dist_x' # Along track distance from the equator crossing to the segment center.
            ]

    D_bulk= dict()
    for var in vars_bulk:
        D_bulk[var] = ATL07[beam+'/sea_ice_segments/'+var]
    dF_bulk = pd.DataFrame(D_bulk)

    #  Nathan says it's the number of seconds since the GPS epoch on midnight Jan. 6, 1980
    delta_time=ATL07[beam+'/sea_ice_segments/delta_time'][:]
    # #Add this value to delta time parameters to compute full gps_seconds
    atlas_epoch=ATL07['/ancillary_data/atlas_sdp_gps_epoch'][:]
    dF_time = get_time_for_track(delta_time,atlas_epoch)
    dF_time['delta_time'] = delta_time
    ### Primary variables of interest
    vars = [
            'across_track_distance', #Across track distance of photons averaged over the sea ice height segment.
            'height_segment_asr_calc', #Computed apparent surface reflectance for the sea ice segment.
            'height_segment_confidence',# # Height segment confidence flag
            'height_segment_fit_quality_flag', # Flag Values: ['-1', '1', '2', '3', '4', '5']
                                            #Flag Meanings: ['invalid', 'best', 'high', 'med', 'low', 'poor']
            'height_segment_height',     # Beam segment height
            'height_segment_length_seg',  # Along track length of segment
            'height_segment_ssh_flag', #Flag for potential leads, 0=sea ice, 1 = sea surface
            'height_segment_surface_error_est', #Error estimate of the surface height
            'height_segment_type',# Flag Values: ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
                                  # Flag Meanings: ['cloud_covered', 'other', 'specular_lead_low_w_bkg', 'specular_lead_low', 'specular_lead_high_w_bkg', 'specular_lead_high', 'dark_lead_smooth_w_bkg', 'dark_lead_smooth'
            'height_segment_w_gaussian', # Width of Gaussian fit
            'height_segment_quality', # Height quality flag, 1 for good fit, 0 for bad
            ]
    D_heights=dict()
    for var in vars:
        D_heights[var] = ATL07[beam+'/sea_ice_segments/heights/' +var][:]
    dF_heights = pd.DataFrame(D_heights)


    vars_env = {
            'mss':'geophysical/height_segment_mss', # Mean sea surface height above WGS-84 reference ellipsoid (range: -105 to 87m), based on the DTU13 model.
            't2m':'geophysical/height_segment_t2m',#Temperature at 2m above the displacement height (K)
            'u2m':'geophysical/height_segment_u2m',#Eastward wind at 2m above the displacement height (m/s-1)
            'v2m':'geophysical/height_segment_v2m',#Northward wind at 2m above the displacement height (m/s-1)
            'n_photons_actual':'stats/n_photons_actual', # Number of photons gathered
            'photon_rate':'stats/photon_rate', #photon_rate
    }

    D_env=dict()
    for var,I in vars_env.items():
        D_env[  var] = ATL07[beam+'/sea_ice_segments/' +I][:]
    dF_env = pd.DataFrame(D_env)


    #Df = pd.concat({k: pd.DataFrame(v).T for k, v in data.items()}, axis=0)
    DF = pd.concat({ 'time': dF_time,  'ref': dF_bulk, 'heights': dF_heights, 'env': dF_env }, axis=1)

    ATL07.close()

    # Filter out high elevation values
    DF = DF[(DF['heights']['height_segment_height']<maxElev)]
    # Reset row indexing
    DF=DF.reset_index(drop=True)
    return DF

def getATL07_height_corrections(fileT, beam='gt1r'):
    """
    This method returns relevant data for wave estimates from ALT 07 tracks.
    returns: Pandas data frame
    """
    # Add in a proper description of the function here

    import h5py
    import pandas as pd
    # Open the file
    ATL07 = h5py.File(fileT, 'r')

    ### bulk positions and statistics
    vars_bulk = [
            'longitude',
            'latitude',
            'height_segment_id',# Height segment ID (10 km segments)
            'seg_dist_x' # Along track distance from the equator crossing to the segment center.
            ]

    D_bulk= dict()
    for var in vars_bulk:
        D_bulk[var] = ATL07[beam+'/sea_ice_segments/'+var]
    dF_bulk = pd.DataFrame(D_bulk)

    #  Nathan says it's the number of seconds since the GPS epoch on midnight Jan. 6, 1980
    delta_time=ATL07[beam+'/sea_ice_segments/delta_time'][:]
    # #Add this value to delta time parameters to compute full gps_seconds
    atlas_epoch=ATL07['/ancillary_data/atlas_sdp_gps_epoch'][:]
    dF_time = get_time_for_track(delta_time,atlas_epoch)

    ### Primary variables of interest
    vars = [
            'height_segment_dac', #
            'height_segment_ib', #
            'height_segment_lpe',# #
            'height_segment_mss', #
            'height_segment_ocean',
            ]
    D_heights=dict()
    for var in vars:
        D_heights[var] = ATL07[beam+'/sea_ice_segments/geophysical/' +var][:]
    dF_heights = pd.DataFrame(D_heights)


    #Df = pd.concat({k: pd.DataFrame(v).T for k, v in data.items()}, axis=0)
    DF = pd.concat({ 'time': dF_time,  'ref': dF_bulk, 'corrections': dF_heights, }, axis=1)

    ATL07.close()
    # Filter out high elevation values
    # Reset row indexing
    DF=DF.reset_index(drop=True)
    return DF
