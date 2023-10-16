
# %%
"""
This file open a ICEsat2 track applied filters and corections and returns smoothed photon heights on a regular grid in an .nc file.
This is python 3.11
"""
import os, sys
exec(open(os.environ["PYTHONSTARTUP"]).read())
exec(open(STARTUP_2021_IceSAT2).read())

import geopandas as gpd

from sliderule import icesat2
from sliderule import sliderule

import shapely
from ipyleaflet import basemaps, Map, GeoData

import ICEsat2_SI_tools.sliderule_converter_tools as sct
import ICEsat2_SI_tools.io as io

import spicke_remover

import h5py, imp, copy

xr.set_options(display_style='text')

# %matplotlib inline
# %matplotlib widget


# %% Select region and retrive batch of tracks

# define file with ID:
track_name, batch_key , ID_flag = '20190219073735_08070210_005_01', 'SH_testSLsinglefile' , False

plot_flag = True


save_path  = mconfig['paths']['work'] +'/'+batch_key+'/B01_regrid/'
MT.mkdirs_r(save_path)

#ID, _, _, _ = io.init_data(track_name, batch_key, True, mconfig['paths']['work'],  )
ATL03_track_name = 'ATL03_'+track_name+'.h5'
#track_name = ID['tracks']['ATL03'][0] +'.h5'

# %% Configure SL Session #
icesat2.init("slideruleearth.io", True) 
asset = 'nsidc-s3'

# %% plot the ground tracks in geographic location
# Generate ATL06-type segments using the ATL03-native photon classification
# Use the ocean classification for photons with a confidence parmeter to 2 or higher (low confidence or better)

params={'srt': 1,  # Ocean classification
 'len': 25,        # 10-meter segments
 'ats':3,          # require that each segment contain photons separated by at least 5 m
 'res':10,         # return one photon every 10 m
 'track': 0,       # return all ground tracks
 'pass_invalid': True,   
 'cnf': 2,         # require classification confidence of 2 or more
 #'iterations':10,  # iterate the fit
}

# YAPC alternatibe
params_yapc={'srt': 1,
 'len': 20,
 'ats':3,
 'res':10,
 'track': 0,
 'pass_invalid': True,
 'cnf': -2,
 'iterations':10,}
#     #   "yapc": dict(knn=0, win_h=6, win_x=11, min_ph=4, score=100),  # use the YAPC photon classifier; these are the recommended parameters, but the results might be more specific with a smaller win_h value, or a higher score cutoff
#   "yapc": dict(knn=0, win_h=3, win_x=11, min_ph=4, score=50),  # use the YAPC photon classifier; these are the recommended parameters, but the results might be more specific with a smaller win_h value, or a higher score cutoff

gdf = icesat2.atl06p(params, asset="nsidc-s3", resources=[ATL03_track_name])

gdf.plot()
# %% main routine for defining the x coordinate and sacing table data

def make_B01_dict(table_data, split_by_beam=True, to_hdf5=False):
    """
    converts a GeoDataFrame from Sliderule to GeoDataFrames for each beam witht the correct columns and names
    inputs:
        table_data: GeoDataFrame with the data
        split_by_beam: True/False. If True the data is split by beam
    returns:
        if split_by_beam:
            table_data: dict of GeoDataFrame with the data for each beam
        else:
            table_data: GeoDataFrame with the data for all beams in one table
    """
    #table_data = copy.copy(B01b_atl06_native)

    table_data.rename(columns= {
                        'n_fit_photons':'N_photos',
                        'w_surface_window_final':'signal_confidence',
                        'across':'y', 
                        #'spot':'beam',
                        }
                        , inplace=True)

    table_data['lons'] = table_data['geometry'].x
    table_data['lats'] = table_data['geometry'].y

    drop_columns = ['cycle','gt','rgt', 'pflags']
    if to_hdf5:
        drop_columns.append('geometry')
    table_data.drop(columns=drop_columns, inplace=True)

    if split_by_beam:
        B01b = dict()
        # this is not tested
        for spot,beam in zip([1,2,3, 4, 5, 6],['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']) : 
            ii = table_data.spot==spot
            B01b[beam] = table_data[ii]
        return B01b
    else:
        return table_data

# define reference point and then define 'x'
table_data = copy.copy(gdf)

# the reference point is defined as the most equatorward point of the polygon. 
# It's distance from the equator is  subtracted from the distance of each photon.
table_data = sct.define_x_coordinate_from_data(table_data)

# fake 'across'
table_data['across'] = table_data['x']*0 +table_data['spot']

# add spike remover

# renames columns and splits beams
Ti = make_B01_dict(table_data, split_by_beam=True, to_hdf5=True) 
#Ti['gt1l'].drop('geometry', axis=1, inplace=True)

segment = track_name.split('_')[1][-2:]
ID_name = sct.create_ID_name(gdf.iloc[0], segment=segment)
print( ID_name )
io.write_track_to_HDF5(Ti, ID_name + '_B01_binned'     , save_path) # regridding heights

print('done')
# %%
