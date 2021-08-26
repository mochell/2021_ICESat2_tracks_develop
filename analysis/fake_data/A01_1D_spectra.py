import os, sys
#execfile(os.environ['PYTHONSTARTUP'])

"""
This file open a ICEsat2 track applied filters and corections and returns smoothed photon heights on a regular grid in an .nc file.
This is python 3
"""
# exec(open(os.environ['PYTHONSTARTUP']).read())
# exec(open(STARTUP_2019_DP).read())

base_path='/Users/Shared/Projects/2021_IceSAT2_tracks/'
sys.path.append(base_path +'modules/')
sys.path.append(base_path +'modules/ICEsat2_SI_tools/')

import matplotlib.pyplot as plt
%matplotlib inline

#import m_general as M
#import m_tools as MT
import numpy as np

import m_general_ph3 as M
import os
from netCDF4 import Dataset
import xarray as xr
import pandas as pd
import h5py

import ICEsat2_SI_tools.convert_GPS_time as cGPS
import ICEsat2_SI_tools.io as io

import imp
import m_spectrum_ph3 as spec
import JONSWAP_gamma
# %%

f = np.arange(0.001, 0.2,0.001)
spec_power = JONSWAP_gamma.JONSWAP_default(f, 2e6, 10)
plt.plot(f, spec_power)

amps = (spec_power * np.gradient(f)) **0.5

2/f[amps.argmax()]

t = np.arange(0, 1000, 0.001)
tt, ww = np.meshgrid(2* np.pi * f, t)
phi = np.random.random(len(amps))*2*np.pi

instance = amps* np.cos( ww * tt + phi )

# %%
M.figure_axis_xy(6, 2, view_scale=0.9)
plt.plot(t, instance.mean(1) )
