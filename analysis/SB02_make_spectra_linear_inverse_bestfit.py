import os, sys
#execfile(os.environ['PYTHONSTARTUP'])

"""
This file open a ICEsat2 track applied filters and corections and returns smoothed photon heights on a regular grid in an .nc file.
This is python 3
"""

exec(open(os.environ['PYTHONSTARTUP']).read())
exec(open(STARTUP_2021_IceSAT2).read())

base_path='/Users/Shared/Projects/2021_IceSAT2_tracks/'
sys.path.append(base_path +'modules/')
sys.path.append(base_path +'modules/ICEsat2_SI_tools/')

import matplotlib.pyplot as plt
%matplotlib inline

#import m_general as M
#import m_tools as MT
import numpy as np

import m_general_ph3 as M

import netCDF4
import datetime
import os
from netCDF4 import Dataset
import xarray as xr
import pandas as pd

import ICEsat2_SI_tools.convert_GPS_time as cGPS
import h5py
import ICEsat2_SI_tools.io as io
import ICEsat2_SI_tools.spectral_estimates as spec

import imp

import lmfit
#import s3fs
# %%

# Python reader based on Pandas. Other reader examples available in readers.py
plot_path = mconfig['paths']['plot'] + '/tests/'


track_name= 'ATL03_20190515060436_07170312_002_02'
load_path   = base_path + 'data/data1/'

# test which beams exist:
all_beams = ['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']
low_beams = ['gt1l',  'gt2l',  'gt3l']
high_beams = ['gt1r',  'gt2r',  'gt3r']

# %%

#Gall= xr.open_dataset(load_path+'/'+track_name +'_filtered_photon_heights.nc')
imp.reload(io)
Gfilt = io.load_pandas_table_dict(track_name + '_B01_corrected', load_path)
Gd = io.load_pandas_table_dict(track_name + '_B01_binned' , load_path)

# %%
Gi =Gd[high_beams[0]]
Gi = Gi[~np.isnan(Gi['heights_c_weighted_mean'])].sort_values('dist')
#np.isnan(Gi['dist']).sum()
Gi = Gi[(Gi['dist'] >= 140000) & (Gi['dist'] <= 170001)]

# %% test and compare to FFT
from scipy.signal import detrend
y =detrend(np.array(Gi['heights_c_weighted_mean']) )

y_gap= np.copy(y)
y_gap[500:2000] = np.nan
nan_mask =np.isnan(y_gap)
data_filled = np.copy(y)
data_filled[nan_mask] = 0

nan_mask.sum()/y_gap.size

plt.plot( Gi['dist'], y, '-')
plt.plot( Gi['dist'], y_gap, '-')

# %%
x= Gi['dist'] - Gi['dist'].iloc[0]
dx = np.median(np.diff(x))
plt.plot( x, y_gap)


y_gap_fill=np.copy(y_gap)
y_gap_fill[nan_mask]= 0


k_fft = np.fft.rfftfreq(x.size, d=dx)
guess1=  M.runningmean(abs(np.fft.rfft(y_gap_fill) ), 30)
guess1[np.isnan(guess1)] =0
plt.plot(k_fft, guess1)
k_max = k_fft[ guess1.argmax() ]
f_max= np.sqrt(9.81 * k_max) / 2/np.pi

import JONSWAP_gamma as spectal_models
spectal_models.pierson_moskowitz_default(f_max, 10)

U= 10
imp.reload(spectal_models)
spectal_models.JONSWAP_default_alt(f_max,U )

# choose wavenumber range
dlambda = (x.size *dx/2)

T_min=5
lambda_min = 9.81 * T_min**2/ (2 *np.pi)

k =np.arange(0, 1/lambda_min,  1/dlambda) * 2*np.pi
#k = k[1:201]

#k_full = np.copy(np.fft.rfftfreq(x.size, d=dx) * 2* np.pi)
#1/np.diff(k_full/2/np.pi).mean()
#k = k_full[1:201]
MM = k.size
NN = x.size

x_gap = x[~nan_mask]

#plt.plot( x_gap , np.cos(np.outer(x_gap, k[0:10])), '.' )
G = np.vstack([ np.cos(np.outer(x_gap, k)).T ,  np.sin(np.outer(x_gap, k)).T ] ).T
G_full = np.vstack([ np.cos(np.outer(x, k)).T ,  np.sin(np.outer(x, k)).T ] ).T

# Model Penalty
# regullization . 1/weight

#plt.plot( np.arange(MM, 0, -1)/MM )
#weights = np.concatenate([ np.arange(MM, 0, -1)/MM ,  np.arange(MM, 0, -1)/MM  ])
weights = np.concatenate([ np.arange(MM, 0, -1)/MM ,  np.arange(MM, 0, -1)/MM  ]) *0+1
#weights = np.concatenate([ np.exp(-4 *np.arange(0, MM)/MM),  np.exp(-4 *np.arange(0, MM)/MM) ])
P = np.diag(1/weights)
P.shape

# uncertainty of the data
R = np.diag( np.array(Gi['heights_c_std'])[~nan_mask] /1.5 )
#np.array(Gi['heights_c_std'])[~nan_mask].std()

from numpy import linalg
inv = linalg.inv
#inv(  (G.T @ inv(R) @ G  ) + inv(P) ).shape
b_hat = inv(  (G.T @ inv(R) @ G  ) + inv(P) ) @ G.T @ inv(R) @ y[~nan_mask]
# %%

M.figure_axis_xy(10,5,  view_scale =0.6)

plt.subplot(2,1, 1)

plt.plot( x, y_gap, '-k', linewidth= 2)
#plt.plot( x, y, '-k', linewidth= 2)

plt.plot( x_gap , (b_hat * G).sum(1), 'r-', markersize= 1 )
plt.plot( x[nan_mask] , (b_hat * G_full).sum(1)[nan_mask], '-', color='orange', linewidth = 0.5 )
plt.xlim(1000,8000)
plt.ylim(-1, 1)

# %
plt.subplot(2,1, 2)

#plt.plot( x, (y - (b_hat * G).sum(1)) , '.k', markersize= 0.8 )
plt.plot( x_gap, y_gap[~nan_mask] - (b_hat * G).sum(1) , '.k', markersize= 0.8 )
rr= y_gap[~nan_mask] - (b_hat * G).sum(1)

y_gap[~nan_mask].var() - (b_hat * G).sum(1).var() - rr.var()

print( ( rr.std() /np.diag(R).std() ).mean() )

plt.xlim(1000,8000)
plt.ylim(-1, 1)


# %%

Z = b_hat[0:MM] - b_hat[MM:] *1j
plt.plot(k, abs(Z) )

# %%
