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

# %%




#plt.xlim(140000, 145000)
#plt.ylim(-1, 1)

#np.diff(np.array(Gi['dist']))
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
dx = np.diff(x).mean()
plt.plot( x, y_gap)


k_full = np.copy(np.fft.rfftfreq(x.size, d=dx) * 2* np.pi)
k = k_full[1:201]
MM = k.size
NN = x.size

x_gap = x[~nan_mask]

#plt.plot( x_gap , np.cos(np.outer(x_gap, k[0:10])), '.' )

G = np.vstack([ np.cos(np.outer(x, k)).T ,  np.sin(np.outer(x, k)).T ] ).T
G_full = np.vstack([ np.cos(np.outer(x, k)).T ,  np.sin(np.outer(x, k)).T ] ).T
G.shape

b = np.zeros([2 *MM])
b[4] = 1
b[300] = 1

b.shape
G.shape
y_gap.shape

plt.plot( (G *b).sum(1) )

# Model Penalty
# regullization . 1/weight
# choose simple decay along diagonal

#plt.plot( np.arange(MM, 0, -1)/MM )
#weights = np.concatenate([ np.arange(MM, 0, -1)/MM ,  np.arange(MM, 0, -1)/MM  ])
weights = np.concatenate([ np.arange(MM, 0, -1)/MM ,  np.arange(MM, 0, -1)/MM  ]) *0+1

#weights = np.concatenate([ np.exp(-4 *np.arange(0, MM)/MM),  np.exp(-4 *np.arange(0, MM)/MM) ])
P = np.diag(1/weights)
P.shape

# uncertainty of the data
R = np.diag( np.array(Gi['heights_c_std'])*0+0.0001 )
R.shape

from numpy import linalg
# R = dvar is either a scalar or a vector n_data long
# estimate is:
# phat = (G'*R^-1*G + P^-1)^-1 * G'*R^-1*y  < -------------------------------

inv = linalg.inv
inv(R).shape
G.T.shape

#inv(  (G.T @ inv(R) @ G  ) + inv(P) ).shape
b_hat = inv(  (G.T @ inv(R) @ G  ) + inv(P) ) @ G.T @ inv(R) @ y#[~nan_mask]

plt.plot(b_hat)
# %%

M.figure_axis_xy(10,5,  view_scale =0.6)

plt.subplot(2,1, 1)

#plt.plot( x_gap, y_gap[~nan_mask], '-k', linewidth= 2)
plt.plot( x, y, '-k', linewidth= 2)

plt.plot( x , (b_hat * G).sum(1), 'r-', markersize= 1 )
#plt.plot( x[nan_mask] , (b_hat * G_full).sum(1)[nan_mask], '-', color='orange', linewidth = 0.5 )
plt.xlim(1000,8000)
plt.ylim(-1, 1)

plt.subplot(2,1, 2)

plt.plot( x, (y - (b_hat * G).sum(1)) , '.k', markersize= 0.8 )
#plt.plot( x_gap, y_gap[~nan_mask] - (b_hat * G).sum(1) , '.k', markersize= 0.8 )
#rr= y - (b_hat * G).sum(1)

rr/ 0.1
plt.xlim(1000,8000)
plt.ylim(-1, 1)


# %%
k.shape

sum1 = (b_hat * G).sum(0)
R = sum1[0:MM]
I = sum1[MM:]
C  =(G**2)[:, 0:MM].sum(0)
S  =(G**2)[:, MM:].sum(0)
Spec_power =  4* np.sqrt(  MM * (R**2/ C + I**2/ S) )
plt.plot( k, Spec_power)
#plt.ylim(0,0.06)


# %%
Z_fft = np.fft.rfft(y)


Z = b_hat[0:MM] - b_hat[MM:] *1j
MM
#Z2 = Z / abs(Z)

plt.plot(k_full, abs(Z_fft)/ k_full.size)
plt.plot(k, abs(Z) )
plt.plot( k, Spec_power , 'k')
#np.median( abs(Z_fft) / (abs(Z)*MM) )

# %% compate phase
plt.plot(Z_fft.real[1:201]/ k_full.size, Z.real, '.')
plt.axis('equal')


(Z_fft.real[1:201]/ k_full.size - (Z.real)).mean()
Z_fft.imag[1:201].sum()/ k_full.size
Z.imag.sum()

Z_fft.real[1:201].sum()/ k_full.size
Z.real.sum()

( (Z_fft.imag[1:201]/ k_full.size) - (Z.imag) ).mean()
# plt.plot( (Z_fft.imag) - (Z.imag*MM))
# plt.plot( (Z_fft.real) - (Z.real*MM))
#plt.plot(  ( (Z_fft.imag) - (Z.imag*MM) ) -  (  (Z_fft.real) - (Z.real*MM) ) )
