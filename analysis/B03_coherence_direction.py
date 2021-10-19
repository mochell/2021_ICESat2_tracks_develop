
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
#import s3fs
# %%
track_name, batch_key, test_flag = io.init_from_input(sys.argv) # loads standard experiment
#track_name, batch_key, test_flag = '20190605061807_10380310_004_01', 'SH_batch01', False
track_name, batch_key, test_flag = '20190601094826_09790312_004_01', 'SH_batch01', False
#track_name, batch_key, test_flag = '20190207111114_06260210_004_01', 'SH_batch02', False


#print(track_name, batch_key, test_flag)
hemis, batch = batch_key.split('_')
#track_name= '20190605061807_10380310_004_01'
ATlevel= 'ATL03'

save_path   = mconfig['paths']['work'] + '/B03_spectra_'+hemis+'/'
save_name   = 'B03_'+track_name

plot_path   = mconfig['paths']['plot'] + '/'+hemis+'/'+batch_key+'/' + track_name + '/B_spectra/'
MT.mkdirs_r(plot_path)
MT.mkdirs_r(save_path)
bad_track_path =mconfig['paths']['work'] +'bad_tracks/'+ batch_key+'/'
# %%

all_beams   = mconfig['beams']['all_beams']
high_beams  = mconfig['beams']['high_beams']
low_beams   = mconfig['beams']['low_beams']
#Gfilt   = io.load_pandas_table_dict(track_name + '_B01_regridded', load_path) # rhis is the rar photon data

load_path   = mconfig['paths']['work'] +'/B01_regrid_'+hemis+'/'
Gd      = io.load_pandas_table_dict(track_name + '_B01_binned' , load_path)  #
load_path   = mconfig['paths']['work'] + '/B02_spectra_'+hemis+'/'
Gpars   = io.load_pandas_table_dict('B02_'+ track_name + '_params' , load_path)  #
Gspec   = xr.open_dataset(load_path + 'B02_'+ track_name + '_LS.nc' )  #
Gspec['Y_model_hat'] = Gspec.Y_model_hat_real + Gspec.Y_model_hat_imag *1j
Gspec = Gspec.drop('Y_model_hat_real').drop('Y_model_hat_imag')

dk = Gspec.k.diff('k').mean().data
Lpoints = Gspec.Lpoints

# %%

k_lim= 0.02
A, B = Gspec.sel(beam= 'gt2r', k=slice(0, k_lim)).Y_model_hat , Gspec.sel(beam= 'gt2l', k=slice(0, k_lim)).Y_model_hat

#(abs(B) - abs(A)).plot()
S_aa = (A*A.conj()).real
S_bb = (B*B.conj()).real
co_spec = (A.conj() *B).rolling(x=5, k=3, center=True, min_periods=2).mean()

abs(co_spec).plot( levels= np.arange(0, 230,10))

phase = np.arctan2( -co_spec.imag, co_spec.real )

phase.plot()
#plt.plot(phase.k, phase, 'k', linewidth=0.8, alpha= 0.6)
phase.sel(x=slice(500000, 600000)).mean('x').rolling(k=20, center=True, min_periods=2).mean().plot()#.hist(bins= 40)
phase.sel(x=slice(500000, 600000)).mean('x').rolling(k=20, center=True, min_periods=2).std().plot()#.hist(bins= 40)

# %%
phase.mean('x').plot()#.hist(bins= 40)

phase.where(phase > np.pi/2, phase + np.pi).plot()


C_ab= abs(co_spec)**2 #/ S_aa/ S_bb
np.log(C_ab).sel(k=slice(0, 0.03)).plot()

C_ab= abs(co_spec)**2 / S_aa/ S_bb
(C_ab).sel(k=slice(0, 0.03)).plot(levels = np.arange(0.5, 1+0.1, 0.05))
#(C_ab.where(~np.isnan(C_ab), 0)).sel(k=slice(0, 0.03)).plot.hist(bins= 30)

x_nan_mask = (~np.isnan(C_ab)).sum('k') == 0

C_ab[:, ~x_nan_mask].min()

phase = np.arctan2( -co_spec.imag, co_spec.real )

plt.plot( phase )
plt.plot(S_bb, 'k', linewidth= 0.5, alpha =0.5)

np.nanmin(C_ab.imag)
np.nanmax(C_ab.imag)
Gspec.beam

# %%
Ai = A.sel(x= 500000, method='nearest')
Bi = B.sel(x= 500000, method='nearest')

M.figure_axis_xy(10, 2.5, view_scale= 0.7)
plt.plot(np.fft.irfft(Ai) )
plt.plot(np.fft.irfft(Bi) )

# %%


co_spec=Ai.imag *Bi.real
C_ab= abs(co_spec)**2 / ( (Ai * Ai.conj()).real * (Bi * Bi.conj()).real  )
plt.plot(C_ab.k, C_ab)
#plt.semilogx(C_ab.k, abs(Ai.imag *Bi.real)**2)


phase = np.arctan2( -co_spec.imag, co_spec.real )
phase[phase <- np.pi/2] = phase[phase < - np.pi/2] + np.pi
phase.plot.hist()
plt.plot(phase.k, phase)
# %%


def get_coherency(d1, d2, Lint, dt=None):
    """
    returns the squared coherence Cxy and its phase
    """
    from scipy import signal

    if dt is None:
        dt =np.diff((dd.time - dd.time[0]).data.astype('m8[D]')).mean().astype('int')
    else:
        dt= dt

    coh = signal.coherence( d1, d2, fs=1/dt, window ='hann',  nperseg=Lint,  noverlap=None, detrend='linear')

    co_spec = signal.csd( d1, d2, fs=1/dt, window ='hann',  nperseg=Lint,noverlap=None, detrend='linear')
    phase = np.arctan2( -co_spec[1].imag, co_spec[1].real )

    #Spec = m_spec.Spectrum(dd, dt=dt, verbose=False)
    #co_spec = signal.csd( d1, d2, fs=1/dt, window ='hann',  nperseg=Lint)
    #np.arctan2(-)
    return {'f' : coh[0], 'coh2' : coh[1], 'coh' : np.sqrt(abs(coh[1])), 'cospec' :co_spec , 'phase':phase}
