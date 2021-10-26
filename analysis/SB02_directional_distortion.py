
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


# %%

f = np.arange(1/1000, 1/5, 1/500)
import JONSWAP_gamma as spectal_models


col.colormaps2(21)
font_for_pres()
F = M.figure_axis_xy(5,3, view_scale= 0.7)

U= 20 # results are incensitive to U
f_max = 1/20
gamma= 1
Jswap= spectal_models.JONSWAP_default_alt(f, f_max, 20 ,  gamma=gamma)

import itertools

clist = itertools.cycle([col.cascade1, col.rascade1, col.cascade2, col.rascade2])

#plt.plot(1/f, Jswap, 'k', label='true', zorder=12, linewidth = 0.5)
for f_max in np.arange(1/25, 1/10, 1/100):
#for U in np.arange(1, 30, 10):

    fold= f_max
    cc = next(clist)
    for alpha in np.arange(-82.5, 85, 7.5):
        #alpha = 80
        f_prime= f * np.sqrt( np.cos(np.pi *alpha/ 180) )
        Jswap= spectal_models.JONSWAP_default_alt(f, f_max, U ,  gamma=gamma)

        # f_max_prime = f_max * np.sqrt( np.cos(np.pi *alpha/ 180) )
        #
        # Jswap_prime= spectal_models.JONSWAP_default_alt(f_prime, f_max, 20 ,  gamma=gamma)

        if alpha == 0:
            lstring= '$T_p$=' + str(np.round(1/f_max, 1)) +'s'
            plt.plot(1/f_prime, Jswap, c=cc, label=lstring, linewidth = 1)
        else:
            plt.plot(1/f_prime, Jswap, c=cc, linewidth = 0.5)

#F.ax.set_yscale('log')
plt.title('Spectal Distortion with oberservation angle ($\\alpha \pm85 ^\circ$)')
F.ax.set_xscale('log')
plt.xlabel("observed Period (T')")
plt.xlim(0, 250)
plt.legend()

# %%


f = np.arange(1/1000, 1/5, 1/500)
import JONSWAP_gamma as spectal_models


col.colormaps2(21)
font_for_pres()
F = M.figure_axis_xy(5,3, view_scale= 0.7)

U= 20 # results are incensitive to U
f_max = 1/20
gamma= 1
Jswap= spectal_models.JONSWAP_default_alt(f, f_max, 20 ,  gamma=gamma)

import itertools

clist = itertools.cycle([col.cascade1, col.rascade1, col.cascade2, col.rascade2])

#plt.plot(1/f, Jswap, 'k', label='true', zorder=12, linewidth = 0.5)
for f_max in np.arange(1/25, 1/10, 1/100):
#for U in np.arange(1, 30, 10):

    fold= f_max
    cc = next(clist)
    for alpha in np.arange(-82.5, 85, 7.5):
        #alpha = 80
        f_prime= f * np.sqrt( np.cos(np.pi *alpha/ 180) )
        Jswap= spectal_models.JONSWAP_default_alt(f, f_max, U ,  gamma=gamma)

        # f_max_prime = f_max * np.sqrt( np.cos(np.pi *alpha/ 180) )
        #
        # Jswap_prime= spectal_models.JONSWAP_default_alt(f_prime, f_max, 20 ,  gamma=gamma)

        if alpha == 0:
            lstring= '$T_p$=' + str(np.round(1/f_max, 1)) +'s'
            plt.plot(1/f_prime, Jswap, c=cc, label=lstring, linewidth = 1)
        else:
            plt.plot(1/f_prime, Jswap, c=cc, linewidth = 0.5)

#F.ax.set_yscale('log')
plt.title('Spectal Distortion with oberservation angle ($\\alpha \pm85 ^\circ$)')
F.ax.set_xscale('log')
plt.xlabel("observed Period (T')")
plt.xlim(0, 250)
plt.legend()
