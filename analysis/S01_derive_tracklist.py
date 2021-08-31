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
import datetime
import os

import ICEsat2_SI_tools.io as io
import ICEsat2_SI_tools.spectral_estimates as spec

import imp

# %%
path = mconfig['paths']['analysis']+'../track_lists/'



with open(path+  'alex_ATL07_filelist.txt', 'r') as f:
    contents = f.readlines()

h5_files= list()
for l in contents:
    if '.h5' in l:
        h5_files.append(l)

file_instances = list()
for h in h5_files:
    #h.split('.')[0].split('_')
    file_instances.append(  h.split('.')[0].split('_')[1:4] )

MT.json_save('Batch02_alex_tracks', path, file_instances)
