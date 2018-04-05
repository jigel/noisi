# # Make synthetic data from modelled data

# NOTE: This Jupyter Notebook should be run in /project_name/source_name until arguments can be put in externally.
# Expects stationlist.csv, modelled correlations

# The aim is to calculate the kernels. To do this, observed correlation data is necessary. This script converts modelled data to synthetic data.

import os
import glob
import shutil
import obspy
import pandas as pd
from obspy import read
from obspy.geodetics import gps2dist_azimuth
import numpy as np
import sys
from glob import glob
from obspy.core import AttribDict



# For dataless conversion set it to true
dataless = True
corr_filt_sectionplot = True

# Get Project name
project_name = os.path.basename(os.path.dirname(os.getcwd()))
# print(project_name)

# first get paths to different files
path_stations = ('../stationlist.csv')
path_model = ('./step_0/corr/')
path_obs = ('./observed_correlations/')

# ABOVE CAN BE CHANGED FOR PYTHON SCRIPT TO RUN WITH INPUT ARGUMENTS, see Laura's code

if dataless:
    # delete files in observed_correlations folder if necessary
    for files in glob(os.path.join(path_obs,'*')):
        os.remove(files)
    # copy files from the synthetic correlations to the observed correlations '/Source_1/observed_correlations/'
    for files in glob(os.path.join(path_model,'*.sac')):
        shutil.copy(files,path_obs)
        print 'Copied:',files


# Rename files as noisi expects different filename
for filename in glob(os.path.join(path_obs,'*.sac*')):
    # make sure they're not renamed if they've already been renamed
    if filename.endswith(project_name + '.sac'): 
        break
    else:
        # get filename without extension
        filename_wo_ext = os.path.splitext(filename)[0]
        ext = os.path.splitext(filename)[1]
        # change -- to . and add project name and extension
        filename_1 = filename_wo_ext.replace('--','.')
        filename_2 = filename_1 + '.' + project_name + ext
        # rename the file
        os.rename(filename,filename_2)
        print 'Renamed:', filename_2


# Check metadata in observed_correlations folder
# load the correlations into a file with obspy
ext = '*.sac'
corrs_path_obs = os.path.join(path_obs,ext) # get all .sac files in directory
st = obspy.read(corrs_path_obs) # load all into one stream
# print(st)
#print st[0].stats


# # Laura's code: assign_geodata.py

# Changed the indir and metafile input so it would run in this notebook. 
# For meta = .. engine = 'python' has been added.

#indir = sys.argv[1]
#metafile = sys.argv[2]
indir = path_obs
metafile = '../stationlist.csv'


print(indir)
traces = glob(indir+'/*.SAC')
traces.extend(glob(indir+'/*.sac'))
print('Found traces:\n')
print(traces[0])
print('...to...')
print(traces[-1])
print('\n')
print('Assign geographical information.\n')
print('Number of traces:')
print(np.size(traces))
print('\n')

meta = pd.read_csv(metafile, engine='python')

for t in traces:
    tr = read(t)
    sta1 = os.path.basename(t).split('.')[1]
    try:
        sta2 = os.path.basename(t).split('--')[1].split('.')[1]
    except IndexError:
        sta2 = os.path.basename(t).split('.')[5]
    print(sta1,sta2)
    lat1 = float(meta[meta['sta']==sta1].iloc[0]['lat'])
    lat2 = float(meta[meta['sta']==sta2].iloc[0]['lat'])
    lon1 = float(meta[meta['sta']==sta1].iloc[0]['lon'])
    lon2 = float(meta[meta['sta']==sta2].iloc[0]['lon'])
    print(lat1,lon1,lat2,lon2)
    
    tr[0].stats.network = os.path.basename(t).split('.')[0]
    tr[0].stats.station = sta1
    tr[0].stats.location = ''
    tr[0].stats.channel = os.path.basename(t).split('.')[3] #os.path.basename(t).split('.')[3].split('--')[0]
    tr[0].stats.sac.stlo = lon1
    tr[0].stats.sac.stla = lat1
    tr[0].stats.sac.evlo = lon2
    tr[0].stats.sac.evla = lat2
    tr[0].stats.sac.kuser0 = meta[meta['sta']==sta2].iloc[0]['net']
    
    tr[0].stats.sac.kevnm = sta2
    tr[0].stats.sac.kuser1 = ''
    try:
        tr[0].stats.sac.kuser2 = os.path.basename(t).split('.')[7] #os.path.basename(t).split('--')[1].split('.')[3]
    except IndexError:
        sta2 = os.path.basename(t).split('.')[7]
    tr[0].stats.sac.user0 = 100.   
    #print lat1 > -90.
    #print lat1 < 90.
    #print type(lat1)
    #print(float(lat1))
    #print lat1,lon1,lat2,lon2
    
    geoinf = gps2dist_azimuth(lat1,lon1,lat2,lon2)
    tr[0].stats.sac.dist = geoinf[0]
    tr[0].stats.sac.az = geoinf[1]
    tr[0].stats.sac.baz = geoinf[2]
    tr[0].stats['distance'] = geoinf[0]   # add stats.distance for section plot
    #print tr[0].stats.keys()

    tr.write(t,format='SAC')
    #tr.plot()



# # Back to my code

# Check the metadata again
ext = '*.sac'
corrs_path_obs = os.path.join(path_obs,ext) # get all .sac files in directory
st = obspy.read(corrs_path_obs) # load all into one stream
print(st)
st.plot()


# Plot to see correlations
if corr_filt_sectionplot: 
    st1 = obspy.Stream()
    st2 = obspy.Stream()
    # need to set stats.distance
    for tr in traces:
        t = read(tr)
        t[0].stats.distance = t[0].stats.sac.dist
        #print t[0].stats.distance
        t_filt = t
        t_filt.filter('bandpass',freqmin=0.02,freqmax=0.05,zerophase = True)
        #t_filt.plot()
        #t_filt.plot(type='section')
        st1 += t_filt
        st2 += t

    st1.plot(type='section')
    #st1.spectrogram(log=True,wlen=50)  # spectrogram plot for fun
    #st2.plot(type='section')


# # Option 1: dataless measurement
# All the traces are changed to 1 as if the measurement is 0. 
# That means: synthetic data is compared with dataless observed data.
# Change all trace.data to 1 if dataless = true or leave if dataless = False

print(st)

if dataless:
    for trace in st:
        size = np.size(trace.data)
        trace.data = np.ones(size)
else:
    pass

# plot dataless data
st.plot()


# # Option 2: compare it to a different source distribution
# Set dataless = False for this
# For this, a second source distribution is used to calculate the adjoint source.
# All above steps are run but the trace.data is not changed.
