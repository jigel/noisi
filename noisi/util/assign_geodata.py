import pandas as pd
from obspy import read
from obspy.geodetics import gps2dist_azimuth

import sys
import os
from glob import glob

indir = sys.argv[1]
metafile = sys.argv[2]

print(indir)
traces = glob(indir+'/*.SAC')
traces.extend(glob(indir+'/*.sac'))
print('Found traces:\n')
print(traces[0])
print('...to...')
print(traces[-1])
print('Assign geographical information.\n')

meta = pd.read_csv(metafile)

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
    tr[0].stats.channel = os.path.basename(t).split('.')[3].split('--')[0]
    tr[0].stats.sac.stlo = lon1
    tr[0].stats.sac.stla = lat1
    tr[0].stats.sac.evlo = lon2
    tr[0].stats.sac.evla = lat2
    tr[0].stats.sac.kuser0 = meta[meta['sta']==sta2].iloc[0]['net']
    
    tr[0].stats.sac.kevnm = sta2
    tr[0].stats.sac.kuser1 = ''
    try:
        tr[0].stats.sac.kuser2 = os.path.basename(t).split('--')[1].split('.')[3]
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

    tr.write(t,format='SAC')
