    # # Make synthetic data from modelled data
    
    
    
    # The aim is to calculate the kernels. To do this, observed correlation data is necessary. This script converts modelled data to synthetic data.
    
    
def make_synth_data(project_name, stations_path, corr_path):
    
    """
    Input: project_name, path to stationlist.csv file, path to correlations,
    """
    
    
    import os
    #import glob
    #import shutil
    import obspy
    import pandas as pd
    from obspy import read
    from obspy.geodetics import gps2dist_azimuth
    import numpy as np
    #import sys
    from glob import glob
    #from obspy.core import AttribDict
    
    
    
    # Get Project name
    #project_name = project_name
    # print(project_name)
    
    # first get paths to different files
    #path_stations = stations_path
    #ath_model = corr_path
    path_obs = corr_path
        
    # Rename files as noisi expects different filename
    for filename in glob(os.path.join(path_obs,'*.sac*')):
        # make sure they're not renamed if they've already been renamed
        if filename.endswith(project_name + '.sac'): 
            print('Already renamed: ',filename)
            continue
        else:
            # get filename without extension
            filename_wo_ext = os.path.splitext(filename)[0]
            ext = os.path.splitext(filename)[1]
            # change -- to . and add project name and extension
            filename_1 = filename_wo_ext.replace('--','.')
            filename_2 = filename_1 + '.' + project_name + ext
            # rename the file
            os.rename(filename,filename_2)
            print('Renamed:', filename_2)
    
    
    # Check metadata in observed_correlations folder
    # load the correlations into a file with obspy
    ext = '*.sac'
    corrs_path_obs = os.path.join(path_obs,ext) # get all .sac files in directory
    st = obspy.read(corrs_path_obs) # load all into one stream
    
    
    # # Laura's code: assign_geodata.py
    
    # Changed the indir and metafile input so it would run in this notebook. 
    # For meta = .. engine = 'python' has been added.
    
    #indir = sys.argv[1]
    #metafile = sys.argv[2]
    indir = path_obs
    metafile = stations_path
    
    
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
        #print(lat1 > -90.)
        #print(lat1 < 90.)
        #print(type(lat1))
        #print(float(lat1))
        #print(lat1,lon1,lat2,lon2)
        
        geoinf = gps2dist_azimuth(lat1,lon1,lat2,lon2)
        tr[0].stats.sac.dist = geoinf[0]
        tr[0].stats.sac.az = geoinf[1]
        tr[0].stats.sac.baz = geoinf[2]
        tr[0].stats['distance'] = geoinf[0]   # add stats.distance for section plot
        #print(tr[0].stats.keys())
    
        tr.write(t,format='SAC')
        #tr.plot()
    
    
    
    # # Back to my code
    
    # Check the metadata again
    ext = '*.sac'
    corrs_path_obs = os.path.join(path_obs,ext) # get all .sac files in directory
    st = obspy.read(corrs_path_obs) # load all into one stream
    print(st)
    #st.plot()
    

