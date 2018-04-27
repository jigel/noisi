from __future__ import print_function
from mpi4py import MPI
import numpy as np
import os
import h5py
import json
import click
from glob import glob
from math import ceil
from scipy.signal.signaltools import fftconvolve
from scipy.fftpack import next_fast_len
from obspy import Trace, read, Stream
from noisi import NoiseSource, WaveField
from noisi.util import geo, natural_keys
from obspy.signal.invsim import cosine_taper
from noisi.util import filter
from scipy.signal import sosfilt
from noisi.util.windows import my_centered, zero_buddy
from noisi.util.geo import geograph_to_geocent
from noisi.util.corr_pairs import *
import matplotlib.pyplot as plt
import instaseis


#ToDo: put in the possibility to run on mixed channel pairs
def paths_input(cp,source_conf,step,ignore_network,instaseis):
    
    inf1 = cp[0].split()
    inf2 = cp[1].split()
    
    conf = json.load(open(os.path.join(source_conf['project_path'],
        'config.json')))
    measr_conf = json.load(open(os.path.join(source_conf['source_path'],
        'measr_config.json')))
    channel = source_conf['channel']
    
    # station names
    if ignore_network:
        sta1 = "*.{}..{}".format(*(inf1[1:2]+[channel]))
        sta2 = "*.{}..{}".format(*(inf2[1:2]+[channel]))
    else:
        sta1 = "{}.{}..{}".format(*(inf1[0:2]+[channel]))
        sta2 = "{}.{}..{}".format(*(inf2[0:2]+[channel]))


    # Wavefield files  
    if instaseis == False:
        if source_conf['preprocess_do']:
            dir = os.path.join(source_conf['source_path'],'wavefield_processed')
            
        else:
            dir = conf['wavefield_path']
    
        wf1 = glob(os.path.join(dir,sta1+'.h5'))[0]
        wf2 = glob(os.path.join(dir,sta2+'.h5'))[0]
    else:
        # need to return two receiver coordinate pairs. 
        # For buried sensors, depth could be used but no elevation is possible,
        # so maybe keep everything at 0 m?
        # lists of information directly from the stations.txt file.
        wf1 = inf1
        wf2 = inf2

    
    # Starting model for the noise source
   
        # The base model contains no spatial or spectral weights.
    nsrc = os.path.join(source_conf['project_path'],
                 source_conf['source_name'],'step_'+str(step),
                 'base_model.h5')
  
    # Adjoint source
    if measr_conf['mtype'] == 'energy_diff':
        adj_src_basicnames = [ os.path.join(source_conf['source_path'],
                 'step_'+str(step),
                 'adjt',"{}--{}.c".format(sta1,sta2)),
                 os.path.join(source_conf['source_path'],
                 'step_'+str(step),
                 'adjt',"{}--{}.a".format(sta1,sta2))]
    else:
        adj_src_basicnames = [os.path.join(source_conf['source_path'],
                 'step_'+str(step),
                 'adjt',"{}--{}".format(sta1,sta2))]


     
    return(wf1,wf2,nsrc,adj_src_basicnames)
    
    
def paths_output(cp,source_conf,step):
    

    id1 = cp[0].split()[0]+cp[0].split()[1]
    id2 = cp[1].split()[0]+cp[1].split()[1]

    if id1 < id2 :
        inf1 = cp[0].split()
        inf2 = cp[1].split()
    else:
        inf2 = cp[0].split()
        inf1 = cp[1].split()

    channel = source_conf['channel']
    sta1 = "{}.{}..{}".format(*(inf1[0:2]+[channel]))
    sta2 = "{}.{}..{}".format(*(inf2[0:2]+[channel]))
    

    kern_basicname = "{}--{}".format(sta1,sta2)
    kern_basicname = os.path.join(source_conf['source_path'],
        'step_'+str(step), 'kern',
        kern_basicname)

    return (kern_basicname)
    
def get_ns(wf1,source_conf,insta):
    
    # Nr of time steps in traces
    if insta:
        # get path to instaseis db
        #ToDo: ugly.
        dbpath = json.load(open(os.path.join(source_conf['project_path'],
            'config.json')))['wavefield_path']
        # open 
        db = instaseis.open_db(dbpath)
        # get a test seismogram to determine...
        stest = db.get_seismograms(source=instaseis.ForceSource(latitude=0.0,
            longitude=0.0),receiver=instaseis.Receiver(latitude=10.,
            longitude=0.0),dt=1./source_conf['sampling_rate'])[0]
        
        nt = stest.stats.npts
        Fs = stest.stats.sampling_rate
    else:
        with WaveField(wf1) as wf1:
            nt = int(wf1.stats['nt'])
            Fs = round(wf1.stats['Fs'],8)
    
    # Necessary length of zero padding for carrying out 
    # frequency domain correlations/convolutions
    n = next_fast_len(2*nt-1)     
    
    # Number of time steps for synthetic correlation
    n_lag = int(source_conf['max_lag'] * Fs)
    if nt - 2*n_lag <= 0:
        click.secho('Resetting maximum lag to %g seconds: Synthetics are too\
 short for a maximum lag of %g seconds.' %(nt//2/Fs,n_lag/Fs))
        n_lag = nt // 2
        
    n_corr = 2*n_lag + 1
    
    return nt,n,n_corr,Fs
        
    


def g1g2_kern(wf1str,wf2str,kernel,adjt,
    src,source_conf,insta):
        
    measr_conf = json.load(open(os.path.join(source_conf['source_path'],
        'measr_config.json')))


    bandpass = measr_conf['bandpass']

    if bandpass == None:
        filtcnt = 1
    elif type(bandpass) == list:
        if type(bandpass[0]) != list:
            filtcnt = 1
        else:
            filtcnt = len(bandpass)    
    
    ntime, n, n_corr, Fs = get_ns(wf1str,source_conf,insta)
    taper = cosine_taper(ntime,p=0.05)

    
########################################################################
# Prepare filenames and adjoint sources
########################################################################   

    filenames = []
    adjt_srcs = []
    adjt_srcs_cnt = 0

    for ix_f in range(filtcnt):
    
        filename = kernel+'.{}.npy'.format(ix_f)
        print(filename)
        filenames.append(filename)
        #if os.path.exists(filename):
         #   continue

        f = Stream()
        for a in adjt:
            adjtfile = a + '*.{}.sac'.format(ix_f)
            adjtfile = glob(adjtfile)
            try:    
                f += read(adjtfile[0])[0]
                f[-1].data = my_centered(f[-1].data,n_corr)
                adjt_srcs_cnt += 1
            except IndexError:
                print('No adjoint source found: {}\n'.format(a))
                break

        if len(f) == len(adjt):
            adjt_srcs.append(f)
        else:
            adjt_srcs.append(None)

    if adjt_srcs_cnt == 0:
        return()

########################################################################
# Compute the kernels
######################################################################## 


    with NoiseSource(src) as nsrc:

        
        ntraces = nsrc.src_loc[0].shape[0]


        if insta:
            # open database
            dbpath = json.load(open(os.path.join(source_conf['project_path'],
                'config.json')))['wavefield_path']
            # open and determine Fs, nt
            db = instaseis.open_db(dbpath)
            # get receiver locations
            lat1 = geograph_to_geocent(float(wf1[2]))
            lon1 = float(wf1[3])
            rec1 = instaseis.Receiver(latitude=lat1,longitude=lon1)
            lat2 = geograph_to_geocent(float(wf2[2]))
            lon2 = float(wf2[3])
            rec2 = instaseis.Receiver(latitude=lat2,longitude=lon2)

        else:
            wf1 = WaveField(wf1str)
            wf2 = WaveField(wf2str)

        kern = np.zeros((filtcnt,ntraces,len(adjt)))

    

        
            
        
        ########################################################################
        # Loop over locations
        ########################################################################            
        for i in range(ntraces):

            # noise source spectrum at this location
            # For the kernel, this contains only the basis functions of the 
            # spectrum without weights; might still be location-dependent, 
            # for example when constraining sensivity to ocean
            S = nsrc.get_spect(i)
            

            if S.sum() == 0.: 
            # The spectrum has 0 phase so only checking absolute value here
                continue

            ####################################################################
            # Get synthetics
            ####################################################################                
            if insta:
            # get source locations
                lat_src = geograph_to_geocent(nsrc.src_loc[1,i])
                lon_src = nsrc.src_loc[0,i]
                fsrc = instaseis.ForceSource(latitude=lat_src,
                    longitude=lon_src,f_r=1.e12)
                
                s1 = np.ascontiguousarray(db.get_seismograms(source=fsrc,
                    receiver=rec1,
                    dt=1./source_conf['sampling_rate'])[0].data*taper)
                s2 = np.ascontiguousarray(db.get_seismograms(source=fsrc,
                    receiver=rec2,
                    dt=1./source_conf['sampling_rate'])[0].data*taper)
                

            else:
                s1 = np.ascontiguousarray(wf1.data[i,:]*taper)
                s2 = np.ascontiguousarray(wf2.data[i,:]*taper)
            
            

            spec1 = np.fft.rfft(s1,n)
            spec2 = np.fft.rfft(s2,n)
            
          
            g1g2_tr = np.multiply(np.conjugate(spec1),spec2)
            c = np.multiply(g1g2_tr,S)

        #######################################################################
        # Get Kernel at that location
        #######################################################################   
            corr_temp = my_centered(np.fft.ifftshift(np.fft.irfft(c,n)),n_corr)
            
        #######################################################################
        # Apply the 'adjoint source'
        #######################################################################
            for ix_f in range(filtcnt):
                f = adjt_srcs[ix_f]

                if f==None:
                    continue
                for j in range(len(f)):
                    delta = f[j].stats.delta
                    kern[ix_f,i,j] = np.dot(corr_temp,f[j].data) * delta
                    print(kern[ix_f,i,j])
           
            
            if i%50000 == 0:
                print("Finished {} source locations.".format(i))


    if not insta:
        wf1.file.close()
        wf2.file.close()

    for ix_f in range(filtcnt):
        filename = filenames[ix_f]
        if kern[ix_f,:,:].sum() != 0:
            np.save(filename,kern[ix_f,:,:]) 
    return()

       
            


def run_kern(source_configfile,step,ignore_network=False):


    # simple embarrassingly parallel run:

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    step = int(step)


    #ToDo think about that configuration decorator
    source_config=json.load(open(source_configfile))
    obs_only = source_config['model_observed_only']
    #ToDo: ugly.
    insta = json.load(open(os.path.join(source_config['project_path'],
        'config.json')))['instaseis']
    
    p = define_correlationpairs(source_config['project_path'])
    if rank == 0:
        print('Nr all possible kernels %g ' %len(p))
    
    # Remove pairs for which no observation is available
    if obs_only:
        directory = os.path.join(source_config['source_path'],
            'observed_correlations')
        p = rem_no_obs(p,source_config,directory=directory)
        if rank == 0:
            print('Nr kernels after checking available observ. %g ' %len(p))
    


    # The assignment of station pairs should be such that one core 
    # has as many occurrences of the same station as possible; 
    # this will prevent that many processes try to access the 
    # same hdf5 file all at once.
    num_pairs = int( ceil(float(len(p))/float(size)) )
    p_p = p[ rank*num_pairs : rank*num_pairs + num_pairs] 
    
    print('Rank number %g' %rank)
    print('working on pair nr. %g to %g of %g.' %(rank*num_pairs,
        rank*num_pairs+num_pairs,len(p)))


    
    for cp in p_p:
        
        try:
            wf1,wf2,src,adjt = paths_input(cp,source_config,
                step,ignore_network,insta)
            print(wf1,wf2,src)
            kernel = paths_output(cp,source_config,step)
            print(kernel)
     
        except:
            print('Could not find input for: %s\
\nCheck if wavefield .h5 file and base_model file are available.' %cp)
            continue

        kern = g1g2_kern(wf1,wf2,kernel,adjt,src,source_config,insta=insta)
    return()

