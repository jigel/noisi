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
try:
    from scipy.fftpack import next_fast_len
except ImportError:
    from noisi.util.scipy_next_fast_len import next_fast_len
from obspy import Trace, read, Stream
from noisi import NoiseSource, WaveField
from noisi.util import geo#, natural_keys
from obspy.signal.invsim import cosine_taper
from noisi.util import filter
try:
    from scipy.signal import sosfilt
except ImportError:
    from obspy.signal._sosfilt import _sosfilt as sosfilt
from noisi.util.windows import my_centered, zero_buddy
from noisi.util.geo import geograph_to_geocent
from noisi.util.corr_pairs import *
import matplotlib.pyplot as plt
import instaseis

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
    # use a one-sided taper: The seismogram probably has a non-zero end, 
    # being cut off whereever the solver stopped running.
    taper = cosine_taper(ntime,p=0.01)
    taper[0:ntime//2] = 1.0

    
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

        adjt_srcs.append(f)
        

    if adjt_srcs_cnt == 0:
        return()
    print(adjt_srcs[0])

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

       
def g1g2_kern_wftype(wf1str,wf2str,kernel,adjt,
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
    # use a one-sided taper: The seismogram probably has a non-zero end, 
    # being cut off whereever the solver stopped running.
    taper = cosine_taper(ntime,p=0.01)
    taper[0:ntime//2] = 1.0

    
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

        adjt_srcs.append(f)
        

    if adjt_srcs_cnt == 0:
        return()
    print(adjt_srcs[0])

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

            