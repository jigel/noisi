import os
import numpy as np
import pandas as pd
from math import log, pi
import click
import copy
import json
from glob import glob
from obspy import read, Trace
from obspy.geodetics import gps2dist_azimuth
import matplotlib.pyplot as plt
#ToDo plot if requested.
from noisi.scripts import measurements as rm
#from noisi.scripts import adjnt_functs as af
from noisi.util.windows import my_centered, snratio
from noisi.util.corr_pairs import get_synthetics_filename
# Get and return measurement as a table or something.
from warnings import warn

def get_station_info(stats):

    sta1 = '{}.{}.{}.{}'.format(stats.network,stats.station,stats.location,
    stats.channel)

    try:
        sta2 = '{}.{}.{}.{}'.format(stats.sac.kuser0.strip(),stats.sac.kevnm.strip(),
        stats.sac.kuser1.strip(),stats.sac.kuser2.strip())
    except AttributeError:
        sta2 = '{}.{}.{}.{}'.format(stats.sac.kuser0.strip(),stats.sac.kevnm.strip(),
        '',stats.sac.kuser2.strip())
    lat1 = stats.sac.stla
    lon1 = stats.sac.stlo
    lat2 = stats.sac.evla
    lon2 = stats.sac.evlo
    dist = stats.sac.dist
    az,baz = gps2dist_azimuth(lat1,lon1,lat2,lon2)[1:]


    return([sta1,sta2,lat1,lon1,lat2,lon2,dist,az,baz])


def measurement(source_config,mtype,step,ignore_network,
    bandpass,step_test,**options):
    
    """
    Get measurements on noise correlation data and synthetics. 
    options: g_speed,window_params (only needed if 
    mtype is ln_energy_ratio or enery_diff)
    """
    step_n = 'step_{}'.format(int(step))


    step_dir = os.path.join(source_config['source_path'],
    step_n)

    if step_test:
        corr_dir = os.path.join(step_dir,'obs_slt')
    else:
        corr_dir = os.path.join(source_config['source_path'],
    'observed_correlations')


    files = [f for f in os.listdir(corr_dir) ]

    files = [os.path.join(corr_dir,f) for f in files]

    synth_dir = os.path.join(step_dir,'corr')


    columns = ['sta1','sta2','lat1','lon1','lat2','lon2','dist','az','baz',
    'syn','syn_a','obs','obs_a','l2_norm','snr','snr_a','nstack']
    measurements = pd.DataFrame(columns=columns)

    _options_ac = copy.deepcopy(options)
    _options_ac['window_params']['causal_side'] = not(options['window_params']['causal_side'])

    # ToDo
    if mtype == 'inst_phase':
         _opt_inst = copy.deepcopy(options)


    if files == []:
        msg = 'No input found!'
        raise ValueError(msg)

    i = 0
    with click.progressbar(files,label='Taking measurements...') as bar:

        for f in bar:



            #======================================================
            # Reading
            #======================================================

            try:
                tr_o = read(f)[0]
            except:
                print('\nCould not read data: '+os.path.basename(f))
                i+=1
                continue
            try:
                synth_filename = get_synthetics_filename(os.path.basename(f),
                    synth_dir,ignore_network=ignore_network)

            except:
                print('\nCould not obtain synthetics filename: ' + \
                    os.path.basename(f))
                i+=1
                continue

            if synth_filename is None:
                    continue
                #sfile = glob(os.path.join(synth_dir,synth_filename))[0]
                #print(synth_filename)
            try:
                tr_s = read(synth_filename)[0]
            except:
                print('\nCould not read synthetics: ' + \
                    synth_filename)
                i+=1
                continue

            #======================================================
            # Filtering
            #======================================================
            print(bandpass)
            if bandpass != None:
                tr_o.taper(0.05)
                tr_o.filter('bandpass',freqmin=bandpass[0],
                    freqmax=bandpass[1],corners=bandpass[2],
                    zerophase=True)
                tr_s.taper(0.05)
                tr_s.filter('bandpass',freqmin=bandpass[0],
                    freqmax=bandpass[1],corners=bandpass[2],
                    zerophase=True)

            #======================================================
            # Assigning stats to synthetics, cutting them to right length
            #======================================================

            tr_s.stats.sac = tr_o.stats.sac.copy() #ToDo: Give the stats to this thing before!
            tr_s.data = my_centered(tr_s.data,tr_o.stats.npts)
            # Get all the necessary information
            info = get_station_info(tr_o.stats)

            #======================================================
            # Weight observed stack by nstack
            #======================================================

            tr_o.data /= tr_o.stats.sac.user0



            #======================================================
            # Measurement
            #======================================================

            # Take the measurement
            func = rm.get_measure_func(mtype)

            # ToDo Change this!!!
            if mtype == 'inst_phase':
                _opt_inst['corr_syn'] = tr_s
                try:
                    msr = func(tr_o,**_opt_inst)
                except:
                    print("** Could not take measurement")
                    print(f)
                    continue

            else:
                try:

                    msr_o = func(tr_o,**options)
                    msr_s = func(tr_s,**options)
                except:
                    print("** Could not take measurement")
                    print(f)
                    continue

            # timeseries-like measurements:
            if mtype in ['square_envelope','windowed_envelope','waveform',
            'windowed_waveform']:
                #l2_so = np.trapz(0.5*(msr_s-msr_o)**2) * tr_o.stats.delta
                l2_so = 0.5*np.dot((msr_s-msr_o),np.transpose(msr_s-msr_o))
                snr = snratio(tr_o,**options)
                snr_a = snratio(tr_o,**_options_ac)
                info.extend([np.nan,np.nan,np.nan,np.nan,
                l2_so,snr,snr_a,tr_o.stats.sac.user0])
            # single value measurements:
            else:

                if mtype == 'energy_diff':
                    l2_so = 0.5*(msr_s-msr_o)**2
                    msr = msr_o[0]
                    msr_a = msr_o[1]
                    snr = snratio(tr_o,**options)
                    snr_a = snratio(tr_o,**_options_ac)
                    l2 = l2_so.sum()/2.
                    info.extend([msr_s[0],msr_s[1],msr,msr_a,
                    l2,snr,snr_a,tr_o.stats.sac.user0])
                elif mtype == 'ln_energy_ratio':
                    l2_so = 0.5*(msr_s-msr_o)**2
                    msr = msr_o
                    snr = snratio(tr_o,**options)
                    snr_a = snratio(tr_o,**_options_ac)
                    info.extend([msr_s,np.nan,msr,np.nan,
                    l2_so,snr,snr_a,tr_o.stats.sac.user0])

                elif mtype == 'inst_phase':
                    snr = snratio(tr_o,**options)
                    snr_a = snratio(tr_o,**_options_ac)
                    info.extend([np.nan,np.nan,np.nan,np.nan,
                    msr,snr,snr_a,tr_o.stats.sac.user0])

            measurements.loc[i] = info

            # step index
            i+=1

        return measurements

def run_measurement(source_configfile,measr_configfile,
    step,ignore_network,step_test):


    # get parameters
    source_config=json.load(open(source_configfile))
    measr_config=json.load(open(measr_configfile))
    mtype = measr_config['mtype']
    bandpass = measr_config['bandpass']
    step_n = 'step_{}'.format(int(step))
    step_dir = os.path.join(source_config['source_path'],
    step_n)

    window_params                   =    {}
    window_params['hw']             =    measr_config['window_params_hw']
    
    window_params['sep_noise']      =    measr_config['window_params_sep_noise']
    window_params['win_overlap']    =    measr_config['window_params_win_overlap']
    window_params['wtype']          =    measr_config['window_params_wtype']
    window_params['causal_side']    =    measr_config['window_params_causal']
    window_params['plot']           =    measr_config['window_plot_measurements']
   


    if bandpass == None:
        bandpass = [None]
    if type(bandpass[0]) != list and bandpass[0] != None:
            bandpass = [bandpass]
            warn('\'Bandpass\' should be defined as list of filters.')

    if type(window_params['hw']) != list:
        window_params['hw'] = [window_params['hw']]
    if len(window_params['hw']) != len(bandpass):
        warn('Using the same window length for all measurements.')
        window_params['hw'] = len(bandpass)*[window_params['hw'][0]]
    if type(measr_config['g_speed']) in [float,int]:
        warn('Using the same group velocity for all measurements.')
        g_speeds = len(bandpass)*[measr_config['g_speed']]
    # ToDo: This is ugly and should be sorted out beforehand but 
    # I am too lazy.
    elif type(measr_config['g_speed']) == list \
    and len(measr_config['g_speed']) == len(bandpass):
        g_speeds = measr_config['g_speed']
            

    #if bandpass is None or type(bandpass[0]) != list:
    #    ms = measurement(source_config,mtype,step,ignore_network,bandpass=bandpass,
    #        step_test=step_test,g_speed=g_speed,window_params=window_params)
    #
    #    filename = '{}.0.measurement.csv'.format(mtype)
    #    ms.to_csv(os.path.join(step_dir,filename),index=None)

    #else:

    hws = window_params['hw'][:]

    for i in range(len(bandpass)):

        g_speed = g_speeds[i]

        window_params['hw'] = hws[i]
        ms = measurement(source_config,mtype,step,ignore_network,bandpass=bandpass[i],
        step_test=step_test,g_speed=g_speed,window_params=window_params)

        filename = '{}.{}.measurement.csv'.format(mtype,i)
        ms.to_csv(os.path.join(step_dir,filename),index=None)
