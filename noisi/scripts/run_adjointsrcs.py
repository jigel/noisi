import os
import numpy as np
from math import log, pi
import click
import json
from scipy.signal import hilbert
from glob import glob
from obspy import read, Trace, Stream
from obspy.geodetics import gps2dist_azimuth

from noisi.scripts import adjnt_functs as af
from noisi.util.corr_pairs import get_synthetics_filename
from noisi.util.windows import my_centered, snratio
from warnings import warn

def get_station_info(stats):

    sta1 = '{}.{}.{}.{}'.format(stats.network,stats.station,stats.location,
    stats.channel)
    sta2 = '{}.{}.{}.{}'.format(stats.sac.kuser0.strip(),stats.sac.kevnm.strip(),
    stats.sac.kuser1.strip(),stats.sac.kuser2.strip())
    lat1 = stats.sac.stla
    lon1 = stats.sac.stlo
    lat2 = stats.sac.evla
    lon2 = stats.sac.evlo
    dist = stats.sac.dist
    az = gps2dist_azimuth(lat1,lon1,lat2,lon2)[2]


    return([sta1,sta2,lat1,lon1,lat2,lon2,dist,az])

def get_essential_sacmeta(sac):

    newsacdict={}
    #==============================================================================
    #- Essential metadata
    #==============================================================================

    newsacdict['user0']   =   sac['user0']
    newsacdict['b']       =   sac['b']
    newsacdict['e']       =   sac['e']
    newsacdict['stla']    =   sac['stla']
    newsacdict['stlo']    =   sac['stlo']
    newsacdict['evla']    =   sac['evla']
    newsacdict['evlo']    =   sac['evlo']
    newsacdict['dist']    =   sac['dist']
    newsacdict['az']      =   sac['az']
    newsacdict['baz']     =   sac['baz']
    newsacdict['kuser0']  =   sac['kuser0']
    try:
        newsacdict['kuser1']  =   sac['kuser1']
    except KeyError:
        newsacdict['kuser1'] = ''
    newsacdict['kuser2']  =   sac['kuser2']
    newsacdict['kevnm']   =   sac['kevnm']

    return newsacdict




def adjointsrcs(source_config,mtype,step,ignore_network,bandpass,**options):

    """
    Get 'adjoint source' from noise correlation data and synthetics.
    options: g_speed,window_params (only needed if mtype is ln_energy_ratio or enery_diff)
    """


    files = [f for f in os.listdir(os.path.join(source_config['source_path'],
    'observed_correlations')) ]
    files = [os.path.join(source_config['source_path'],
    'observed_correlations',f) for f in files]


    step_n = 'step_{}'.format(int(step))
    synth_dir = os.path.join(source_config['source_path'],
    step_n,'corr')
    adj_dir = os.path.join(source_config['source_path'],
    step_n,'adjt')



    if files == []:
        msg = 'No input found!'
        raise ValueError(msg)

    #i = 0
    hws = options['window_params']['hw'][:]
    g_speed = options['g_speed'][:]

    with click.progressbar(files,label='Determining adjoint sources...') as bar:

        for f in bar:

            # read data
            try:
                tr_o = read(f)[0]
            except:
                print('\nCould not read data: '+os.path.basename(f))
                #i+=1
                continue

            # read synthetics
            try:
                synth_filename = get_synthetics_filename(os.path.basename(f),
                    synth_dir,ignore_network=ignore_network)
                if synth_filename is None:
                    continue
                #sname = glob(os.path.join(synth_dir,synth_filename))[0]
                print(synth_filename)
                tr_s = read(synth_filename)[0]

            except:
                print('\nCould not read synthetics: '+os.path.basename(f))
                #i+=1
                continue

            # Add essential metadata
            tr_s.stats.sac = get_essential_sacmeta(tr_o.stats.sac)

            # Check sampling rates.
            if round(tr_s.stats.sampling_rate,6) != round(tr_o.stats.
                sampling_rate,6):
                print("Sampling Rates (Hz):\n")
                print(tr_s.stats.sampling_rate)
                print(tr_o.stats.sampling_rate)
                msg = 'Sampling rates of data and synthetics must match.'
                raise ValueError(msg)

            # Waveforms must have same nr of samples.
            tr_s.data = my_centered(tr_s.data,tr_o.stats.npts)

            func = af.get_adj_func(mtype)

            # ugly...sorry

            # Bandpasses
            for j in range(len(bandpass)):

                options['window_params']['hw'] = hws[j]
                options['g_speed'] = g_speed[j]

                tr_o_filt = tr_o.copy()
                tr_s_filt = tr_s.copy()
                # trace with the taper and filter for envelope misfit:
                tr_t_filt = tr_s.copy()
                tr_t_filt.data = np.ones(tr_t_filt.stats.npts)


                bp = bandpass[j]
                tr_o_filt.taper(0.1)
                tr_s_filt.taper(0.1)
                tr_t_filt.taper(0.1) # filtering artifacts will occur 
                # right now.

                if bp != None:
                    tr_o_filt.filter('bandpass',freqmin=bp[0],freqmax=bp[1],
                        corners=bp[2],zerophase=True)
                    tr_s_filt.filter('bandpass',freqmin=bp[0],freqmax=bp[1],
                        corners=bp[2],zerophase=True)
                    tr_t_filt.filter('bandpass',freqmin=bp[0],freqmax=bp[1],
                        corners=bp[2],zerophase=True)
                    tr_t_filt.taper(0.1)


                if mtype == 'square_envelope':
                    options['taper_filter'] = tr_t_filt
                # Get the adjoint source
                data, success = func(tr_o_filt,tr_s_filt,**options)
                if not success:
                    continue

                adj_src = Stream()

                if isinstance(data,list):

                    adj_src += Trace(data=data[0])
                    adj_src += Trace(data=data[1])
                    brnchs = ['c','a']
                    for k in range(2):
                        adjtrc = adj_src[k]
                        adjtrc.stats.sampling_rate = tr_s.stats.sampling_rate
                        adjtrc.stats.sac = tr_s.stats.sac.copy()
                        # Save the adjoint source
                        file_adj_src = os.path.join(adj_dir,
                            os.path.basename(synth_filename).
                            rstrip('sac')+'{}.{}.sac'.format(brnchs[k],j))
                        adjtrc.write(file_adj_src,format='SAC')


                else:
                    adj_src += Trace(data=data)
                    for adjtrc in adj_src:
                        adjtrc.stats.sampling_rate = tr_s.stats.sampling_rate
                        adjtrc.stats.sac = tr_s.stats.sac.copy()
                        # Save the adjoint source
                        file_adj_src = os.path.join(adj_dir,
                        os.path.basename(synth_filename).
                            rstrip('sac')+'{}.sac'.format(j))
                        adjtrc.write(file_adj_src,format='SAC')
    return()



def run_adjointsrcs(source_configfile,measr_configfile,step,ignore_network):

    source_config=json.load(open(source_configfile))
    measr_config=json.load(open(measr_configfile))

    g_speed = measr_config['g_speed']
    mtype = measr_config['mtype']
    bandpass = measr_config['bandpass']

    if bandpass == None:
        bandpass = [None]

    if type(bandpass[0]) != list and bandpass[0] != None:
            bandpass = [bandpass]
            warn('\'Bandpass\' should be defined as list of filters.')

    window_params = {}
    window_params['hw'] = measr_config['window_params_hw']
    if type(window_params['hw']) != list:
        window_params['hw'] = [window_params['hw']]
    if len(window_params['hw']) != len(bandpass):
        warn('Using the same window length for all measurements.')
        window_params['hw'] = len(bandpass)*[window_params['hw'][0]]
    if type(measr_config['g_speed']) in [float,int]:
        warn('Using the same group velocity for all measurements.')
        g_speed = len(bandpass)*[measr_config['g_speed']]
    elif type(measr_config['g_speed']) == list \
    and len(measr_config['g_speed']) == len(bandpass):
        g_speed = measr_config['g_speed']

    window_params['sep_noise'] = measr_config['window_params_sep_noise']
    window_params['win_overlap'] = measr_config['window_params_win_overlap']
    window_params['wtype'] = measr_config['window_params_wtype']
    window_params['causal_side'] = measr_config['window_params_causal']
    window_params['plot'] = False 
    adjointsrcs(source_config,mtype,step,ignore_network=ignore_network,g_speed=g_speed,
        bandpass=bandpass,window_params=window_params)
