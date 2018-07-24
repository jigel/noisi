
# coding: utf-8


import numpy as np
from obspy.geodetics import gps2dist_azimuth
from obspy.signal.invsim import cosine_taper
import matplotlib.pyplot as plt
import h5py
from noisi import WaveField
import json
from glob import glob
import os
from scipy.fftpack import next_fast_len
from scipy.signal import hann


##################################################################
# USER INPUT
##################################################################
# path to project
projectpath = '../'
sourcepath = '.'


# geography - Add anything else than a homogeneous distribution by setting to "True" the following:
only_ocean = False
gaussian_blobs = False
no_background = False
params_gaussian_blobs = [{'center':(0.,0.),'sigma_radius_m':500000.,'rel_weight':0.00000001}]

#spectra
params_gaussian_spectra = [{'central_freq':0.001,'sigma_freq':0.0002,'weight':10.}]

###############################################################################

grd  = np.load(os.path.join(projectpath,'sourcegrid.npy'))
ntraces = np.shape(grd)[-1]

config = json.load(open(os.path.join(projectpath,'config.json')))
source_config = json.load(open(os.path.join(sourcepath,'source_config.json')))

if source_config['preprocess_do']:
    ext = '*.h5'
    wavefield_path = os.path.join(sourcepath,'wavefield_processed')
else:
    ext = '*.h5'
    wavefield_path = config['wavefield_path']


wfs = glob(os.path.join(wavefield_path,ext))
if wfs != []:
    with WaveField(wfs[0]) as wf:
        df = wf.stats['Fs']
        nt = wf.stats['nt']
        
else:
    df = float(raw_input('Sampling rate in Hz?\n'))
    nt = int(raw_input('Nr of time steps?\n'))

# The number of points for the fft is larger due to zeropadding --> apparent higher frequency sampling\n",
n = next_fast_len(2*nt-1)    
freq = np.fft.rfftfreq(n,d=1./df)
    
taper = cosine_taper(len(freq),0.01)


def get_distance(grid,location):
    def f(lat,lon,location):
        return abs(gps2dist_azimuth(lat,lon,location[0],location[1])[0])
    dist = np.array([f(lat,lon,location) for lat,lon in zip(grid[1],grid[0])])
    return dist
    # Use Basemap to figure out where ocean is
def get_ocean_mask():
    from mpl_toolkits.basemap import Basemap
    m = Basemap(rsphere=6378137,resolution='c',projection='cea',lat_0=0.,
                lon_0=0.,llcrnrlat=-90.,urcrnrlat=90.,llcrnrlon=-180.,urcrnrlon=180.)
    (x,y) = m(grd[0],grd[1])
    #ocean_mask = map(lambda (x,y): not m.is_land(x,y),zip(x,y))
    return ocean_mask


#########################
# Create the source distr
#########################

# geography
num_bases = 1
if gaussian_blobs:
    num_bases += len(params_gaussian_blobs)

basis1 = np.zeros((num_bases,ntraces))

# homogeneous layer
basis1[0,:] = np.ones(ntraces) 

if only_ocean:
    ocean_mask = np.array(get_ocean_mask()).astype(int)
    basis1[0,:] *= ocean_mask

    # superimposed Gaussian blob(s)
if gaussian_blobs:
    i = 1
    for blob in params_gaussian_blobs:
        dist = get_distance(grd,blob['center'])
        basis1[i,:] = np.exp(-(dist)**2/(2*blob['sigma_radius_m']**2))
        
        if only_ocean:
            basis1[i,:] *= ocean_mask
        i+=1


# spectra
basis2 = np.zeros((len(params_gaussian_spectra),len(freq)))
# 'sort of hum gaussian'
i = 0
for spec in params_gaussian_spectra:
    basis2[i,:] = taper*np.exp(-(freq-spec['central_freq'])**2/(2*spec['sigma_freq']**2))
# This normalization means different integrals...
    basis2[i,:] /= np.max(np.abs(basis2[0,:]))
    i+=1


######################
# set the weights
#####################
# geography
weights1 = np.ones(np.shape(basis1)[0])

if gaussian_blobs:
    i = 1
    for blob in params_gaussian_blobs:
        weights1[i] = blob['rel_weight']
        i+=1
    if no_background:
        weights1[0] = 0.
#


from noisi.util import plot



distr = np.dot(weights1,basis1)
plot.plot_grid(grd[0],grd[1],distr,outfile = os.path.join(sourcepath,'geog_distr_startingmodel.png'))


plt.figure()
plt.semilogx(freq,basis2[0,:])
plt.xlabel('Frequency (Hz)')
plt.ylabel('Source power (scaled)')
plt.savefig(os.path.join(sourcepath,'freq_distr_startingmodel.png'))


# Save to an hdf5 file

with h5py.File(os.path.join(sourcepath,'step_0','starting_model.h5'),'w') as fh:
    fh.create_dataset('coordinates',data=grd.astype(np.float32))
    fh.create_dataset('frequencies',data=freq.astype(np.float32))
    fh.create_dataset('distr_basis',data=basis1.astype(np.float32))
    fh.create_dataset('distr_weights',data=weights1.astype(np.float32))
    fh.create_dataset('spect_basis',data=basis2.astype(np.float32))


# Save the 'base model' to an hdf5 file.

basis1_b = np.ones(basis1.shape)
weights1_b = np.ones(weights1.shape)
with h5py.File(os.path.join(sourcepath,'step_0','base_model.h5'),'w') as fh:
    fh.create_dataset('coordinates',data=grd.astype(np.float32))
    fh.create_dataset('frequencies',data=freq.astype(np.float32))
    fh.create_dataset('distr_basis',data=basis1_b.astype(np.float32))
    fh.create_dataset('distr_weights',data=weights1_b.astype(np.float32))
    fh.create_dataset('spect_basis',data=basis2.astype(np.float32))





