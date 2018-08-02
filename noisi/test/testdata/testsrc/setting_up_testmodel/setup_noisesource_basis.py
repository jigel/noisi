##################################################################
# USER INPUT
##################################################################
# path to project and source model
projectpath = '../'
sourcepath = '.'
# Resolution of the coastlines (only relevant for ocean distributions)
# (see basemap documentation)
# Use coarser for global and finer for regional models
coastres = 'i' 
# sampling rate of synthetic Green's function in Hz
sampling_rate = 1.0
# length of synthetic seismograms
n_samples = 3600

################
# geography
################
# list distributions: 'homogeneous', 'ocean','gaussian_blob', 'from_file'
distribution_types = [
'homogeneous'
]
# parameters for homogeneous, ocean: none
# parameters for gaussian blob: center (lat,lon), sigma_radius_m, only_ocean
# parameters for from_file: filename
distribution_params = [
None
]

################
# spectra
################
# list spectra for the above distributions. 'gaussian','from_file'
spectrum_types = ['gaussian']
# parameters for gaussian: mean, standard deviation in Hz
spectrum_params = [ {'mean':0.15,'std':0.02,'weight':5.}]
###############################################################################


import numpy as np
import os
import matplotlib.pyplot as plt
import h5py
from noisi.my_classes.basisfunction import BasisFunction
from noisi.util.source_masks import get_source_mask
from noisi.util.plot import plot_grid
try:
    from scipy.fftpack import next_fast_len
except ImportError:
    from noisi.borrowed_functions.scipy_next_fast_len import next_fast_len
from obspy.signal.invsim import cosine_taper
import json

n = next_fast_len(2*n_samples-1)    
freq = np.fft.rfftfreq(n,d=1./sampling_rate)
print(freq.shape)
taper = cosine_taper(len(freq),0.005)

grd  = np.load(os.path.join(projectpath,'sourcegrid.npy'))
source_config = json.load(open(os.path.join(sourcepath,'source_config.json')))
bfunc_type = source_config['spectra_decomposition']
bfunc_K = source_config['spectra_nr_parameters']


b = BasisFunction(bfunc_type,bfunc_K,N=len(freq))





spectrum_coefficients = np.zeros((len(spectrum_types),bfunc_K))
geographic_weights = np.zeros((len(spectrum_types),grd.shape[-1]))


def gauss_spectrum(sparams):
    spec = taper*np.exp(-(freq-sparams['mean'])**2/
        (2*sparams['std']**2))
    return spec / np.max(np.abs(spec)) * sparams['weight']


for ix_spec in range(len(spectrum_types)):
    
    # get the spectrum
    if spectrum_types[ix_spec] == 'gaussian':
        spectrum = gauss_spectrum(spectrum_params[ix_spec])
    elif spectrum_types[ix_spec] == 'from_file':
        spectrum = np.load(spectrum_params[ix_spec])
    
    # decompose the spectra in the chosen basis
    coeff = b.coeff(spectrum)
    spectrum_coefficients[ix_spec,:] = coeff


# get the geographic weights
for ix_geo in range(len(distribution_types)):

    if distribution_types[ix_geo] =='gaussian_blob':

        geographic_weights[ix_geo,:] = get_source_mask('gaussian',grd,
            coastres,distribution_params[ix_geo])
        print(geographic_weights[ix_geo])

    elif distribution_types[ix_geo] in ['ocean','homogeneous']:

        geographic_weights[ix_geo,:] = get_source_mask(
            distribution_types[ix_geo],grd,coastres)
        
    
    else:
        print(distributions)
        raise NotImplementedError('Unknown geographical distributions. \
            Must be \'gaussian\', \'homogeneous\' or \'ocean\'.')


# get the weighted sum for each location and save

with h5py.File(os.path.join(sourcepath,'step_0','starting_model.h5'),'w') as fh:
    fh.create_dataset('coordinates',data=grd.astype(np.float64))
    fh.create_dataset('frequencies',data=freq.astype(np.float64))
    fh.create_dataset('model',data=np.zeros((grd.shape[-1],bfunc_K)),
        dtype=np.float32)


    for ix_loc in range(grd.shape[-1]):

        for ix_spec in range(len(spectrum_types)):
            #print(geographic_weights[ix_spec,ix_loc])
            fh['model'][ix_loc,:] += geographic_weights[ix_spec,ix_loc] *\
             spectrum_coefficients[ix_spec,:]

        

    fh.flush()
    fh['model'].attrs['spectral_basis'] = bfunc_type
    fh.create_dataset('surf_areas',data=np.ones(grd.shape[-1]))

with h5py.File(os.path.join(sourcepath,'step_0','base_model.h5'),'w') as fh:
    fh.create_dataset('coordinates',data=grd.astype(np.float64))
    fh.create_dataset('frequencies',data=freq.astype(np.float64))
    fh.create_dataset('model',data=np.empty((grd.shape[-1],bfunc_K)),
        dtype=np.float32)


    for ix_loc in range(grd.shape[-1]):

        for ix_spec in range(len(spectrum_types)):
        
            fh['model'][ix_loc,:] += spectrum_coefficients[ix_spec,:]

    fh.flush()
    fh['model'].attrs['spectral_basis'] = bfunc_type
    fh.create_dataset('surf_areas',data=np.ones(grd.shape[-1]))
# plot
for ix_spec in range(len(spectrum_types)):
    spec = np.zeros(freq.shape)
    for i in range(bfunc_K):
        spec += b.basis_vector(i,len(freq)) \
        * spectrum_coefficients[ix_spec,i]
    plt.plot(freq,spec,linewidth=2)

plt.xlabel('Frequency (Hz)')
plt.ylabel('Source power (scaled)')
plt.savefig(os.path.join(sourcepath,'freq_distr_startingmodel.png'))
#
#plt.plot_grid(grd[0],grd[1],colored_by_frequency,
#    normalize=False,sequential=True,cmap='viridis')
