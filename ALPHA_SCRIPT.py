#### Alpha Script
# Script to run the forward model automatically
# Input: alpha_config.json

# PLOTS MIGHT HAVE TO BE REMOVED TO RUN WITHOUT HAVING TO INTERACT

import os
import io
import json
import numpy as np
from pandas import read_csv
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import shutil
from glob import glob
from noisi.util.setup_data_gaussgrid import setup_data_gaussgrid
from noisi.scripts.source_grid_gauss import gauss_grid
from noisi.util.setup_noisesource_new import setup_noisesource_new
from noisi.util.plot_with_azimuth_sel import plot_section
from noisi.util.make_synthetic_data import make_synth_data
from noisi.util.plot_cartopy import plot_gradient

# get the main directory
main_path = os.getcwd()
print(main_path)

# load the alpha_config.json file
with io.open('./alpha_config.json','r') as fh:
    alpha_config = json.load(fh)
    
# get input
n_cores = alpha_config['n_cores']  
project_name = alpha_config['project_name'] 
data_path = alpha_config['data_path'] 
t_data = alpha_config['t_data'] 
f_dom = alpha_config['f_dom']
data_thresh = alpha_config['data_thresh'] 
gamma_thresh = alpha_config['gamma_thresh'] 
extent = alpha_config['extent']
gamma_const = alpha_config['gamma_const']
back_grid_phi = alpha_config['back_grid_phi'] 
data_grid_phi = alpha_config['data_grid_phi']
stationlist_path = alpha_config['stationlist_path'] 
wavefield_path = alpha_config['wavefield_path']

############## Setup new project ##############
os.system('noisi setup_project ' + project_name)
print('New project created:', project_name)

project_path = os.path.join(main_path,project_name)
print('Path to new project: ', project_path)

# copy stationlist file to project and rename stationlist.csv
os.system ('cp {} {}'.format(stationlist_path,os.path.join(project_path,'stationlist.csv')))
print ('Copied stationlist file to project directory.')

# Change current directory to project directory
os.chdir(project_path)
print('Changed to project directory:',os.getcwd())



########## Setup the grid  #################
# Need to make variables for how dense additional grids are

sigma,beta,phi_ini,phi_max,lat_0,lon_0,n,gamma = setup_data_gaussgrid(project_path,data_path,t_data,f_dom,data_thresh,gamma_thresh,gamma_const,back_grid_phi,data_grid_phi,extent,plot=False)

plot = False
only_ocean = True
dense_antipole = False

# Adjust so it can be written to config file
lat_0 = [float(i) for i in lat_0]
lon_0 = [float(i) for i in lon_0]

# Write to config.json file
config_path = os.path.join(project_path,'config.json')


with io.open(config_path,'r+') as fh:
        conf = json.loads(fh.read())
        
conf['gauss_grid'] = True      
conf['gauss_sigma'] = sigma
conf['gauss_beta'] = beta
conf['gauss_phi_ini'] = phi_ini
conf['gauss_phi_max'] = phi_max
conf['gauss_lat_0'] = lat_0
conf['gauss_lon_0'] = lon_0
conf['gauss_n'] = n
conf['gauss_gamma'] = gamma
conf['gauss_plot'] = plot
conf['gauss_dense_antipole'] = dense_antipole
conf['gauss_only_ocean'] = only_ocean

# Set below to true if voronoi cell surface area is to be calculated
conf['voronoi_surface_area'] = True

# change instaseis and wavefield
conf['instaseis'] = True
conf['wavefield_path'] = os.path.abspath(wavefield_path)

with io.open(config_path,'w') as fh:
    cf = json.dumps(conf,sort_keys=False, indent=4, separators=(",", ": "))
    fh.write(cf)


########  calculate grid with noisi ###########
print('Computing grid...')
os.system('noisi setup_sourcegrid .')
print('Grid computed and saved as sourcegrid.npy')
plt.close()

######## Plot and save sourcegrid #########
sourcegrid_path = os.path.join(project_path,'sourcegrid.npy')
grid = np.load(sourcegrid_path)

plt.figure(figsize=(25,10))
ax = plt.axes(projection=ccrs.Mollweide())
ax.coastlines()
plt.scatter(grid[0],grid[1],s=1,c='k',transform=ccrs.Geodetic())
plt.title('Final grid with {} gridpoints'.format(np.size(grid[0])))
plt.savefig(os.path.join(project_path,"sourcegrid.png"))
plt.show(block=False)

######## Setup a homogeneous source called homo_source  ######
source_homo = "homo_source"

source_homo_path = os.path.join(project_path,source_homo)
os.system ('noisi setup_source ' + source_homo_path)

print('New source created: ', source_homo_path)


# Get number of stations for mpi
stationlist = read_csv(stationlist_path)
station_n = np.shape(stationlist)[0]
print('Number of stations: ',station_n)


############# Convert wavefield  ##################
# need arguments: source_config, config, sourcegrid, stationlist, output folder

source_config_path = os.path.join(source_homo_path,'source_config.json')
wavefield_from_instaseis_path = os.path.join(source_homo_path,'wavefield_from_instaseis.py')

print('Converting wavefield from instaseis...')
os.system('mpirun -np {} python {} {} {} {} {} {}'.format(station_n,wavefield_from_instaseis_path,source_config_path,config_path,sourcegrid_path,stationlist_path,project_path))
print('Done.')


wavefield_processed_path = os.path.join(project_path,'wavefield_processed')

# Change config.json file
with io.open(config_path,'r+') as fh:
        conf = json.loads(fh.read())

# change instaseis and wavefield
conf['instaseis'] = False
conf['wavefield_path'] = os.path.abspath(wavefield_processed_path)

with io.open(config_path,'w') as fh:
    cf = json.dumps(conf,sort_keys=False, indent=4, separators=(",", ": "))
    fh.write(cf)


# Change source_config file
source_config_path = os.path.join(source_homo_path,'source_config.json')

with io.open(source_config_path,'r+') as fh:
        conf = json.loads(fh.read())

# change instaseis and wavefield
conf['max_lag'] = 2500
conf['preprocess_do'] = False
conf ['project_name'] = project_name
conf ['project_path'] = project_path
with io.open(source_config_path,'w') as fh:
    cf = json.dumps(conf,sort_keys=False, indent=4, separators=(",", ": "))
    fh.write(cf)


####### setup noise source for homogeneous model   #########

print('Setting up noisesource distribution...')
setup_noisesource_new(project_path,source_homo_path)
print('Done.')


######## calculate correlations for homogeneous source #########
print("Computing correlations...")
os.system("mpirun -np {} noisi correlation {} {}".format(n_cores,source_homo_path,0))
print("All correlations computed.")


##### plot correlations and save the file ####
corr_homo_path = os.path.join(source_homo_path,"step_0/corr")

# make azimuth selection
traces_homo =plot_section(corr_homo_path,stationlist_path,bandpass = None,comp = 'BHZ',fmt = 'SAC',az_selection = [0,180], 
                      scale = 1., resol = 1,plot=False)

#plot
maxlag = (traces_homo[0].stats.npts-1) / 2.0

traces_homo.plot(type='section',orientation='horizontal',
reftime = traces_homo[0].stats.starttime + maxlag,scale=1.,outfile=os.path.join(project_path,'homo_correlations.png'))


########## Setup the Data Source  #########

# Setup a homogeneous source called homo_source
source_data = "data_source"

source_data_path = os.path.join(project_path,source_data)
os.system ('noisi setup_source ' + source_data_path)


print('New source created: ', source_data_path)


# Change source_config file
source_config_path = os.path.join(source_data_path,'source_config.json')

with io.open(source_config_path,'r+') as fh:
        conf = json.loads(fh.read())

# change instaseis and wavefield
conf['max_lag'] = 2500
conf['preprocess_do'] = False
conf ['project_name'] = project_name
conf ['project_path'] = project_path
with io.open(source_config_path,'w') as fh:
    cf = json.dumps(conf,sort_keys=False, indent=4, separators=(",", ": "))
    fh.write(cf)



########## Setup noise source with data   #############
setup_noisesource_new(project_path,source_data_path,data_path = data_path, t_data = t_data, f_dom = f_dom,outfile=True)
# save file with noise source distribution

##########  calculate correlations for homogeneous source #########
print("Computing correlations...")
os.system("mpirun -np {} noisi correlation {} {}".format(n_cores,source_data_path,0))
print("All correlations computed.")

# plot correlations and save the file
corr_data_path = os.path.join(source_data_path,"step_0/corr")

# make azimuth selection
traces_data =plot_section(corr_data_path,stationlist_path,bandpass = None,comp = 'BHZ',fmt = 'SAC',az_selection = [0,180], 
                      scale = 1., resol = 1,plot=False)

#plot
maxlag = (traces_data[0].stats.npts-1) / 2.0

traces_data.plot(type='section',orientation='horizontal', reftime = traces_data[0].stats.starttime + maxlag,scale=1.,outfile=os.path.join(project_path,'data_correlations.png'))



############ copy correlations to observed correlations in homo_source
corr_obs_path = os.path.join(source_homo_path,"observed_correlations")

for files in glob(os.path.join(corr_data_path,'*.sac')):
        shutil.copy(files,corr_obs_path)
        print('Copied:',files)

# make synthetic data
make_synth_data(project_name,stationlist_path,corr_obs_path)
print('Synthetic data created.')



######## noisi measurement ###############
print('Computing adjoint sources...')
os.system('noisi measurement {} {}'.format(source_homo_path,0))
print('Done.')


############ Compute kernels  ##############
print('Computing kernels...')
os.system('mpirun -np {} noisi kernel {} {}'.format(n_cores,source_homo_path,0))
print('Done.')


########## Compute gradient #############
print('Computing gradient...')
os.system('noisi gradient {} {}'.format(source_homo_path,0))
print('Done.')


###### Plot Gradient
plot_gradient(project_path,source_homo_path,extent=extent,gridlines = True)
# save file
plt.savefig(os.path.join(project_path,"data_gradient.png"))

######### Go back to main directory #####
os.chdir(main_path)
print('Directory changed to:',os.getcwd())

print('===================')
print('Alpha Script done. ')
print('===================')
