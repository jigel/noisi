#### Alpha Script
# Script to run the forward model automatically
# Input: alpha_config.json

# PLOTS MIGHT HAVE TO BE REMOVED TO RUN WITHOUT HAVING TO INTERACT

import os
import io
import json
import sys
import numpy as np
from pandas import read_csv
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import shutil
from glob import glob
from noisi.scripts.source_grid import setup_sourcegrid as setup_sgrid
from noisi.scripts.source_grid_gauss import setup_sourcegrid_gauss as setup_sgrid_gauss
from noisi.util.setup_data_gaussgrid import setup_data_gaussgrid
from noisi.util.setup_noisesource_new import setup_noisesource_new
from noisi.util.plot_with_azimuth_sel import plot_section
from noisi.util.make_synthetic_data import make_synth_data
from noisi.util.plot_cartopy import plot_gradient
#from noisi.main import setup_proj
#from noisi.main import setup_sourcegrid
#from noisi.main import setup_source
#from noisi.main import correlation
#from noisi.main import kernel
#from noisi.main import gradient
from subprocess import call
from noisi.util.setup_new import setup_proj
import time

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

if os.path.exists(project_name):
    print('Project exists already, must give it a new name.')

else:
    setup_proj(project_name)
    print('New project created:', project_name)

project_path = os.path.join(main_path,project_name)
print('Path to new project: ', project_path)

# copy stationlist file to project and rename stationlist.csv
call('cp {} {}'.format(stationlist_path,os.path.join(project_path,'stationlist.csv')),shell=True)
print ('Copied stationlist file to project directory.')


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
# calculate grid with noisi
print('Computing grid...')
conf = json.load(open(os.path.join(project_path,'config.json')))

if conf['gauss_grid']:
    setup_sgrid_gauss(os.path.join(project_path,'config.json'))
else:
    setup_sgrid(os.path.join(project_path,'config.json'))
print('Grid computed and saved as sourcegrid.npy')

######## Plot and save sourcegrid #########
sourcegrid_path = os.path.join(project_path,'sourcegrid.npy')
grid = np.load(sourcegrid_path)

print('Number of gridpoints: ', np.size(grid[0]))

plt.figure(figsize=(25,10))
ax = plt.axes(projection=ccrs.Mollweide())
ax.coastlines()
plt.scatter(grid[0],grid[1],s=1,c='k',transform=ccrs.Geodetic())
plt.title('Final grid with {} gridpoints'.format(np.size(grid[0])))
plt.savefig(os.path.join(project_path,"sourcegrid.png"))
plt.show(block=False)
plt.close()

######## Setup a homogeneous source called homo_source  ######
# Setup a homogeneous source called homo_source
source_homo = "homo_source"
source_homo_path = os.path.join(project_path,source_homo)
source_model = source_homo_path
    
if not os.path.exists('config.json'):
    print('No config file for project found \
    (detailing e.g. source grid). Run setup_project first.')

if os.path.exists(source_model):
    print('Source exists already, must give it a new name.')
    
else:

    os.makedirs(os.path.join(source_model,'step_0'))
    os.mkdir(os.path.join(source_model,'observed_correlations'))

    for d in ['adjt','grad','corr','kern']:
        os.mkdir(os.path.join(source_model,'step_0',d))

    from noisi import _ROOT
    
    with io.open(os.path.join(_ROOT,'config','source_config.json'),'r') as fh:
        conf = json.loads(fh.read())
        conf['date_created'] = str(time.strftime("%Y.%m.%d"))
        conf['project_name'] = os.path.basename(os.getcwd())
        conf['project_path'] = os.getcwd()
        conf['source_name'] = source_model
        conf['source_path'] = os.path.abspath(source_model)

    with io.open(os.path.join(source_model,'source_config.json'),'w') as fh:
        cf = json.dumps(conf,sort_keys=True, indent=4, separators=(",", ": "))
        fh.write(cf)

    with io.open(os.path.join(_ROOT,'config','measr_config.json'),'r') as fh:
        conf = json.loads(fh.read())
        conf['date_created'] = str(time.strftime("%Y.%m.%d"))

    with io.open(os.path.join(source_model,'measr_config.json'),'w') as fh:
        cf = json.dumps(conf,sort_keys=True, indent=4, separators=(",", ": "))
        fh.write(cf)

    #from noisi import _ROOT
    #os.system('cp {} {}'.format(os.path.join(_ROOT,'jnotebks/\
    #setup_noisesource.ipynb'),
    #source_model))  
    #os.system('cp {} {}'.format(os.path.join(_ROOT,'util/setup_noisesource.py'),
    #source_model))
    #os.system('cp {} {}'.format(os.path.join(_ROOT,
    #    'util/wavefield_from_instaseis.py'),source_model))
    #print("Copied default source_config.json and measr_config.json \
    #to source model directory, please edit. \
    #Please run setup_noisesource.ipynb or setup_noisesource.py after editing to \
    #create starting model.")

print('New source created: ', source_homo_path)


# Get number of stations for mpi
stationlist = read_csv(stationlist_path)
station_n = np.shape(stationlist)[0]
print('Number of stations: ',station_n)


############# Convert wavefield  ##################
# need arguments: source_config, config, sourcegrid, stationlist, output folder

source_config_path = os.path.join(source_homo_path,'source_config.json')
# copy wavefield_from_instaseis.py
call('cp {} {}'.format(os.path.join(main_path,'noisi/util/wavefield_from_instaseis.py'),source_homo_path),shell=True)
print('Copied wavefield_from_instaseis.py file to source directory.')


wavefield_from_instaseis_path = os.path.join(source_homo_path,'wavefield_from_instaseis.py')

print('Converting wavefield from instaseis...')
call('mpirun -np {} python {} {} {} {} {} {}'.format(n_cores,wavefield_from_instaseis_path,source_config_path,config_path,sourcegrid_path,stationlist_path,project_path),shell=True)
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
call("mpirun -np {} noisi correlation {} {}".format(n_cores,source_homo_path,0),shell=True)
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

## EXIT if config 'compute' = 'correlation'
if alpha_config['compute'] == 'correlation':
    os.chdir(main_path)
    print('Directory changed to:',os.getcwd())
    sys.exit('##### Correlations computed. Exiting. #####')

########## Setup the Data Source  #########

source_data = "data_source"
source_data_path = os.path.join(project_path,source_homo)
source_model = source_homo_path


    
if not os.path.exists('config.json'):
    print('No config file for project found \
    (detailing e.g. source grid). Run setup_project first.')

if os.path.exists(source_model):
    print('Source exists already, must give it a new name.')
    
else:

    os.makedirs(os.path.join(source_model,'step_0'))
    os.mkdir(os.path.join(source_model,'observed_correlations'))

    for d in ['adjt','grad','corr','kern']:
        os.mkdir(os.path.join(source_model,'step_0',d))

    from noisi import _ROOT
    
    with io.open(os.path.join(_ROOT,'config','source_config.json'),'r') as fh:
        conf = json.loads(fh.read())
        conf['date_created'] = str(time.strftime("%Y.%m.%d"))
        conf['project_name'] = os.path.basename(os.getcwd())
        conf['project_path'] = os.getcwd()
        conf['source_name'] = source_model
        conf['source_path'] = os.path.abspath(source_model)

    with io.open(os.path.join(source_model,'source_config.json'),'w') as fh:
        cf = json.dumps(conf,sort_keys=True, indent=4, separators=(",", ": "))
        fh.write(cf)

    with io.open(os.path.join(_ROOT,'config','measr_config.json'),'r') as fh:
        conf = json.loads(fh.read())
        conf['date_created'] = str(time.strftime("%Y.%m.%d"))

    with io.open(os.path.join(source_model,'measr_config.json'),'w') as fh:
        cf = json.dumps(conf,sort_keys=True, indent=4, separators=(",", ": "))
        fh.write(cf)

    #from noisi import _ROOT
    #os.system('cp {} {}'.format(os.path.join(_ROOT,'jnotebks/\
    #setup_noisesource.ipynb'),
    #source_model))  
    #os.system('cp {} {}'.format(os.path.join(_ROOT,'util/setup_noisesource.py'),
    #source_model))
    #os.system('cp {} {}'.format(os.path.join(_ROOT,
    #    'util/wavefield_from_instaseis.py'),source_model))
    #print("Copied default source_config.json and measr_config.json \
    #to source model directory, please edit. \
    #Please run setup_noisesource.ipynb or setup_noisesource.py after editing to \
    #create starting model.")

print('New source created: ', source_homo_path)


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
call("mpirun -np {} noisi correlation {} {}".format(n_cores,source_data_path,0),shell=True)
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
call('noisi measurement {} {}'.format(source_homo_path,0),shell=True)
print('Done.')


############ Compute kernels  ##############
print('Computing kernels...')
call('mpirun -np {} noisi kernel {} {}'.format(n_cores,source_homo_path,0),shell=True)
print('Done.')

## EXIT if config 'compute' = 'kernel'
if alpha_config['compute'] == 'kernel':
    os.chdir(main_path)
    print('Directory changed to:',os.getcwd())
    sys.exit('##### Kernel computed. Exiting. #####')

########## Compute gradient #############
print('Computing gradient...')
call('noisi gradient {} {}'.format(source_homo_path,0),shell=True)
print('Done.')


###### Plot Gradient
plot_gradient(project_path,source_homo_path,extent=extent,gridlines = True)
# save file
plt.savefig(os.path.join(project_path,"data_gradient.png"))

## EXIT if config 'compute' = 'gradient'
if alpha_config['compute'] == 'gradient':
    os.chdir(main_path)
    print('Directory changed to:',os.getcwd())
    
    sys.exit('##### Gradient computed. Exiting. #####')

