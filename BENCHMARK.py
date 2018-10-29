# code to do a benchmark test for the grid

import os
import io
import json
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from pandas import read_csv
from noisi.util.setup_noisesource_new import setup_noisesource_new
from noisi.util.plot_with_azimuth_sel import plot_section



# get main path
main_path = os.getcwd()
print(main_path)


# homogeneous grids input

gridpoints_homo = [1018,2543,5092,7499,10071,20125,30048,40748,50993,60079,70771,81857,101004]

phi_ini_homo = [5.4,3.4,2.4,1.98,1.71,1.21,0.99,0.85,0.76,0.7,0.645,0.6,0.54]

sigma_homo = [10]*13
beta_homo = [5]*13
phi_max_homo = phi_ini_homo
lat_0_homo = [0]*13
lon_0_homo = [0]*13
n_homo = [1000]*13
gamma_homo = [0]*13


# multiple gaussian grids input

gridpoints_gauss = [1080,2510,5106,7331,10330,20840,29984,40095,50000,61968,73311,81861,103365]

sigma_gauss = [[10,10],[10,10],[15,15],[16,16],[26,26],[20,20],[23,23],[27,27],[28,28],[31,31],[35,35],[34,34],[3,3]]
beta_gauss = [[5,5],[5,5],[5,5],[5,5],[5,5],[5,5],[5,5],[5,5],[5,5],[5,5],[5,5],[5,5],[5,5]]
phi_ini_gauss = [[3.5,3.5],[2,2],[1.2,1.2],[1,1],[0.8,0.8],[0.5,0.5],[0.4,0.4],[0.3,0.3],[0.3,0.3],[0.2,0.2],[0.1,0.1],[0.1,0.1],[0.1,0.1]]
phi_max_gauss = [[10,10],[5,5],[4,4],[3,3],[4,4],[4,4],[3,3],[3,3],[2,2],[2,2],[2,2],[1.5,1.5],[1,1]]
lat_0_gauss = [[50,-16],[50,-16],[50,-16],[50,-16],[50,-16],[50,-16],[50,-16],[50,-16],[50,-16],[50,-16],[50,-16],[50,-16],[50,-16]]
lon_0_gauss = [[-30,70],[-30,70],[-30,70],[-30,70],[-30,70],[-30,70],[-30,70],[-30,70],[-30,70],[-30,70],[-30,70],[-30,70],[-30,70]]
n_gauss = [[200,200],[200,200],[200,200],[200,200],[200,200],[400,400],[400,400],[400,400],[400,400],[400,400],[400,400],[400,400],[400,400]]
gamma_gauss = [[0,40],[0,40],[0,40],[0,40],[0,40],[0,40],[0,40],[0,40],[0,40],[0,40],[0,40],[0,40],[0,40]]


# station list input
stationlist_path = os.path.abspath("./stationlist_atlantic_26.csv")

# wavefield input
wavefield_path = os.path.abspath("../axisem/SOLVER/Noisi_10s_repacked/")


# Need to create a project for each, then a source, convert wavefield, then calculate correlations
for i in range(0,np.size(gridpoints_homo)):
    print('============== Homogeneous Grids =================')
    print('============= Working on Grid {} of {} ================'.format(i+1,np.size(gridpoints_homo)))
    
    # new project
    project_name = "BENCHMARK_Homo_" + str(gridpoints_homo[i])
    os.system('noisi setup_project ' + project_name)
    print('New project created:', project_name)
    project_path = os.path.join(main_path,project_name)
    print('Path to new project: ', project_path)
    
    # move to project directory
    os.chdir(project_path)
    print('Changed to project directory:',os.getcwd())
    
    # setup grid
    
    # Write to config.json file
    config_path = os.path.join(project_path,'config.json')
    
    with io.open(config_path,'r+') as fh:
        conf = json.loads(fh.read())
        
    conf['gauss_grid'] = True      
    conf['gauss_sigma'] = [sigma_homo[i]]
    conf['gauss_beta'] = [beta_homo[i]]
    conf['gauss_phi_ini'] = [phi_ini_homo[i]]
    conf['gauss_phi_max'] = [phi_max_homo[i]]
    conf['gauss_lat_0'] = [lat_0_homo[i]]
    conf['gauss_lon_0'] = [lon_0_homo[i]]
    conf['gauss_n'] = [n_homo[i]]
    conf['gauss_gamma'] = [gamma_homo[i]]
    conf['gauss_plot'] = False
    conf['gauss_dense_antipole'] = False
    conf['gauss_only_ocean'] = True

    # Set below to true if voronoi cell surface area is to be calculated
    conf['voronoi_surface_area'] = True

    # change instaseis and wavefield
    conf['instaseis'] = True
    conf['wavefield_path'] = os.path.abspath(wavefield_path)

    with io.open(config_path,'w') as fh:
        cf = json.dumps(conf,sort_keys=False, indent=4, separators=(",", ": "))
        fh.write(cf)

    # calculate grid with noisi
    print('Computing grid...')
    os.system('noisi setup_sourcegrid ' + project_path)
    print('Grid computed and saved as sourcegrid.npy')
    
    #plot and save sourcegrid
    sourcegrid_path = os.path.join(project_path,'sourcegrid.npy')
    grid = np.load(sourcegrid_path)

    plt.figure(figsize=(25,10))
    ax = plt.axes(projection=ccrs.Mollweide())
    ax.coastlines()
    plt.scatter(grid[0],grid[1],s=1,c='k',transform=ccrs.Geodetic())
    plt.title('Final grid with {} gridpoints'.format(np.size(grid[0])))
    plt.savefig(os.path.join(project_path,"sourcegrid.png"))
    plt.show(block=False)
    
    # Setup a homogeneous source called homo_source
    source_homo = "homo_source"

    source_homo_path = os.path.join(project_path,source_homo)
    os.system ('noisi setup_source ' + source_homo_path)

    print('New source created: ', source_homo_path)

    # Get number of stations for mpi
    stationlist = read_csv(stationlist_path)
    station_n = np.shape(stationlist)[0]
    print('Number of stations: ',station_n)
    
    
    # Convert wavefield 
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
        
    # setup noise source for homogeneous model

    print('Setting up noisesource distribution...')
    setup_noisesource_new(project_path,source_homo_path)
    print('Done.')
    
    
    # calculate correlations for homogeneous source
    print("Computing correlations...")
    os.system("mpirun -np 8 noisi correlation {} {}".format(source_homo_path,0))
    print("All correlations computed.")
    
    # plot correlations and save the file
    corr_homo_path = os.path.join(source_homo_path,"step_0/corr")

    # make azimuth selection
    traces_homo =plot_section(corr_homo_path,stationlist_path,bandpass = None,comp = 'BHZ',fmt = 'SAC',az_selection = [0,180], 
                          scale = 1., resol = 1,plot=False)
    
    #plot
    maxlag = (traces_homo[0].stats.npts-1) / 2.0
    
    traces_homo.plot(type='section',orientation='horizontal',
    reftime = traces_homo[0].stats.starttime + maxlag,scale=1.,outfile=os.path.join(project_path,'homo_correlations.png'))


    # Go back to main directory
    os.chdir(main_path)
    print('Back to main directory:',os.getcwd())
    
    
    
print('=======================================')  
print('======= HOMOGENEOUS GRIDS DONE ========')
print('=======================================')  

    
# Need to create a project for each, then a source, convert wavefield, then calculate correlations

for i in range(0,np.size(gridpoints_gauss)):
    print('============= Gaussian Grids ===============')
    print('============= Working on Grid {} of {} ================'.format(i+1,np.size(gridpoints_gauss)))
    
    # new project
    project_name = "BENCHMARK_Gauss_" + str(gridpoints_gauss[i])
    os.system('noisi setup_project ' + project_name)
    print('New project created:', project_name)
    project_path = os.path.join(os.getcwd(),project_name)
    print('Path to new project: ', project_path)
        
    # move to project directory
    os.chdir(project_path)
    print('Changed to project directory:',os.getcwd())
    
    # setup grid
    
    # Write to config.json file
    config_path = os.path.join(project_path,'config.json')
    
    with io.open(config_path,'r+') as fh:
        conf = json.loads(fh.read())
        
    conf['gauss_grid'] = True      
    conf['gauss_sigma'] = [sigma_gauss[i]]
    conf['gauss_beta'] = [beta_gauss[i]]
    conf['gauss_phi_ini'] = [phi_ini_gauss[i]]
    conf['gauss_phi_max'] = [phi_max_gauss[i]]
    conf['gauss_lat_0'] = [lat_0_gauss[i]]
    conf['gauss_lon_0'] = [lon_0_gauss[i]]
    conf['gauss_n'] = [n_gauss[i]]
    conf['gauss_gamma'] = [gamma_gauss[i]]
    conf['gauss_plot'] = False
    conf['gauss_dense_antipole'] = False
    conf['gauss_only_ocean'] = True

    # Set below to true if voronoi cell surface area is to be calculated
    conf['voronoi_surface_area'] = True

    # change instaseis and wavefield
    conf['instaseis'] = True
    conf['wavefield_path'] = os.path.abspath(wavefield_path)

    with io.open(config_path,'w') as fh:
        cf = json.dumps(conf,sort_keys=False, indent=4, separators=(",", ": "))
        fh.write(cf)

    # calculate grid with noisi
    print('Computing grid...')
    os.system('noisi setup_sourcegrid ' + project_path)
    print('Grid computed and saved as sourcegrid.npy')
    
    #plot and save sourcegrid
    sourcegrid_path = os.path.join(project_path,'sourcegrid.npy')
    grid = np.load(sourcegrid_path)

    plt.figure(figsize=(25,10))
    ax = plt.axes(projection=ccrs.Mollweide())
    ax.coastlines()
    plt.scatter(grid[0],grid[1],s=1,c='k',transform=ccrs.Geodetic())
    plt.title('Final grid with {} gridpoints'.format(np.size(grid[0])))
    plt.savefig(os.path.join(project_path,"sourcegrid.png"))
    plt.show(block=False)
    
    # Setup a homogeneous source called homo_source
    source_homo = "homo_source"

    source_homo_path = os.path.join(project_path,source_homo)
    os.system ('noisi setup_source ' + source_homo_path)

    print('New source created: ', source_homo_path)

    # Get number of stations for mpi
    stationlist = read_csv(stationlist_path)
    station_n = np.shape(stationlist)[0]
    print('Number of stations: ',station_n)
    
    
    # Convert wavefield 
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
        
    # setup noise source for homogeneous model

    print('Setting up noisesource distribution...')
    setup_noisesource_new(project_path,source_homo_path)
    print('Done.')
    
    
    # calculate correlations for homogeneous source
    print("Computing correlations...")
    os.system("mpirun -np 8 noisi correlation {} {}".format(source_homo_path,0))
    print("All correlations computed.")
    
    # plot correlations and save the file
    corr_homo_path = os.path.join(source_homo_path,"step_0/corr")

    # make azimuth selection
    traces_homo =plot_section(corr_homo_path,stationlist_path,bandpass = None,comp = 'BHZ',fmt = 'SAC',az_selection = [0,180], 
                          scale = 1., resol = 1,plot=False)
    #plot
    maxlag = (traces_homo[0].stats.npts-1) / 2.0
    traces_homo.plot(type='section',orientation='horizontal',
    reftime = traces_homo[0].stats.starttime + maxlag,scale=1.,outfile=os.path.join(project_path,'homo_correlations.png'))


    
    # Go back to main directory
    os.chdir(main_path)
    print('Back to main directory:',os.getcwd())
    
print('=======================================')  
print('======== GAUSSIAN GRIDS DONE ==========')
print('=======================================')     
    
    
print('======== BENCHMARK DONE =========')
