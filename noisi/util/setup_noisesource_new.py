# New setup_sourcegrid script   
import numpy as np
from obspy.geodetics import gps2dist_azimuth
from obspy.signal.invsim import cosine_taper
import matplotlib.pyplot as plt
import h5py
from noisi.my_classes.wavefield import WaveField
import json
from glob import glob
import os
from scipy.fftpack import next_fast_len
from netCDF4 import Dataset
#from obspy.core import UTCDateTime
import cartopy.crs as ccrs
#from noisi.util import plot

def setup_noisesource_new(project_path,source_path, data_path = None,t_data= None,f_dom = None,outfile = False):    

    ######################
    # Define the bases!
    ######################
    # geography - Add anything else than a homogeneous distribution by setting to "True" the following:
    only_ocean = True
    gaussian_blobs = False
    params_gaussian_blobs = [{'center':(0.,0.),'sigma_radius_m':500000.,'rel_weight':2.}]
    
    #spectra
    params_gaussian_spectra = [{'central_freq':0.075,'sigma_freq':0.1,'weight':10.}]
                              #{'central_freq':0.075,'sigma_freq':0.1,'weight':10.}]
                              
            
    sourcegrid_path = os.path.join(project_path,'sourcegrid.npy')              
    
    grd  = np.load(sourcegrid_path)
    ntraces = np.shape(grd)[-1]
    print("Number of gridpoints: ", ntraces)
    
    if data_path is not None:
        data = Dataset(data_path,'r')
        #p2l = data.variables['p2l'][t_data,f_data,:,:] # pick frequency band and time here
        longitude = data.variables['longitude'][:]
        latitude = data.variables['latitude'][:]
        f = data.variables['f'][:]
        
        
        # sum everything up to dominant period of wavefield
        f_max = np.argmin(np.abs(f-f_dom))
        
        p2l_freq = data.variables['p2l'][t_data,f_max-2:f_max+1,:,:]
        
        #p2l_freq_sum = np.sum(p2l_freq,axis=0)
        #p2l_freq_10 = 10**p2l_freq_sum
        p2l_freq_sum_10 = np.sum(10**p2l_freq,axis=0)  # Data is converted before it's summed, this way the maximum value isn't as big
        
        # Turn latitude and longitude into one grid
        grid_data_lon = []
        grid_data_lat = []
        grid_data_p2l = []
        
        for i in range(0,np.size(latitude)):
            for j in range(0,np.size(longitude)):
                grid_data_lon.append(longitude[j])
                grid_data_lat.append(latitude[i])
                grid_data_p2l.append(p2l_freq_sum_10[i,j])
        
                
        grid_data = [grid_data_lon,grid_data_lat,grid_data_p2l]
        
    
        # create frequency spectrum
        params_gaussian_spectra = [{'central_freq':f[f_max],'sigma_freq':0.1,'weight':10.}]
                                  #{'central_freq':0.075,'sigma_freq':0.1,'weight':10.}]
            
        plt.figure(figsize=(25,10))
        ax = plt.axes(projection=ccrs.PlateCarree())
        plt.contourf(longitude, latitude, p2l_freq_sum_10,60, transform=ccrs.PlateCarree(),cmap=plt.get_cmap('jet'))
        plt.colorbar()
        plt.title('Dataset_path transformed up to frequency {} Hz'.format(np.around(f[f_max],decimals=5)))
        ax.coastlines()
        plt.show(block=False)
    
    
    config = json.load(open(os.path.join(project_path,'config.json')))    
    source_config = json.load(open(os.path.join(source_path,'source_config.json')))
    # Find the sampling rate and the nr. of points in synthetic seismograms
    synthetics_path = config['wavefield_path']
    
    if source_config['preprocess_do']:
        ext = '*.h5_proc'
    else:
        ext = '*.h5'
    #print(ext)
    wfs = glob(os.path.join(synthetics_path,ext))
    with WaveField(wfs[0]) as wf:
        df = wf.stats['Fs']
        nt = wf.stats['nt']
        # The number of points for the fft is larger due to zeropadding --> apparent higher frequency sampling\n",
        n = next_fast_len(2*nt-1)
        print('df =',df,'nt = ',nt,'n = ',n)
        freq = np.fft.rfftfreq(n,d=1./df)
        print('freq_min = ',np.min(freq),'freq_max = ',np.max(freq))
        taper = cosine_taper(len(freq),0.05)
        #print(len(freq))
        
        
        
        
    def get_distance(grid,location):
        def f(lat,lon,location):
            return abs(gps2dist_azimuth(lat,lon,location[0],location[1])[0])
        dist = np.array([f(lat,lon,location) for lat,lon in zip(grid[1],grid[0])])
        return dist
        # Use cartopy to figure out where ocean is
    def get_ocean_mask():
        print('Getting ocean mask...')
        from noisi.util.cartopy_is_land import is_land
        print("Latitude from {} to {},\n\
        Longitude from {} to {}".format(
        round(np.min(grd[1]),2),
        round(np.max(grd[1]),2),
        round(np.min(grd[0]),2),
        round(np.max(grd[0]),2)))
        ocean_mask = [not is_land(lon,lat) for (lon,lat) in grd.transpose()]
        
        return np.array(ocean_mask)
    
    
    
    #########################
    # Create the source distr
    #########################
    
    # geography
    num_bases = 1
    if gaussian_blobs:
        num_bases += len(params_gaussian_blobs)
    
    basis1 = np.zeros((num_bases,ntraces))
    #print(np.shape(basis1))
    # homogeneous layer
    basis1[0,:] = np.ones(ntraces) 
    if only_ocean:
        basis1[0,:] *= np.array(get_ocean_mask()).astype(int)
        # superimposed Gaussian blob(s)
    if gaussian_blobs:
        i = 1
        for blob in params_gaussian_blobs:
            dist = get_distance(grd,blob['center'])
            basis1[i,:] = np.exp(-(dist)**2/(2*blob['sigma_radius_m']**2))
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
    
    #print(np.shape(basis1))
    
    
    
    ######################
    # set the weights
    #####################
    # geography
    weights1 = np.ones(np.shape(basis1)[0])
    #print(np.shape(weights1))
    if gaussian_blobs:
        i = 1
        for blob in params_gaussian_blobs:
            weights1[i] = blob['rel_weight']
            i+=1
    #print weights1
    # spectra --- much harder to assign manually, since we need weights for every location. just assigning ones.\n",
    weights2 = np.ones((np.shape(grd)[-1],np.shape(basis2)[0]))
    i=0
    for spec in params_gaussian_spectra:
        weights2[:,i] *= spec['weight']
    #print weights2
    
    if data_path is not None:
         # Sample grid onto Gauss Grid
        print('Sampling grid...')
        
        p2l_new = []
        lat_new = []
        lon_new = []
        
        for k in range(0,np.size(grd[0])):
        
            if k%1000 == 0:
                print('At Gridpoint: ',k)
            dist_var = np.sqrt((np.array(grid_data[1])-grd[1][k])**2+(np.array(grid_data[0])-grd[0][k])**2)
        
            # Append interpolated grid to new variables
            p2l_new.append(grid_data_p2l[np.nanargmin(dist_var)])
            lat_new.append(grd[1][k])
            lon_new.append(grd[0][k])
        
        print('Done.')
        
        # Get rid of all nan's in the data
        p2l_final = np.nan_to_num(p2l_new)
        
        grd = np.array([grd[0],grd[1]])
        with h5py.File('grid_sampled.h5','w') as fh:
            fh.create_dataset('coordinates',data=grd.astype(np.float32))
            fh.create_dataset('data',data=p2l_final.astype(np.float32))       
    
        basis1 = np.array([p2l_final])
    
    
    distr = np.dot(weights1,basis1)
    #plot.plot_grid(grd[0],grd[1],distr,mode='srcdots')
    plt.figure(figsize=(25,10))
    ax = plt.axes(projection=ccrs.PlateCarree())
    plt.scatter(grd[0], grd[1],s=1, c=distr, transform=ccrs.PlateCarree(),cmap=plt.get_cmap('jet'))
    plt.colorbar()
    if data_path is not None:
        plt.title('Noisesource distribution for t = {} of {} with frequencies summed up to {}'.format(t_data,data_path,f_dom))
    else:
        plt.title('Noisesource distribution for {}'.format(source_path))
    ax.coastlines()
    if outfile:
        plt.savefig(os.path.join(source_path,'noisesource_dist.png'))
    plt.show(block=False)
    
    plt.semilogx(freq,np.dot(weights2[0,:],basis2))
    plt.title('Frequency spectrum of noise source distribution')
    plt.show(block=False)
    
    
    ##################
    # Calculate Voronoi cell surface area
    ##################
    
    # Set config["voronoi_surface_area"]: true if voronoi surface areas are to be calculated
    # Otherwise all surface areas will be set to one. 

    if config['voronoi_surface_area']:
        from noisi.util.voronoi_surface_area import get_voronoi_surface_area
    
        if only_ocean:
            from noisi.scripts.source_grid_gauss import gauss_grid
            
            n_grids = np.size(config['gauss_sigma'])
              
            only_ocean_false = False
                
            plot_false = False
    
            # compute full grid with only_ocean = False
            grd_full = gauss_grid(config['gauss_sigma'],config['gauss_beta'],config['gauss_phi_ini'],config['gauss_phi_max'],
                          config['gauss_lat_0'],config['gauss_lon_0'],config['gauss_n'],config['gauss_gamma'],
                          plot=plot_false,dense_antipole=config['gauss_dense_antipole'],
                          only_ocean=only_ocean_false)
            
            # Calculate voronoi cells for the whole grid
            grd, voronoi_areas = get_voronoi_surface_area(grd_full)
            
            # Plot
            plt.figure(figsize=(25,10))
            ax = plt.axes(projection=ccrs.Mollweide())
            ax.coastlines()
            plt.scatter(grd[0],grd[1],s=1,c='k',transform=ccrs.Geodetic())
            plt.title('Intial Gauss Grid with all gridpoints')
            plt.show(block=False) 
    
            # Remove points on land again
            from noisi.util.cartopy_is_land import is_land
            grid_onlyocean_lon = []
            grid_onlyocean_lat = []
            voronoi_areas_onlyocean = []
    
            for i in range(0,np.size(grd[0])):
                if not is_land(grd[0][i],grd[1][i]):
                    grid_onlyocean_lon.append(grd[0][i])
                    grid_onlyocean_lat.append(grd[1][i])
                    voronoi_areas_onlyocean.append(voronoi_areas[i])
                else:
                    continue
             
            plt.figure(figsize=(25,10))
            ax = plt.axes(projection=ccrs.Mollweide())
            ax.coastlines()
            plt.scatter(grid_onlyocean_lon,grid_onlyocean_lat,s=1,c=voronoi_areas_onlyocean,transform=ccrs.Geodetic(),cmap=plt.get_cmap('seismic_r'))
            plt.title('Final Gauss Grid with gridpoints removed')
            plt.colorbar()
            plt.show(block=False) 
            
            print('Gridpoints and voronoi cells on land removed.')
            grd = np.asarray([grid_onlyocean_lon,grid_onlyocean_lat])
            surf_areas = np.asarray(voronoi_areas_onlyocean)
            print('Final number of gridpoints:', int(np.size(grd)/2))
    
        else:
            grd, voronoi_areas = get_voronoi_surface_area(grd)
            surf_areas = voronoi_areas
    else:
        surf_areas = np.ones(ntraces)


    # Export to files
    # Save to an hdf5 file
    
    with h5py.File(os.path.join(source_path,'sourcemodel.h5'),'w') as fh:
        fh.create_dataset('coordinates',data=grd.astype(np.float32))
        fh.create_dataset('frequencies',data=freq.astype(np.float32))
        fh.create_dataset('distr_basis',data=basis1.astype(np.float32))
        fh.create_dataset('distr_weights',data=weights1.astype(np.float32))
        fh.create_dataset('spect_basis',data=basis2.astype(np.float32))
        fh.create_dataset('spect_weights',data=weights2.astype(np.float32))
        fh.create_dataset('surf_areas',data=surf_areas.astype(np.float32))
        
        
    # Save sourcemodel as starting_model
    
    with h5py.File(os.path.join(source_path,'step_0','starting_model.h5'),'w') as fh:
        fh.create_dataset('coordinates',data=grd.astype(np.float32))
        fh.create_dataset('frequencies',data=freq.astype(np.float32))
        fh.create_dataset('distr_basis',data=basis1.astype(np.float32))
        fh.create_dataset('distr_weights',data=weights1.astype(np.float32))
        fh.create_dataset('spect_basis',data=basis2.astype(np.float32))
        fh.create_dataset('spect_weights',data=weights2.astype(np.float32))
        fh.create_dataset('surf_areas',data=surf_areas.astype(np.float32))
    
    
        
    # Create base_model
    
    basis1_b = np.ones(basis1.shape)
    with h5py.File(os.path.join(source_path,'step_0','base_model.h5'),'w') as fh:
        fh.create_dataset('coordinates',data=grd.astype(np.float32))
        fh.create_dataset('frequencies',data=freq.astype(np.float32))
        fh.create_dataset('distr_basis',data=basis1_b.astype(np.float32))
        fh.create_dataset('distr_weights',data=weights1.astype(np.float32))
        fh.create_dataset('spect_basis',data=basis2.astype(np.float32))
        fh.create_dataset('surf_areas',data=surf_areas.astype(np.float32))



    
    
    

 
        
        
        
