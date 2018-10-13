# automatically setup gaussian grid with Ardhuin data

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from noisi.scripts.source_grid_gauss import spherical_distance_degrees


def setup_data_gaussgrid(project_path,data_path,t_data,f_dom,data_thresh,gamma_thresh,gamma_const,back_grid_phi,data_grid_phi,extent = None, plot=False):
    """
    This function automatically sets up the gaussian grid
    
    Input: 
    project_path, data_path
    
    t_data :: set index for time step in data (check with panoply)
    f_dom :: all frequencies up to dominant frequency f_dom are summed#
    data_thresh :: dense grid will be added for points above this threshold. Data is normalised so 1 is max. 
    gamma_thresh :: radius in which points will be removed. (e.g. 10)
    back_grid_phi :: background gridpoint distance in degrees (e.g. 4)
    data_grid_phi :: additional data gridpoint distance in degrees (e.g. 0.5)
    gamma_const :: multiplies gamma_thresh with a constant which is then the gamma (circle radius in degrees which will be replaced) for the additional grid (e.g. 3)
    
    extent :: [lon_min,lon_max,lat_min,lat_max] gives area in which grid will be made denser
    plot :: True/False
    
    Output:
    sigma,beta,phi_ini,phi_max,lat_0,lon_0,n,gamma 
    """
    
    data = Dataset(data_path,'r')
    #p2l = data.variables['p2l'][t_data,f_data,:,:] # pick frequency band and time here
    longitude = data.variables['longitude'][:]
    latitude = data.variables['latitude'][:]
    f = data.variables['f'][:]
    
    
    # sum everything up to dominant period of wavefield
    #f_dom = 1/10
    f_max = np.argmin(np.abs(f-f_dom))
    
    p2l_freq = data.variables['p2l'][t_data,:f_max,:,:]
    
    #p2l_freq_sum = np.sum(p2l_freq,axis=0)
    #p2l_freq_10 = 10**p2l_freq_sum
    
    # Data is converted before it's summed, this way the maximum value isn't as big
    p2l_freq_sum_10 = np.sum(10**p2l_freq,axis=0)  
    
    # Turn latitude and longitude into one grid
    grid_data_lon = []
    grid_data_lat = []
    grid_data_p2l = []
    
    for i in range(0,np.size(latitude)):
        for j in range(0,np.size(longitude)):
            grid_data_lon.append(longitude[j])
            grid_data_lat.append(latitude[i])
            grid_data_p2l.append(p2l_freq_sum_10[i,j])
    
            
    #grid_data = [grid_data_lon,grid_data_lat,grid_data_p2l]
    
    # let's normalise the data
    grid_data_p2l_norm = np.nan_to_num(grid_data_p2l)/np.nanmax(grid_data_p2l)
    grid_data_norm = [grid_data_lon,grid_data_lat,grid_data_p2l_norm]

    # plot
    if plot:
        plt.figure(figsize=(25,10))
        ax = plt.axes(projection=ccrs.PlateCarree())
        plt.scatter(grid_data_lon, grid_data_lat,s=1, c=grid_data_p2l_norm, transform=ccrs.PlateCarree(),cmap=plt.get_cmap('jet'),vmax=data_thresh)
        plt.colorbar()
        plt.title('Dataset_path up to frequency {} Hz'.format(np.around(f[f_max],decimals=5)))
        ax.coastlines()
        plt.show()
    
    

    lat_ini = []
    lon_ini = []
    
    # iterate over all grid points
    # data is in grid_data
    
    for i in range(0,np.size(grid_data_norm[2])):
        if grid_data_norm[2][i] >= data_thresh:
            lat_ini.append(grid_data_norm[1][i])
            lon_ini.append(grid_data_norm[0][i])
        else:
            continue
            
    print("Number of gridpoints above threshold: ", np.size(lat_ini))
    
    
    ## Remove all points above threshold that are not within area
    if extent is not None:
        lon_min = extent[0]
        lon_max = extent[1]
        lat_min = extent[2]
        lat_max = extent[3]
        
        lat_ini_cut = []
        lon_ini_cut = []

        for i in range(0,np.size(lat_ini)):
            if lat_min < lat_ini[i] < lat_max and lon_min < lon_ini[i] < lon_max:
                lat_ini_cut.append(lat_ini[i])
                lon_ini_cut.append(lon_ini[i])
            else:
                continue
                
        lat_ini = lat_ini_cut
        lon_ini = lon_ini_cut

    
    if plot:
        plt.figure(figsize=(25,10))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.coastlines()
        ax.set_extent([-180,180,-90,90])
        plt.scatter(lon_ini,lat_ini,s=1,c='k')
        plt.show()
        
    # Now remove grids that are not really necessary
    
    lat_new = []
    lon_new = []
    
    lat_var = []
    lon_var = []
    
    j = 0
    
    while j < np.size(lat_ini):
        if j == 0:
            for i in range(0,np.size(lat_ini)):
                lat_var.append(lat_ini[0])
                lon_var.append(lon_ini[0])
                dist_var = spherical_distance_degrees(lat_ini[0],lon_ini[0],lat_ini[i],lon_ini[i])
                if dist_var > gamma_thresh:
                    lat_var.append(lat_ini[i])
                    lon_var.append(lon_ini[i])
                else:
                    continue
        else:
            try:
                lat_new = lat_var
                lon_new = lon_var
                lat_var = []
                lon_var = []
                lat_var.append(lat_new[j])
                lon_var.append(lon_new[j])
                for i in range(0,np.size(lat_new)):
                    dist_var = spherical_distance_degrees(lat_new[j],lon_new[j],lat_new[i],lon_new[i])
                    if dist_var > gamma_thresh:
                        lat_var.append(lat_new[i])
                        lon_var.append(lon_new[i])
                    else:
                        continue
            except:
                break          
        j += 1
        
                
             
    print("Number of additional grids:", np.size(lat_new))
    
    if plot:
        plt.figure(figsize=(25,10))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.coastlines()
        ax.set_extent([-180,180,-90,90])
        plt.scatter(lon_new,lat_new,s=1,c='k')
        plt.show()
        
    # Now create variables for gaussian grids
    
    # Background grid
    sigma = [20]
    beta = [5]
    phi_ini = [back_grid_phi]
    phi_max = [back_grid_phi]
    lat_0 = [0]
    lon_0 = [0]
    n = [200]
    gamma = [0]
    #plot = False
    #dense_antipole = False
    #only_ocean = True
    
    for i in range(0,np.size(lat_new)):
        lat_0.append(lat_new[i])
        lon_0.append(lon_new[i])
        sigma.append(20)
        beta.append(5)
        phi_ini.append(data_grid_phi)
        phi_max.append(data_grid_phi)
        n.append(200)
        gamma.append(gamma_const*gamma_thresh)
        
    
    return sigma,beta,phi_ini,phi_max,lat_0,lon_0,n,gamma
