import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import json
import os
import io
import warnings
warnings.filterwarnings("ignore")

def gauss_grid_one(sigma,beta,phi_ini,phi_max,lat_0,lon_0,n,plot=True,dense_antipole = True,only_ocean=False):
    """
    This function creates a Gaussian grid. Input parameters:
    :sigma (greater than 2) = standard deviation, i.e. size of the area of denser grid points
    :beta = steepness of the drop of to the maximum distance
    :phi_ini = initial distance between grid points, in degrees
    :phi_max = maximum distance between grid points, in degrees
    :lat_0 = latitude of point of interest, -90° to 90°
    :lon_0 = longitude of point of interest, -180° to 180°
    :n = number of circles
    :plot = True/False
    :dense_antipole = True/False
    :only_ocean = True/False
    Returns: list of longitudes and latitudes
    """
    
    # Error messages
    if lat_0 < -90 or lat_0 > 90:
        msg = 'lat_0 has to be between -90° and 90°'
        raise ValueError(msg)
    
    if lon_0 < -180 or lon_0 > 180:
        msg = 'lon_0 has to be between -180 and 180'
        raise ValueError(msg)
        
    if phi_ini > 90 or phi_max > 90:
        msg = 'phi_ini and phi_max should not be larger than 90'
        raise ValueError(msg)
    
    if phi_ini > phi_max:
        msg = 'phi_ini should be smaller than phi_max'
        raise ValueError(msg)

    # Step 1: Gauss
    # Calculate radii of the circles
    phi_max = phi_max - phi_ini
    
    lat = np.linspace(0,90,n)
    dphi1 = phi_max*(1 - np.exp(-(lat/sigma)**beta))+phi_ini

    phi = []
    dphi = []
    phi_0 = 0
    
    if dense_antipole:
        for i in range(0,np.size(dphi1)):
            phi_0 += dphi1[i]
            phi.append(phi_0)
            dphi.append(dphi1[i])
            # Change condition so that if the distance between equator and previous circle is greater than that before the point is removed
            if phi_0 > 90:
                if dphi[i] > dphi[i-1]:
                    if 90-phi[i-1] < dphi[i-1]:
                        phi = phi[:-2]  # removes last two entries of phi
                        dphi = dphi[:-2] # removes last two dphi
                        phi_0 = 90
                        phi.append(phi_0)
                        dphi.append(90-phi[i-2])
                        break
                    else:
                        phi = phi[:-1]  # removes last entry of phi since it would be bigger than 90
                        dphi = dphi[:-1] # removes last phi
                        break
                elif dphi[i] <= dphi[i-1]: 
                    if 90-phi[i-1] < dphi[i-1]:
                        phi = phi[:-2]  # removes last two entries of phi
                        dphi = dphi[:-2] # removes last two dphi
                        phi_0 = 90
                        phi.append(phi_0)
                        dphi.append(90-phi[i-2])
                        break
                    else:
                        phi = phi[:-1]  # removes last entry of phi since it would be bigger than 90
                        dphi = dphi[:-1] # removes last phi
                        break       
                else:
                    phi = phi[:-1]  # removes last entry of phi since it would be bigger than 90
                    dphi = dphi[:-1] # removes last phi
                    break
    else:
        for i in range(0,np.size(dphi1)):
            phi_0 += dphi1[i]
            phi.append(phi_0)
            dphi.append(dphi1[i])            
            if phi_0 > 180:
                phi = phi[:-1]  # removes last entry of phi since it would be bigger than 180
                dphi = dphi[:-1] # removes last phi
                break                

    # phi now consists of the latitudinal angles up to 90° over which we should loop

    # Step 2: Longitudinal angles (azmiuth)
    # We want these angles to be close to the latitudinal distance between the circles

    dtheta = dphi
    N = np.around(np.divide(360,dtheta))
    
    # To get the angle we now use 2*Pi/N
    theta = np.divide(2*np.pi,N)
    
    ## We now have the latitudes and the angle over which we have to loop to get the longitudes. 
    # Step 3: Loop

    lat_final1 = [0]
    lon_final1 = [0]
    

    for i in range(0,np.size(phi)-1):
        for j in range(0,int(N[i])):
            # first need to calculate all longitudes
            lon_final1.append(np.rad2deg(j*theta[i]))
            lat_final1.append(phi[i])
            
    # need to shift it
    lon_final1 = np.subtract(lon_final1,180)
    lat_final1 = np.subtract(lat_final1,90)
    
    
    if dense_antipole:
        # Now flip the longitude to make the other hemisphere
        lat_final2 = [0]
        lon_final2 = [0] 

        for i in range(0,np.size(phi)):    # size(phi) - 1 since we don't want two sets of points around the equator
            for j in range(0,int(N[i])):
                # first need to calculate all longitudes
                lon_final2.append(np.rad2deg(j*theta[i]))
                lat_final2.append(phi[i])

        # Shift coordinates and flip it
        lon_final2 = np.subtract(lon_final2,180)
        lat_final2 = np.subtract(lat_final2,90)
        lat_final2 = np.multiply(lat_final2,-1)
        lon_final2 = lon_final2

        # Combine the two
        lat_final = np.concatenate((lat_final1,lat_final2))
        lon_final = np.concatenate((lon_final1,lon_final2))
    else:
        lat_final = lat_final1
        lon_final = lon_final1
    
    # Calculate grid point distance of densest grid area in m for latitude
    dlat = 111132.954 - 559.822 * np.cos(2*lat_0) + 1.175*np.cos(4*lat_0)
    dx_min = dphi[0]*dlat
    dx_max = dphi[-1]*dlat
        
    #print('Number of Gridpoints:',np.size(lat_final))
    print('Minimum dx in m:',np.around(dx_min,3),'m which is',np.around(dphi[0],3),'°')
    print('Maximum dx in m:',np.around(dx_max,3),'m which is',np.around(dphi[-1],3),'°')

    # ROTATION TO AREA OF INTEREST
    # We have the variables lon_final, lat_final
    # lat_0 and lon_0 give the point to which the center should be shifted in degrees

    # Step 1: Convert lat_final & lon_final from degrees to radians
    if dense_antipole:
        lat_final_rad = np.deg2rad(lat_final)
        lon_final_rad = np.deg2rad(lon_final)
    else:
        lat_final_rad = -np.deg2rad(lat_final)
        lon_final_rad = -np.deg2rad(lon_final)

    # Rotation around y-axis and z-axis
    theta_rot = np.deg2rad(90-lat_0)
    phi_rot = np.deg2rad(lon_0)

    # Convert grid from spherical to cartesian 
    x_cart = np.cos(lon_final_rad)*np.cos(lat_final_rad)
    y_cart = np.sin(lon_final_rad)*np.cos(lat_final_rad)
    z_cart = np.sin(lat_final_rad)

    # Shift the coordinates, it's the two rotation matrices for z and y multiplied.
    x_cart_new = np.cos(theta_rot)*np.cos(phi_rot)*x_cart + np.cos(theta_rot)*np.sin(phi_rot)*y_cart + np.sin(theta_rot)*z_cart
    y_cart_new = -np.sin(phi_rot)*x_cart + np.cos(phi_rot)*y_cart
    z_cart_new = -np.sin(theta_rot)*np.cos(phi_rot)*x_cart - np.sin(theta_rot)*np.sin(phi_rot)*y_cart + np.cos(theta_rot)*z_cart

    # Convert cartesian back to spherical coordinates
    # Brute force for longitude because rotation matrix does not seem to rotate the longitudes
    lon_final_rot = []
    
    for i in range(0,np.size(lon_final_rad)):
        lon_final_rot.append(np.rad2deg(np.arctan2(y_cart_new[i],x_cart_new[i])+phi_rot))
        if lon_final_rot[i] > 180:
            lon_final_rot[i] = lon_final_rot[i] - 360
        elif lon_final_rot[i] < -180:   
            lon_final_rot[i] = lon_final_rot[i] + 360
            
    lat_final_rot = np.rad2deg(np.arcsin(z_cart_new))

    
    if only_ocean:
        from noisi.util.cartopy_is_land import is_land
        grid_onlyocean_lon = []
        grid_onlyocean_lat = []

        for i in range(0,np.size(lat_final_rot)):
            if not is_land(lon_final_rot[i],lat_final_rot[i]):
                grid_onlyocean_lon.append(lon_final_rot[i])
                grid_onlyocean_lat.append(lat_final_rot[i])
            else:
                continue
        
        lon_final_final = grid_onlyocean_lon
        lat_final_final = grid_onlyocean_lat
    else:
        lon_final_final = lon_final_rot
        lat_final_final = lat_final_rot
        
    print('Number of gridpoints:',np.size(lat_final_final))
    
    if plot:
        plt.figure(figsize=(25,10))
        ax = plt.axes(projection=ccrs.Mollweide())
        ax.coastlines()
        plt.scatter(lon_final_final,lat_final_final,s=1,c='k',transform=ccrs.Geodetic())
        plt.title('Centre at %0.2f ° latitude and %0.2f ° longitude with %i gridpoints' %(lat_0,lon_0,np.size(lat_final_final)))
        plt.show()

        
    return list((lon_final_final,lat_final_final))


def spherical_distance_degrees(lat1,lon1,lat2,lon2):
    """
    Calculate ths distance between two points in degrees
    Input: latitude 1, longitude 1, latitude 2, longitude 2
    """
    
    dlat = lat2 - lat1
    dlon = lon2 - lon1

    hav = np.sin(np.deg2rad(dlat)/2)**2 + np.cos(np.deg2rad(lat1)) * np.cos(np.deg2rad(lat2)) * np.sin(np.deg2rad(dlon)/2)**2
    dist_deg = np.rad2deg(2 * np.arctan2(np.sqrt(hav), np.sqrt(1-hav)))

    return dist_deg


def gauss_grid(sigma,beta,phi_ini,phi_max,lat_0,lon_0,n,gamma,plot,dense_antipole,only_ocean):
    """
    This function creates several gaussian grids. To turn them into one grid the variable gamma has to be defined. 
    The input should be arrays with the variables for the gaussian grid plus gamma.
    The first grid does not need a gamma, but a random value should be given. This grid will be used as the base grid.
    Example input: [sigma1,sigma2],[beta1,beta2],[phi_ini1,phi_ini2],....,[gamma1,gamma2],[True,True],[False,True],[True,True]
    
    CAREFUL: When additional grids overlap, gridpoints are NOT removed. 
    
    :sigma,beta,phi_ini,phi_max,lat_0,lon_0,n,gamma = array of integers
    :plot,dense_antipole,only_ocean = array of True/False
    """
    
    # first check the number of grids given
    n_grids = np.size(sigma)
    
    if n_grids == 1:
        plot = True
        grid = gauss_grid_one(sigma,beta,phi_ini,phi_max,lat_0,lon_0,n,plot,dense_antipole,only_ocean)
        final_grid_lon = grid[0]
        final_grid_lat = grid[1]
    else:
        
        all_grids = []

        # Compute the different grids
        for i in range(0,n_grids):
            print('GRID: ',i+1)
            grid_one = gauss_grid_one(sigma[i],beta[i],phi_ini[i],phi_max[i],lat_0[i],lon_0[i],n[i],
                                  plot[i],dense_antipole[i],only_ocean[i])
            all_grids.append(grid_one)
            # different grids are now stored in all_grids


        final_grid_lon = []
        final_grid_lat = []

        # loop over all grids
        ## this only really works correctly if all grids have the same length otherwise it doesn't iterate over all points
        for i in range(0,np.size(all_grids[0][0])):
            if n_grids == 2:
                dist_deg = []
                dist_deg_var = spherical_distance_degrees(lat_0[1],lon_0[1],all_grids[0][1][i],all_grids[0][0][i])
                dist_deg.append(dist_deg_var)

                if dist_deg[0] <= gamma[1]:
                    continue
                else:
                    final_grid_lat.append(all_grids[0][1][i])
                    final_grid_lon.append(all_grids[0][0][i])

            else:
                dist_deg = []
                leave_out = False # to skip the points that are in the areas

                for j in range(0,n_grids-1):
                    dist_deg_var = spherical_distance_degrees(lat_0[j+1],lon_0[j+1],all_grids[0][1][i],all_grids[0][0][i])
                    dist_deg.append(dist_deg_var)
                    if dist_deg_var <= gamma[j+1]:
                        leave_out = True
                    else:
                        continue

                if leave_out:
                    continue
                else:
                    final_grid_lat.append(all_grids[0][1][i])
                    final_grid_lon.append(all_grids[0][0][i])

        plt.figure(figsize=(25,10))
        ax = plt.axes(projection=ccrs.Mollweide())
        ax.coastlines()
        plt.scatter(final_grid_lon,final_grid_lat,s=1,c='k',transform=ccrs.Geodetic())
        plt.title('Gauss Grid with gridpoints removed')
        plt.show() 


        # Now loop over the non-primary grids and add them to the grid
        # This seems to work
        dist_deg = []
        if n_grids == 2:
            for i in range(0,np.size(all_grids[1][0])):
                dist_deg = spherical_distance_degrees(lat_0[1],lon_0[1],all_grids[1][1][i],all_grids[1][0][i])
                if dist_deg <= gamma[1]:
                    final_grid_lat.append(all_grids[1][1][i])
                    final_grid_lon.append(all_grids[1][0][i])
                else:
                    continue
        else:
            for j in range(1,n_grids):
                for i in range(0,np.size(all_grids[j][0])):
                    dist_deg = spherical_distance_degrees(lat_0[j],lon_0[j],all_grids[j][1][i],all_grids[j][0][i])
                    if dist_deg <= gamma[j]:
                        final_grid_lat.append(all_grids[j][1][i])
                        final_grid_lon.append(all_grids[j][0][i])
                    else:
                        continue


        # plot
        plt.figure(figsize=(25,10))
        ax = plt.axes(projection=ccrs.Mollweide())
        ax.coastlines()
        plt.scatter(final_grid_lon,final_grid_lat,s=1,c='k',transform=ccrs.Geodetic())
        plt.title('Final grid with {} gridpoints'.format(np.size(final_grid_lon)))
        plt.show()
        
        print('Final number of gridpoints:',np.size(final_grid_lon))

    return list((final_grid_lon,final_grid_lat))
    


def create_sourcegrid_gauss(config):
        
    cfile = open(config,'r')
    config = json.load(cfile)
    cfile.close()
    # ToDo: Pretty printing of all dictionaries such as this one.
    print(config)
    
    #ToDo read extra parameters into configuration
    grid = gauss_grid(config['gauss_sigma'],config['gauss_beta'],config['gauss_phi_ini'],config['gauss_phi_max'],
                      config['gauss_lat_0'],config['gauss_lon_0'],config['gauss_n'],gamma=config['gauss_gamma'],
                      plot=config['gauss_plot'],dense_antipole=config['gauss_dense_antipole'],
                      only_ocean=config['gauss_only_ocean'])
   
    sources = np.zeros((2,len(grid[0])))
    #sources[0,:] = ids
    sources[0:2,:] = grid
    return sources


def setup_sourcegrid_gauss(configfile):
        
    sourcegrid = create_sourcegrid_gauss(configfile)
    
    with io.open(configfile,'r') as fh:
        config = json.load(fh)
    grid_filename = os.path.join(config['project_path'],'sourcegrid.npy')    
    
    # write to .npy    
    np.save(grid_filename,sourcegrid)
    print('Sourcegrid saved as sourcegrid.npy')
