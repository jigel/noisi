import numpy as np
import pandas as pd
from noisi.borrowed_functions.voronoi_polygons import getVoronoiCollection
from noisi.borrowed_functions.voronoi_surface_area import calculate_surface_area_of_a_spherical_Voronoi_polygon
from noisi.borrowed_functions.voronoi_polygons import xyzToSpherical
import warnings
warnings.filterwarnings("ignore")


def get_voronoi_surface_area(grd):
    """
    Computes the spherical voronoi cells and calculates their surface areas.
    Input: grid with longitude and latitude 
    Output: grid (since the order might change) and voronoi surface areas corresponding to each point
    
    Functions from:
        https://github.com/tylerjereddy/spherical-SA-docker-demo/blob/master/docker_build/demonstration.py
        https://github.com/MITHaystack/scikit-discovery/blob/master/skdiscovery/visualization/spherical_voronoi.py#L40
    """
    # convert grid into panda dataframe
    gridpd = {'lat': grd[1], 'lon': grd[0]}
    grid_data = pd.DataFrame(data=gridpd)
    
    # Calculate the vertices for the voronoi cells
    voronoi = getVoronoiCollection(data=grid_data,lat_name='lat',lon_name='lon',full_sphere=True)
    
    # Calculate the surface area for each voronoi cell
    voronoi_lat = []
    voronoi_lon = []
    voronoi_area = []
    
    for i in range(0,np.size(voronoi.points,0)):
        P_cart = xyzToSpherical(x=voronoi.points[i,0],y=voronoi.points[i,1],z=voronoi.points[i,2])
        voronoi_lat.append(P_cart[0])
        voronoi_lon.append(P_cart[1])
        vert_points = voronoi.vertices[voronoi.regions[i]]
        area = calculate_surface_area_of_a_spherical_Voronoi_polygon(vert_points,6371)
        voronoi_area.append(area)
        if i%1000 == 0:
            print('%g of %g voronoi cell surface areas calculated.' %(i,np.size(voronoi.points,0)),flush=True)
        
    # Reassign grd so that everything is in the right order
    grd = np.asarray([voronoi_lon,voronoi_lat])
    voronoi_area = np.asarray(voronoi_area)
    print('All voronoi cell surface areas calculated.')
    
    return grd, voronoi_area
