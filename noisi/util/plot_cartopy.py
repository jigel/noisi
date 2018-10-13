## plotting using cartopy

import cartopy.crs as ccrs
#import cartopy
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from pandas import read_csv
import os

def plot_gradient(project_path, source_path,gridlines = False,extent=None):
    """
    This function uses cartopy to plot the gradient. 
    
    Input: path to source, extent = [lon_min,lon_max,lat_min,lat_max]
    
    """
    
    #source_path = source_path
    #project_path = project_path
    
    grad_file = np.load(os.path.join(source_path,'step_0/grad/grad_all.npy')).T
    grd = np.load(os.path.join(project_path,'sourcegrid.npy'))
    stationlist = read_csv(os.path.join(project_path,'stationlist.csv'))
    source_name = os.path.basename(os.getcwd())
    
    # get all stations for plotting
    stations_lat = []
    stations_lon = []
    
    for i in range(0,np.size(stationlist,0)):
        lat = stationlist.at[i,'lat']
        lon = stationlist.at[i,'lon']
        stations_lat.append(lat)
        stations_lon.append(lon)
    
    # normalise the data for plotting
    v=np.max(np.abs(grad_file))
    
    # plot
    plt.figure(figsize=(25,10))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    if gridlines:
        ax.gridlines(draw_labels=True)
    if extent is not None:
        ax.set_extent(extent)
    plt.scatter(grd[0],grd[1],s=30,c=grad_file[:,0],transform=ccrs.Geodetic(),cmap=plt.get_cmap('seismic'),norm=colors.Normalize(-v,v))
    plt.colorbar()
    plt.scatter(stations_lon,stations_lat,c='k',marker='^')
    plt.title('Gradient for {} with {} stations'.format(source_name,np.size(stations_lat)),y=1.05)
    
    