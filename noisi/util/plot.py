# plotting on the map
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import matplotlib.tri as tri    
import numpy as np
    
def plot_grid(map_x,map_y,map_z,stations=[],v=1.0,globe=False,outfile=None,title=None,shade='flat',
    sequential=False,v_min=None,normalize=True,coastres='c',proj='cyl'):
    
    lat_0  = 0.5*(map_y.max()-map_y.min())
    m = Basemap(rsphere=6378137,resolution=coastres,projection=proj,lat_0=0.,           lon_0=0.,llcrnrlat=np.min(map_y),urcrnrlat=np.max(map_y),
    llcrnrlon=np.min(map_x),urcrnrlon=np.max(map_x))
    if globe:
        map_x = np.append(map_x,map_x[0])
        map_y = np.append(map_y,map_y[0])
        map_z = np.append(map_z,map_z[0])
    triangles = tri.Triangulation(map_x,map_y)
    # tripcolor plot.
    plt.figure()
    plt.subplot(111)
    plt.gca().set_aspect('equal')
    if title is not None:
        plt.title(title)

    if sequential:
        cm = plt.cm.magma
        if v_min == None:
            v_min = 0.
    else:
        cm = plt.cm.bwr
        v_min =-v

    if normalize:
        map_z /= np.max(np.abs(map_z))
   
    plt.tripcolor(triangles,map_z,shading=shade, vmin=v_min,vmax=v,cmap=cm)
    m.colorbar(location='bottom',pad=0.4)
    m.drawcoastlines(linewidth=0.5)
    #m.drawparallels(np.arange(-90.,120.,30.),labels=[1,0,0,0]) # draw parallels
    #m.drawmeridians(np.arange(-180,210,60.),labels=[0,0,0,1]) # draw meridians
    d_lon = abs(map_x.max()-map_x.min()) / 5.
    d_lat = abs(map_y.max()-map_y.min()) / 5.

    parallels = np.arange(np.min(map_y),np.max(map_y),d_lat).astype(int)
    meridians = np.arange(np.min(map_x),np.max(map_x),d_lon).astype(int)
    m.drawparallels(parallels,labels=[1,0,0,0]) # draw parallels
    m.drawmeridians(meridians,labels=[0,0,0,1])

    #draw station locations
    for sta in stations:
        m.plot(sta[0],sta[1],'rv',latlon=True)
    
    if outfile is None:
        plt.show()
    else:
        plt.savefig(outfile,format='png')
    
def plot_sourcegrid(gridpoints,**kwargs):

    plt.figure()
    plt.subplot(111)
    m = Basemap(rsphere=6378137,**kwargs)
    m.drawcoastlines()
    
    m.plot(gridpoints[0],gridpoints[1],'+',markersize=10.,latlon=True)
    plt.show()
    

def plot_window(correlation, window, measurement):
    
    
    maxlag = correlation.stats.npts * correlation.stats.delta
    lag = np.linspace(-maxlag,maxlag,correlation.stats.npts)
    
    plt.plot(lag,correlation.data/np.max(np.abs(correlation.data)))
    plt.plot(lag,window/np.max(np.abs(window)),'--')
    plt.title(correlation.id)
    plt.text(0,-0.75,'Measurement value: %g' %measurement)
    plt.xlabel('Correlation Lag in seconds.')
    plt.ylabel('Normalized correlation and window.')
    
    plt.show()