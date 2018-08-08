# plotting on the map
import matplotlib as mpl
mpl.rcParams['font.size'] = 14
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import matplotlib.tri as tri
import numpy as np
import time


def plot_grid(map_x,map_y,map_z,stations=[],v=None,globe=False,
    outfile=None,title=None,shade='flat',cmap=None,
    sequential=False,v_min=None,normalize=False,
    coastres='c',proj='cyl',
    lat_0=None,lon_0=None,lonmin=None,lonmax=None,
    latmin=None,latmax=None,mode='interp',resol=1,alpha=1.0):

    if lat_0 is None:
        lat_0  = 0.5*(map_y.max()-map_y.min())
    if lon_0 is None:
        lon_0 = 0.5*(map_x.max()-map_x.min())

    if lonmin == None:
        lonmin = np.min(map_x)
    if lonmax == None:
        lonmax = np.max(map_x)
    if latmax == None:
        latmax = np.max(map_y)
    if latmin == None:
        latmin = np.min(map_y)


    if resol != 1:
        map_x = map_x[::resol]
        map_y = map_y[::resol]
        map_z = map_z[::resol]


    if not proj == 'ortho':
        m = Basemap(rsphere=6378137,resolution=coastres,
        projection=proj,lat_0=lat_0,lon_0=lon_0,
        llcrnrlat=latmin,urcrnrlat=latmax,
        llcrnrlon=lonmin,urcrnrlon=lonmax)
    else:
        m = Basemap(rsphere=6378137,resolution=coastres,
        projection=proj,lat_0=lat_0,lon_0=lon_0)

    plt.figure(figsize=(11,9))
    plt.subplot(111)

    if title is not None:
        plt.title(title)



    if normalize:
        map_z /= np.max(np.abs(map_z))


    if v is None:
        v = np.max(map_z)

    if sequential:
        cm = plt.cm.magma
        if v_min == None:
            v_min = 0.
    else:
        cm = plt.cm.bwr
        v_min =-v

    if cmap is not None:
        cm = cmap




    print('max. value on map: %g' %map_z.max())
    if mode == 'interp':
        pass
        # triangulate first, then project,
        # and use plt.tripcolor to put it on the map.
        triangles = tri.Triangulation(map_x,map_y)
        (triangles.x,triangles.y) = m(triangles.x,triangles.y)
        #if it doesn't work, use pcolor mode

        plt.tripcolor(triangles,map_z,shading=shade, vmin=v_min,
            vmax=v,cmap=cm,alpha=alpha,linewidth=0.0,edgecolor='none')

        cbar = m.colorbar(location='bottom',pad=0.3)

    elif mode == 'pcolor':
        mx, my = m(map_x,map_y)
        m.pcolor(mx,my,map_z,cmap=cm,
            tri=True,shading=shade,vmin=v_min,
            vmax=v,alpha=alpha,linewidth=0.,edgecolor='none')
        cbar = m.colorbar(location='bottom',pad=0.3)


    elif mode == 'srclocs':

        indx = map_x % 3 <= 0.5
        map_x = map_x[indx].astype(int)
        map_y = map_y[indx].astype(int)

        indy = map_y % 3 <= 0.5

        mx,my = m(map_x[indy],map_y[indy])
        m.scatter(mx,my,marker='o',c='0.5',lw=1.,s=np.ones(len(mx))*0.2)

    elif mode == 'srcdots':


        colors = cm(map_z)
        indx = abs(map_z) > 0.4*np.max(map_z)
        sizes = np.ones(len(map_x))*1

        mx,my = m(map_x,map_y)


        #m.scatter(mx[indx],my[indx],marker='o',c=colors[indx],s=sizes[indx])
        scplt = m.scatter(mx,my,marker='o',c=map_z,cmap=cm,s=sizes)
        cbar = m.colorbar(scplt,location='bottom',pad=0.3)

    if normalize and v==1.:
        cbar.set_ticks([-1.0,-0.5,0.,0.5,1.0])
    elif normalize and v!=1.:
        cbar.set_ticks([-v,-v/2.,0.,v/2.,v])


    if globe:
       m.drawcoastlines(linewidth=1.,color='0.5')
    else:
       m.drawcoastlines(linewidth=1.0)


    if globe:
       #pass
       m.drawparallels(np.arange(-90.,120.,30.),labels=[1,0,0,0],color='0.5') # draw parallels
       m.drawmeridians(np.arange(-180,210,60.),labels=[0,0,0,1],color='0.5') # draw meridians

    else:
        if not proj == 'ortho':
           d_lon = round(abs(lonmax-lonmin) / 3.)
           d_lat = round(abs(latmax-latmin) / 3.)
           parallels = np.arange(latmin,latmax,d_lat).astype(int)
           meridians = np.arange(lonmin,lonmax,d_lon).astype(int)
           m.drawparallels(parallels,labels=[1,0,0,0]) # draw parallels
           m.drawmeridians(meridians,labels=[0,0,0,1])

    #draw station locations
    for sta in stations:
        m.plot(sta[0],sta[1],'^',color='r',markersize=15,markeredgecolor='0.5',latlon=True)
        #m.plot(sta[0],sta[1],'^',color='lime',markersize=5,markeredgecolor='0.5',latlon=True)
    if outfile is None:
        plt.show()
    else:
        plt.savefig(outfile,dpi=300.)
        plt.close()

def plot_sourcegrid(gridpoints,**kwargs):

    plt.figure()
    plt.subplot(111)
    m = Basemap(rsphere=6378137,**kwargs)
    m.drawcoastlines()

    m.plot(gridpoints[0],gridpoints[1],'go',markersize=10.,latlon=True)
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

def plot_envelope(corr,env):
    max_lag = (corr.stats.npts-1)/2 * corr.stats.delta
    lag = np.linspace(-max_lag,max_lag,corr.stats.npts)
    plt.plot(lag,corr.data/np.max(np.abs(corr.data)),'k')
    plt.plot(lag,env/np.max(np.abs(env)),'r',linewidth=2.)
    plt.grid()
    plt.xlabel('Lag (s)')
    plt.ylabel('Normalized correlation / normalized envelope')
    plt.show()
