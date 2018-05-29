import numpy as np
import matplotlib.pyplot as plt
from meshpy.tet import build, Options, MeshInfo
import pyvtk
from pyvtk import PointData, Scalars
import os
import sys
import h5py

infile = sys.argv[1]
try: 
    sourcegrid = sys.argv[2]
except:
    pass
outfilename = os.path.splitext(infile)[0]

#============================================================
#- Read input.
#============================================================

print('Read data ...')
try:
    m = np.load(infile)
    grd = np.load(sourcegrid)
except:
    n = h5py.File(infile,'r')
    m = n['distr_basis'][:]
    grd = n['coordinates'][:]

lat = grd[1]
lon = grd[0]


#============================================================
#- Triangulation.
#============================================================

for i in range(0,m.shape[0]):
    S = m[i,:]
    outfile = outfilename+'.'+str(i)+'.vtk'
    print('Compute Delauney triangulation ...')
    
    x=6371.0*np.cos(lat*np.pi/180.0)*np.cos(lon*np.pi/180.0)
    y=6371.0*np.cos(lat*np.pi/180.0)*np.sin(lon*np.pi/180.0)
    z=6371.0*np.sin(lat*np.pi/180.0)
    
    pts=np.array((x, y, z)).T
    mesh_info=MeshInfo()
    mesh_info.set_points(pts)
    opts=Options("Q")
    mesh=build(mesh_info, options=opts)
    elements=mesh.elements
    
    #============================================================
    #- Write vtk file.
    #============================================================
    
    print('Write vtk file ...')
       
    vtkElements = pyvtk.VtkData(pyvtk.UnstructuredGrid(pts, tetra=elements), 
        PointData(Scalars(S, 'grad_PSD_ZZ')), "Mesh")
    vtkElements.tofile(outfile)
