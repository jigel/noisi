import numpy as np
from math import sqrt, pi
import sys
from mpi4py import MPI
from warnings import warn
from noisi.util.plot import plot_grid
# Try yet another: sort of Gaussian convolution, but determining the distance
# in cartesian coordinates.


def get_distance(gridx,gridy,gridz,x,y,z):
    #def distance_function(x1,y1,z1,x2,y2,z2):
    #    return sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    #dist = np.array([distance_function()])
    xd = gridx - x
    yd = gridy - y
    zd = gridz - z
    
    return np.sqrt(np.power(xd,2)+np.power(yd,2)+np.power(zd,2))



def smooth_gaussian(values,coords,rank,size,sigma,r=6371000.,threshold=1e-9):

	# coords format: (lon,lat)

	# step 1: Cartesian coordinates of map
	theta = np.deg2rad(-coords[1] + 90.) 
	phi = np.deg2rad(coords[0] + 180.)

	x = r*np.sin(theta) * np.cos(phi)
	y = r*np.sin(theta) * np.sin(phi)
	z = r*np.cos(theta)


	v_smooth = np.zeros(values.shape)


	a = 1./(sigma*sqrt(2.*pi))
	
	for i in range(rank,len(values),size):
		
		xp,yp,zp = (x[i],y[i],z[i])
		dist = get_distance(x,y,z,xp,yp,zp)
		weight = a * np.exp(-(dist)**2/(2*sigma**2))
		#print(weight.max())
		# I just had an idea for 'sparsity' here; test this:

		idx = weight >= threshold

		if idx.sum() == 0:
			warn('No weights above threshold, reset threshold.')
			v_smooth[i] = 0.

		else:
			v_smooth[i] = np.sum(np.multiply(weight[idx],values[idx])) / idx.sum()
		

	return v_smooth


def apply_smoothing_sphere(rank,size,values,coords,sigma,cap=95,threshold=1.e-12):


	sigma = float(sigma)
	cap = float(cap)
	threshold = float(threshold)

	# clip
	perc_up = np.percentile(values,cap,overwrite_input=False)
	perc_dw = np.percentile(values,100-cap,overwrite_input=False)
	values = np.clip(values,perc_dw,perc_up)

	

	# get the smoothed map; could use other functions than Gaussian here
	v_s = smooth_gaussian(values,coords,rank,size,sigma,threshold=threshold)
	

	
	comm.barrier()
	
	# collect the values
	print('Gathering...')
	v_s_all = comm.gather(v_s,root=0)
	# rank 0: save the values
	if rank == 0:
		
		print('Gathered.')
		v_s = np.zeros(v_s.shape)
		for i in range(size):

			v_s += v_s_all[i]
		
		return(v_s)

def test_gauss_smoothing(sourcegrid,map):
	#
	grd = np.load(sourcegrid)[:,0:10000]
	v = np.ones(grd.shape[1])
	ihalf = grd.shape[1] // 2
	v[ihalf:] = 10
	np.save('temp_coord.npy',grd)
	np.save('temp_vals.npy',v)
	plot_grid(grd[0],grd[1],v)

	smooth_map = apply_smoothing_sphere('temp_vals.npy',
		'test','temp_coord.npy',500000)
	print(smooth_map.shape)

	plot_grid(grd[0],grd[1],smooth_map)


if __name__=='__main__':

	# initialize parallel comm
	comm = MPI.COMM_WORLD
	size = comm.Get_size()
	rank = comm.Get_rank()

	# pass in: input_file, output_file, coord_file, sigma
	# open the files
	inputfile = sys.argv[1]
	outputfile = sys.argv[2]
	coordfile = sys.argv[3]
	

	#sigma = float(sys.argv[4])
	sigma = sys.argv[4].split(',')
	for ixs in range(len(sigma)):
		sigma[ixs] = float(sigma[ixs])

	cap = float(sys.argv[5])
	
	try:
		thresh = float(sys.argv[6])
	except IndexError:
		thresh = 1.e-12

	coords = np.load(coordfile)
	values = np.array(np.load(inputfile),ndmin=2)
	smoothed_values = np.zeros(values.shape)
	

	for i in range(values.shape[0]):

		array_in = values[i,:]
		try:
			sig = sigma[i]
		except IndexError:
			sig = sigma[-1]

		v = apply_smoothing_sphere(rank,size,array_in,\
			coords,sig,cap,threshold=thresh)
		if rank == 0:
			smoothed_values[i,:] = v
			print(np.isnan(smoothed_values).sum())

	if rank == 0:
		np.save(outputfile,smoothed_values)

	
