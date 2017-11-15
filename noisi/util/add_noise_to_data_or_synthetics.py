import numpy as np
from glob import glob
from obspy import read
from pathlib import Path
import sys

if __name__ == '__main__':

	directory = sys.argv[1]

	perc_of_max = float(sys.argv[2])

	try:
		format = sys.argv[3]
	except IndexError:
		format = 'sac'

	try:
		o = sys.argv[4]
	except IndexError:
		o = input('Are you sure you want to add noise? [n]/yes:\n')



	traces = glob(directory +'/*.'+format.upper())
	traces += glob(directory + '/*.'+format.lower())
	
	if o != 'yes':
		sys.exit("Nothing added.")
	else:
		for t in traces:
			tr = read(t)[0]
			tr.data += np.random.randn(len(tr.data))*tr.data.max()*perc_of_max
			tr.write(t,format=format)

