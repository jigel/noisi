import numpy as np
from glob import glob
from obspy import read
from pathlib import Path
import sys

if __name__ == '__main__':

	directory = sys.argv[1]
	try:
		format = sys.argv[2]
	except IndexError:
		format = 'sac'

	traces = glob(directory +'/*.'+format.upper())
	traces += glob(directory + '/*.'+format.lower())
	
	o = input('Are you sure you want to add noise? [n]/yes:\n')
	
	if o != 'yes':
		sys.exit("You changed your mind.")
	else:
		for t in traces:
			tr = read(t)[0]
			tr.data += np.random.randn(len(tr.data))*tr.data.max()*0.3
			tr.write(t,format=format)

