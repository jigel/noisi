import sys
from pandas import read_csv
import json
import os
import numpy as np
from glob import glob

# evaluate measurement
source_dir = sys.argv[1]
step = sys.argv[2]

# get weights for freq. bands
measr_config = json.load(open(os.path.join(source_dir,
	'measr_config.json')))
weights = measr_config['weights']
mtype = measr_config['mtype']

# for each freq. band,...
step_tests = glob(os.path.join(source_dir,'step_'+step,'steptest_*'))

result = np.zeros((len(step_tests),2))

for i in range(len(step_tests)):
	for j in range(len(weights)):
		# what step length?
		result[i,0]=(float(step_tests[i].split('_')[-1]))

		# load the measurements
		filename = '{}.{}.measurement.csv'.format(mtype,j)
		msr_file = os.path.join(step_tests[i],filename)
		dat = read_csv(msr_file)

	
		# take the average
		result[i,1] += dat.l2_norm.mean()/len(weights) * weights[j]


np.save('result_step_length_test.npy',result)

