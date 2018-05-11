import os
import numpy as np

def test_smoothing():
	
	os.mkdir('test/testdata/testsrc/step_0/grad')
	os.system('python util/smoothing.py test/testdata/testsrc/step_0/grad_archived/grad_all.npy \
		test/testdata/testsrc/step_0/grad/grad_smooth.npy test/testdata/sourcegrid.npy \
		10000.0 95 1e-16')

	# assert the results are the same
	# ToDo: path
	#n1 = NoiseSource('test/testdata/testsrc/step_1_archived/starting_model.h5')
	#n2 = NoiseSource('test/testdata/testsrc/step_1/starting_model.h5')

	grad_old = np.load('test/testdata/testsrc/step_0/grad_archived/grad_smooth.npy')
	grad = np.load('test/testdata/testsrc/step_0/grad/grad_smooth.npy')
	
	# Additional step: Rounding. The smoothed gradient is of the order 10^-16. If precision errors occur, the test fails. So, round first with a precision of 20 digits after 0. 
	grad = np.around(grad,decimals=20)
	grad_old = np.around(grad_old,decimals=20)
	
	assert (grad_old[0,:] == grad[0,:]).sum() == len(grad[0,:])

	# remove stuff
	os.system('rm -rf test/testdata/testsrc/step_0/grad/')
	
