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
	
	assert ((grad - grad_old)/grad_old*100.).max() < 1.e-10 

	# remove stuff
	os.system('rm -rf test/testdata/testsrc/step_0/grad/')
	
