import os
import numpy as np

def test_gradient():
	# copy data
	os.system('cp test/testdata/testsrc/step_0/kern_archived/NET.STA1..CHA--NET.STA2..CHA.npy \
		test/testdata/testsrc/step_0/kern/NET.STA1..CHA--NET.STA2..CHA.0.npy')
	os.system('cp test/testdata/testsrc/step_0/ln_energy_ratio.measurement_archived.csv\
		test/testdata/testsrc/step_0/ln_energy_ratio.0.measurement.csv')
	os.system('cp test/testdata/testsrc/step_0/starting_model_archived.h5\
		test/testdata/testsrc/step_0/starting_model.h5')

	# run forward model
	os.system('noisi gradient test/testdata/testsrc/ 0')

	# assert the results are the same
	# ToDo: path
	
	g1 = np.load('test/testdata/testsrc/step_0/grad_archived/grad_all.npy')
	g2 = np.load('test/testdata/testsrc/step_0/grad/grad_all.npy')[0,:]

	assert (g1 == g2).sum() == len(g1)
	
	
	# remove stuff
	os.system('rm test/testdata/testsrc/step_0/grad/*')
	os.system('rm test/testdata/testsrc/step_0/kern/*')
	os.system('rm test/testdata/testsrc/step_0/starting_model.h5')
	os.system('rm test/testdata/testsrc/step_0/ln_energy_ratio.0.measurement.csv')