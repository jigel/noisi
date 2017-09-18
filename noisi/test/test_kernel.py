import os
import numpy as np

def test_kernel():
	# copy the correlation
	os.system('cp test/testdata/testsrc/step_0/adjt_archived/*.sac \
		test/testdata/testsrc/step_0/adjt/')
	os.system('cp -r test/testdata/testsrc/wavefield_processed_archived\
		test/testdata/testsrc/wavefield_processed')

	# run forward model
	os.system('noisi kernel test/testdata/testsrc/ 0')

	# assert the results are the same
	# ToDo: path
	k1 = np.load('test/testdata/testsrc/step_0/kern/NET.STA1..CHA--NET.STA2..CHA.npy')
	k2 = np.load('test/testdata/testsrc/step_0/kern_archived/NET.STA1..CHA--NET.STA2..CHA.npy')
	assert (k1 == k2).sum() == len(k1)

	# remove stuff
	os.system('rm test/testdata/testsrc/step_0/adjt/*')
	os.system('rm -rf test/testdata/testsrc/step_0/kern/*')
	os.system('rm -rf test/testdata/testsrc/wavefield_processed')