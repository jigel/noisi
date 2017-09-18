import os
from obspy import read

def test_adjoint():
	# copy the correlation
	os.system('cp test/testdata/testsrc/step_0/corr_archived/*.sac \
		test/testdata/testsrc/step_0/corr/')

	# run forward model
	os.system('noisi measurement test/testdata/testsrc/ 0')

	# assert the results are the same
	# ToDo: path
	
	tr1 = read('test/testdata/testsrc/step_0/adjt/NET.STA1..CHA--NET.STA2..CHA.sac')[0]
	tr2 = read('test/testdata/testsrc/step_0/adjt_archived/NET.STA1..CHA--NET.STA2..CHA.sac')[0]
	
	assert (tr1.data == tr2.data).sum() == len(tr2.data)
	assert tr1.stats.sampling_rate == tr2.stats.sampling_rate
	
	
	# remove stuff
	os.system('rm test/testdata/testsrc/step_0/corr/*')
	os.system('rm -rf test/testdata/testsrc/step_0/adjt/*')