import os
from obspy import read

def test_forward():

	srcdir = os.path.join('test','testdata','testsrc')
	os.mkdir('test/testdata/testsrc/step_0/corr/')
	# copy the preprocessed wavefields
	os.system('mkdir '+os.path.join(srcdir,'wavefield_processed'))
	os.system('cp test/testdata/testsrc/step_0/starting_model_archived.h5\
		test/testdata/testsrc/step_0/starting_model.h5')
	os.system('cp test/testdata/testsrc/wavefield_processed_archived/*.h5 \
		test/testdata/testsrc/wavefield_processed')
        
	
	# run forward model
	os.system('noisi correlation %s 0' %srcdir)

	# assert the results are the same
	# ToDo: path
	tr1 = read('test/testdata/testsrc/step_0/corr/NET.STA1..CHA--NET.STA2..CHA.sac')[0]
	tr2 = read('test/testdata/testsrc/step_0/corr_archived/NET.STA1..CHA--NET.STA2..CHA.sac')[0]
	
	assert ((tr1.data - tr2.data)/tr1.data*100.).max() < 1.e-6 
	assert tr1.stats.sampling_rate == tr2.stats.sampling_rate

	# remove the resulting data and the preprocessed wavefields
	os.system('rm -rf test/testdata/testsrc/wavefield_processed/')
	os.system('rm -rf test/testdata/testsrc/step_0/corr/')
	os.system('rm test/testdata/testsrc/step_0/starting_model.h5')
