import os
from noisi import NoiseSource

def test_update():
	# copy data
	os.mkdir('test/testdata/testsrc/step_0/corr')
	os.mkdir('test/testdata/testsrc/step_0/grad')
	
	os.system('cp test/testdata/testsrc/step_0/corr_archived/NET.STA1..CHA--NET.STA2..CHA.sac \
		test/testdata/testsrc/step_0/corr/NET.STA1..CHA--NET.STA2..CHA.sac')

	os.system('cp test/testdata/testsrc/step_0/grad_archived/grad_all.npy\
		test/testdata/testsrc/step_0/grad/grad_all.npy')
	os.system('cp test/testdata/testsrc/step_0/starting_model_archived.h5\
		test/testdata/testsrc/step_0/starting_model.h5')
	os.system('cp test/testdata/testsrc/step_0/ln_energy_ratio.measurement_archived.csv\
		test/testdata/testsrc/step_0/ln_energy_ratio.0.measurement.csv')
	

	# run forward model
	os.system('./test/testdata/testsrc/update.sh')

	# assert the results are the same
	# ToDo: path
	n1 = NoiseSource('test/testdata/testsrc/step_1_archived/starting_model.h5')
	n2 = NoiseSource('test/testdata/testsrc/step_1/starting_model.h5')

	assert (n1.distr_basis == n2.distr_basis).sum() == len(n1.distr_basis[0,:])
	
	
	# remove stuff
	os.system('rm -rf test/testdata/testsrc/step_0/grad')
	os.system('rm -rf test/testdata/testsrc/step_0/corr')
	os.system('rm -rf test/testdata/testsrc/step_1')
	os.system('rm test/testdata/testsrc/step_0/starting_model.h5')
	os.system('rm test/testdata/testsrc/step_0/ln_energy_ratio.0.measurement.csv')
