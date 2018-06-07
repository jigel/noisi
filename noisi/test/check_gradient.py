import numpy as np
import matplotlib.pyplot as plt
from obspy import Trace, read
from noisi.scripts import measurements as rm
from noisi.scripts import adjnt_functs as af
from noisi import WaveField, NoiseSource
import os
import h5py

# necessary data must be located in the test/testdata directory.

# *********************************************
# input:
# *********************************************
steps = np.arange(-8, 0.5, 0.3)
mtype = 'ln_energy_ratio' # only ln_energy_ratio can be used.
g_speed = 3300.
window_params                   =    {}
window_params['bandpass'] 		=    None
window_params['hw']             =    20
window_params['sep_noise']      =    0.
window_params['win_overlap']    =    False
window_params['wtype']          =    'hann'
window_params['causal_side']    =    True
window_params['plot']           =    False
# *********************************************
# *********************************************

# only for testing the test:
# def l2_simple(tr_1,tr_2):
# 	mf = np.sum(0.5 * (tr_1.data - tr_2.data) **2)
# 	adstf = (tr_1.data - tr_2.data)
# 	return mf,adstf

# preparations:
os.mkdir('test/testdata/testsrc/step_0/corr')
os.system('cp -R test/testdata/testsrc/wavefield_processed_archived \
test/testdata/testsrc/wavefield_processed')
os.system('cp test/testdata/config_archived.json \
test/testdata/config.json')
os.system('cp test/testdata/testsrc/measr_config_archived.json \
	test/testdata/testsrc/measr_config.json')
os.system('cp test/testdata/testsrc/source_config_archived.json \
test/testdata/testsrc/source_config.json')


m_a_options = {'g_speed':g_speed,'window_params':window_params}
m_func = rm.get_measure_func(mtype) 

wf = WaveField('test/testdata/wavefield_vel/NET.STA1..CHA.h5')
nlocs = wf.stats['ntraces']


# create perturbation
d_q = 2 * (np.random.rand(nlocs,) - 0.5)



# evaluate original misfit and load original gradient
m_a_options = {'g_speed':g_speed,'window_params':window_params}
m_func = rm.get_measure_func(mtype)

# open the files....
obs = read('test/testdata/testsrc/observed_correlations/*.sac')[0]
syn = read('test/testdata/testsrc/step_0/corr_archived/*.sac')[0]
syn.stats.sac = {}
syn.stats.sac['dist'] = obs.stats.sac.dist
msr_o = m_func(obs,**m_a_options)
msr_s = m_func(syn,**m_a_options)

# unperturbed misfit
j = 0.5*(msr_s-msr_o)**2
# unperturbed gradient
grad = np.load('test/testdata/testsrc/step_0/grad_archived/grad_all.npy')

# left hand side of test 3: gradient * dq = change of misfit wrt q
grad_dq = np.dot(grad,d_q)

dcheck = []
# loop:
for step in steps:
# add perturbation to archived model --> current model
	os.system('cp test/testdata/testsrc/step_0/starting_model_archived.h5 test/testdata/testsrc/step_0/starting_model.h5')

	n = h5py.File('test/testdata/testsrc/step_0/starting_model.h5')
	
	n['distr_basis'][:] += 10.**step * d_q
	
	n.flush()
	n.close()
# run correlation

	os.system('noisi correlation test/testdata/testsrc 0')

# evaluate misfit and add to list.
	syn = read('test/testdata/testsrc/step_0/corr/*.sac')[0]
	syn.stats.sac = {}
	syn.stats.sac['dist'] = obs.stats.sac.dist
	msr_sh = m_func(syn,**m_a_options)

	jh = 0.5 * (msr_sh - msr_o)**2
	djdqh = (jh - j) / (10.**step) 
	
	dcheck.append(abs(grad_dq - djdqh) / abs(grad_dq))

# remove the current synthetics
	os.system('rm test/testdata/testsrc/step_0/corr/*')
	
# plot results

# plot
plt.semilogy(steps,dcheck)
plt.title("Check for gradient")
plt.show()

# clean up...
os.system('rmdir test/testdata/testsrc/step_0/corr')
os.system('rm test/testdata/config.json')
os.system('rm test/testdata/testsrc/source_config.json')
os.system('rm test/testdata/testsrc/measr_config.json')

os.system('rm -rf test/testdata/testsrc/wavefield_processed')


