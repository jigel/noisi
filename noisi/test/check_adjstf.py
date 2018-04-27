import numpy as np
import matplotlib.pyplot as plt
from obspy import Trace
from noisi.scripts import measurements as rm
from noisi.scripts import adjnt_functs as af


# more or less copying Korbi's test with my measurement and adjoint source 


# *********************************************
# input:
# *********************************************
#scale = 1e20 #ununsed
steps = np.arange(-14, 0, 0.1)
mtype = 'ln_energy_ratio'#'ln_energy_ratio'
sacdict = {'dist':1e6}
g_speed = 3700.
window_params                   =    {}
window_params['hw']             =    200
window_params['sep_noise']      =    1.
window_params['win_overlap']    =    False
window_params['wtype']          =    'hann'
window_params['causal_side']    =    False
window_params['plot']           =    False
# *********************************************
# *********************************************

# only for testing the test:
# def l2_simple(tr_1,tr_2):
# 	mf = np.sum(0.5 * (tr_1.data - tr_2.data) **2)
# 	adstf = (tr_1.data - tr_2.data)
# 	return mf,adstf

m_a_options = {'g_speed':g_speed,'window_params':window_params}
m_func = rm.get_measure_func(mtype)
a_func = af.get_adj_func(mtype)


# create observed data, synthetics and perturbation
c_obs = 2 * (np.random.rand(2401,) - 0.5)
c_ini = 2 * (np.random.rand(2401,) - 0.5)
d_c = 2 * (np.random.rand(2401,) - 0.5)
# form traces (measurement script works with obspy trace objects, not pure arrays)
c_obs = Trace(data=c_obs)
c_obs.stats.sampling_rate = 1.0
c_obs.stats.sac = sacdict

c_syn = Trace(data=c_ini)
c_syn.stats.sampling_rate = 1.0
c_syn.stats.sac = sacdict



# obtain a measurement and an adjoint source time function
# for the unperturbed measurement
msr_o = m_func(c_obs,**m_a_options)
msr_s = m_func(c_syn,**m_a_options)
data, success = a_func(c_obs,c_syn,**m_a_options)

if mtype == 'energy_diff':
	data = data[0]
	msr_s = msr_s[0]
	msr_o = msr_o[0]
data *= (msr_s-msr_o)

j = 0.5*(msr_s-msr_o)**2

# testing the test:
# j,data = l2_simple(c_syn,c_obs)

# left hand side of test 1: adjt source time function * du = change of misfit wrt u
djdc = np.dot(data,d_c) 

# right hand side of test 1: Finite difference approx of misfit change for different steps


dcheck = []
d_ch = c_syn.copy()


for step in steps:
	d_ch.data = c_ini + 10. ** step * d_c
	msr_sh = m_func(d_ch,**m_a_options)
	if mtype == 'energy_diff':	
		msr_sh = msr_sh[0]

	jh = 0.5 * (msr_sh - msr_o)**2
	# testing the test:
	# jh, dn = l2_simple(d_ch,c_obs)
	djdch = (jh - j) / (10.**step) 
	dcheck.append(abs(djdc - djdch) / abs(djdc))
	
# plot
plt.semilogy(steps,dcheck)
plt.title("Check for adjoint source time function")
plt.show()


