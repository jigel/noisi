import numpy as np
import pandas as pd
import os
import json

from glob import glob
from math import isnan
import sys
from noisi.util.corr_pairs import glob_obs_corr
from noisi.my_classes.noisesource import NoiseSource
from warnings import warn
####################################
# ToDo: more fancy and more secure with click or argparse
source_model = sys.argv[1]
oldstep = int(sys.argv[2])
grad_file = sys.argv[3]
grad_old = sys.argv[4]
update_mode = sys.argv[5]#'conjgrad'# steepest, conjgrad
min_snr = float(sys.argv[6])#min_snr = 5.0
min_stck = int(sys.argv[7])#min_stck = 320.
nr_msr = int(sys.argv[8])#nr_msr = 300
step_length = float(sys.argv[9])#step_length = 
mode = sys.argv[10] # 'max', 'random'
# Give as part per hundred, e.g 0.1 for 10%
perc_step_length = None
# include those data points in the test which are at or above this
# fraction of maximum misfit:
perc_of_max_misfit = 0.6666
# Only if the following is set to True, a small subset (nr_msr) 
# of data will be selected and copied and their misfit evaluated
# for a step length test. Otherwise, only the update of the source model
# is performed. 

####################################


def _update_steepestdesc(model,
	neg_grad,
	step_length=None,
	perc_step_length=None,
	project=False,
	smooth=False):

	if step_length is not None and perc_step_length is not None:
		raise ValueError('Only one of step_length and perc_step_length can\
			be specified.')

# just in case:
	#os.system('cp {} {}'.format(model,model+'.bak'))

	

# smooth the model
	if smooth:
		raise NotImplementedError('Sorry, not implemented yet.')

# project the model
	if project:
		raise NotImplementedError('Sorry, not implemented yet.')
		# if implemented, here should be a projection of the new kernel
		# onto the distr_basis functions, thereby yielding the new distr_weights

	else:
	# the spectrum remains unchanged.
	# assuming there is one basis only, this will be updated with the new kernel
		if step_length is not None:	
			#src_model.model['distr_basis'][:] += neg_grad * step_length
			descent_direction = neg_grad * step_length
		elif perc_step_length is not None:
			#src_model.model['distr_basis'][:] += neg_grad/np.max(np.abs(neg_grad)) * perc_step_length
			descent_direction = neg_grad/np.max(np.abs(neg_grad)) * perc_step_length

# write to file
# close the underlying h5 file	
	#src_model.model.close()

	return(descent_direction)

def _update_conjugategrad(
	model,
	neg_grad,
	old_grad,
	old_upd,
	updatename,
	step_length
	):
# just in case:
	#os.system('cp {} {}'.format(model,model+'.bak'))

	#src_model = NoiseSource(model)

	# determine beta
	norm_neggrad = np.linalg.norm(-1.*neg_grad,ord=2,axis=(1))
	norm_oldgrad = np.linalg.norm(old_grad,ord=2,axis=(1))

	beta = np.power(norm_neggrad,2) / np.power(norm_oldgrad,2)
	print(beta)
	print(neg_grad[0:10])
	beta_update = np.dot(beta,old_upd)
	print(beta_update[0:10])

	upd = neg_grad + beta_update
	np.save(updatename,upd)

	# save the update so it can be used to determine the next step
	#np.save(updatename,upd)

	#src_model.model['distr_basis'][:] += step_length * upd
	descent_direction = step_length * upd
	

	print(step_length*upd[0:20])
	
	
	return(descent_direction,upd)



def _prepare_test_steplength(msrfile,source_config,newdir):
	
	obs_dir = os.path.join(source_config['source_path'],'observed_correlations')


	# Read in the csv files of measurement.
	for mfile in msrfile:
		if 'data' not in locals():
			data = pd.read_csv(mfile)
		else:
			# We get the addition of both datasets, which means that l2_norms of all
			# measurements are added up and the stations pairs with max overall misfit are chosen
			data.l2_norm += pd.read_csv(mfile).l2_norm
			data.nstack += pd.read_csv(mfile).nstack


	# Get a set of n randomly chosen station pairs. Criteria: minimum SNR, 
	# ---> prelim_stations.txt
	
	# this makes no sense for combined measure
	#data_accept = data[(data.snr >= min_snr)]
	#if len(data_accept) == 0:
#		raise ValueError('No data match selection criteria.')

	#data_accept = data_accept[(data_accept.snr_a >= min_snr)]
	#if len(data_accept) == 0:
#		raise ValueError('No data match selection criteria.')

	#data_accept = data_accept[(data_accept.nstack >= min_stck)]
	data_accept = data[(data.nstack >= min_stck)]
	if len(data_accept) == 0:
		raise ValueError('No data match selection criteria.')
	
	data_accept = data_accept[~(data_accept.l2_norm.apply(np.isnan))]
	if len(data_accept) == 0:
		raise ValueError('No data match selection criteria.')

	data_accept = data_accept[~(data_accept.snr.apply(np.isnan))]
	if len(data_accept) == 0:
		raise ValueError('No data match selection criteria.')

	data_accept = data_accept[~(data_accept.snr_a.apply(np.isnan))]
	if len(data_accept) == 0:
		raise ValueError('No data match selection criteria.')
	

	# select data...
	if mode =='random':
		data_select = data_accept.sample(n=nr_msr)
	elif mode == 'max':
		data_select1 = data_accept.sort_values(by='l2_norm',na_position='first')
		data_select = data_select1.iloc[-nr_msr:]
	
	print(data_select)

	#data_select = pd.concat([data_select1,data_select2])
	
	#stafile = open(os.path.join(newdir,'stations_slt.txt'),'w')
	#stafile.write("Station pairs to be used for step lenght test:\n")

	inffile = open(os.path.join(newdir,'step_length_test_info.txt'),'w')
	inffile.write('Parameters:\n')
	inffile.write('source_model: %s\n' %source_model)
	inffile.write('old step: %s\n' %oldstep)
	inffile.write('min_snr %g\n' %min_snr)
	inffile.write('min_stck %g\n' %min_stck)

	if step_length is not None:
		inffile.write('step_length %g\n' %step_length)
	elif perc_step_length is not None:
		inffile.write('step_length as fraction of max. weight%g\n' %perc_step_length)
	inffile.write('-'*40)
	inffile.write("\nStation pairs to be used for step lenght test:\n")

	cum_misf = 0.0
	# Take care of the test set for the step length test
	
	for i in data_select.index:
		
		sta1 = data_select.at[i,'sta1'].split('.')[0:4]
		sta2 = data_select.at[i,'sta2'].split('.')[0:4]
		
		#lat1 = data_select.at[i,'lat1']
		#lat2 = data_select.at[i,'lat2']
		#lon1 = data_select.at[i,'lon1']
		#lon2 = data_select.at[i,'lon2']

		misf = data_select.at[i,'l2_norm']

		cum_misf += misf
		# synthetics in the old directory?
		#synth_filename = os.path.join(datadir,'corr','{}--{}.sac'.format(sta1,sta2))
		#print(synth_filename)
		# copy the relevant observed correlation, oh my
		obs_dir = os.path.join(source_config['source_path'],'observed_correlations')
		obs_correlations = glob_obs_corr('{}.{}.{}.{}'.format(*sta1),
			'{}.{}.{}.{}'.format(*sta2),obs_dir,ignore_network=True)
		

		if len(obs_correlations) > 0:

			# Use preferentially '', '00' channels.
			obs_correlations.sort()
			corr = obs_correlations[0]

			#sta1 = sta1.split('.')
			#sta2 = sta2.split('.')
			#stafile.write('{} {} {} {}\n'.format(*(sta1[0:2]+[lat1]+[lon1])))
			#stafile.write('{} {} {} {}\n'.format(*(sta2[0:2]+[lat2]+[lon2])))

			#inffile.write('{} {}, {} {} L2 misfit: {}\n'.format(*(sta1[0:2]+sta2[0:2]+[misf])))
			inffile.write('{}, {} L2 misfit: {}\n'.format(sta1,sta2,misf))

			os.system('cp {} {}'.format(corr,os.path.join(newdir,'obs_slt')))
		


	inffile.write('-'*40)
	inffile.write('\nCumulative misfit: %g\n' %cum_misf)
	inffile.write('-'*40)
	inffile.close()
	#stafile.close()

	return()

############ Preparation procedure #################################################
#prepare_test_steplength = False
# where is the measurement database located?
source_model = os.path.join(source_model,'source_config.json')
source_config=json.load(open(source_model))
datadir = os.path.join(source_config['source_path'],'step_' + str(oldstep))
msrfile = os.path.join(datadir,"{}.*.measurement.csv".format(source_config['mtype']))
msrfile = glob(msrfile)

# Initialize the new step directory
newstep = int(oldstep) + 1
newdir = os.path.join(source_config['source_path'],'step_' + str(newstep))

if not os.path.exists(newdir):
	newdir = os.path.join(source_config['source_path'],'step_' + str(newstep))
	os.mkdir(newdir)
	os.mkdir(os.path.join(newdir,'obs_slt'))
	os.mkdir(os.path.join(newdir,'corr'))
	os.mkdir(os.path.join(newdir,'adjt'))
	os.mkdir(os.path.join(newdir,'grad'))
	os.mkdir(os.path.join(newdir,'kern'))
	_prepare_test_steplength(msrfile,source_config,newdir)

os.system('cp {} {}'.format(os.path.join(datadir,'base_model.h5'),newdir))
os.system('cp {} {}'.format(os.path.join(datadir,'starting_model.h5'),newdir))


# Set up a prelim_sourcemodel.h5: 
# Contains starting model + step length * (-grad) for steepest descent
# This would be the point to project to some lovely basis functions..
grad = grad_file

neg_grad = -1. * np.load(grad)
old_grad = np.load(grad_old)

neg_grad = np.array(neg_grad,ndmin=2)
old_grad = np.array(old_grad,ndmin=2)

new_sourcemodel = os.path.join(newdir,'starting_model.h5')
new_update = os.path.join(newdir,'grad','update.npy')
old_upd = os.path.join(datadir,'grad','update.npy')

src_model = NoiseSource(new_sourcemodel,w='r+')

if not os.path.exists(old_upd):
	old_upd = -1. * old_grad.copy()
else:
	old_upd = np.load(old_upd)

if update_mode == 'steepest':

	descent_direction = _update_steepestdesc(new_sourcemodel,neg_grad,step_length=step_length,
	perc_step_length=perc_step_length,project=False,smooth=False)

elif update_mode == 'conjgrad':
	descent_direction, update = _update_conjugategrad(new_sourcemodel,neg_grad,old_grad,
	old_upd,new_update,step_length)
	np.save(new_update,update)


src_model.model['distr_basis'][:] += descent_direction


if src_model.model['distr_basis'][:].min() < 0.:
		warn('Step length leads to negative source model...reset values to be >=0.')
		src_model.model['distr_basis'][:] = src_model.model['distr_basis'][:].clip(0.0)

src_model.model.close()
# (outside of this script) forward model selected correlations
# (outside of this script) evaluate misfit for selected correlations
