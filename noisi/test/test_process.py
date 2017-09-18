import h5py
import os

def test_process():
	# run preprocessing
	os.system('noisi preprocess_synthetics test/testdata/testsrc/')
	# compare
	f1 = h5py.File('test/testdata/testsrc/wavefield_processed/NET.STA1..CHA.h5')
	f2 = h5py.File('test/testdata/testsrc/wavefield_processed_archived/NET.STA1..CHA.h5')

	assert (f1['data'][:][0]-f2['data'][:][0]).sum() < 1.e-16
	
	f1 = h5py.File('test/testdata/testsrc/wavefield_processed/NET.STA2..CHA.h5')
	f2 = h5py.File('test/testdata/testsrc/wavefield_processed_archived/NET.STA2..CHA.h5')

	assert f1['stats'].attrs['reference_station'] == f2['stats'].attrs['reference_station']

	# delete preprocessed directory
	os.system('rm -rf test/testdata/testsrc/wavefield_processed/')
