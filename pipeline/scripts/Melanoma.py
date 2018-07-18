#################################################################
#################################################################
############### Melanoma Analysis - Python Support ############
#################################################################
#################################################################
##### Author: Julia Zhao
##### Affiliation: Ma'ayan Laboratory,
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
import h5py
import numpy as np
import pandas as pd
import math

##### 2. Custom modules #####
# Pipeline running


#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####

#######################################################
#######################################################
########## S1. Zichen's supporting functions
#######################################################
#######################################################

#############################################
########## 1. Load read counts
#############################################

def load_read_counts(h5_file, gsms=[]):
	'''Load data from h5 file using a list of GSMs.
	'''
	f = h5py.File(h5_file, 'r')
	mat = f['data']['expression']

	gsms = [str.encode(x) for x in gsms]

	all_gsms = f['meta']['Sample_geo_accession']

	sample_mask = np.in1d(all_gsms, gsms)
	if sample_mask.sum() == 0:
		raise RuntimeError('None of the gsms exists in the h5 file')
	else:
		sample_ids = all_gsms[sample_mask]
		genes = f['meta']['genes']
		# to prevent MongoDB error
		genes = map(lambda x:x.decode('utf-8').replace('.', '_'), genes)		
		# Retrieve gene by sample matrixpytho
		expr_df = pd.DataFrame(mat[sample_mask, :].T, index=genes, columns=sample_ids)

		# Filter out non-expressed genes
		expr_df = expr_df.loc[expr_df.sum(axis=1) > 0, :] 

		# Fix column names
		expr_df.columns = [x.decode('utf-8') for x in expr_df.columns]

		return expr_df
		
#############################################
########## 2. Compute CPMs
############################################# 
def compute_CPMs(expr_df, CPM_cutoff=0.3, at_least_in_persent_samples=10):
	'''Convert expression counts to CPM.
	'''
	n_samples = expr_df.shape[1]
	at_least_in_n_samples = int(math.ceil(at_least_in_persent_samples/100. * n_samples))

	expr_df = (expr_df * 1e6) / expr_df.sum(axis=0)
	# Filter out lowly expressed genes
	mask_low_vals = (expr_df > CPM_cutoff).sum(axis=1) > at_least_in_n_samples
	expr_df = expr_df.loc[mask_low_vals, :]
	return expr_df
