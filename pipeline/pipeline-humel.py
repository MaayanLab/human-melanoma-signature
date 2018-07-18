#################################################################
#################################################################
############### Melanoma Signature Profile ################
#################################################################
#################################################################
##### Author: Julia Zhao
##### Affiliation: Ma'ayan Laboratory,
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
from ruffus import *
import sys, os, json, operator 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import functools

##### 2. Custom modules #####
# Pipeline running (import supporting codes)
sys.path.append('pipeline/scripts')
from Melanoma import *
# import geode, RNAseq, sklearn.metrics.pairwise

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####
# State the paths of the infiles
archs_matrix = 'rawdata/human_matrix.h5'
metadata_file = 'rawdata/m2.csv'
melanoma_masterfile = 'rawdata/mela_rna_3rep_short.csv'

#######################################################
#######################################################
########## S1. Prepare melanoma expression data
#######################################################
#######################################################
# We process the ARCHS4 human gene expression data, to exact RNA sequencing data of the samples of interest

#############################################
########## 1. Read data from ARCHS4
#############################################

# Make a directory that will host all runned results
@mkdir('s1-expression_data.dir')
# Define infiles & outfile
# Merge two starting files as a 1 list of infiles, and export this single output to the result folder
@merge([archs_matrix, metadata_file],
	   's1-expression_data.dir/melanoma-counts.txt')

def getCounts(infiles, outfile):
	# Define infiles
	archs_matrix, metadata_file = infiles
	# Read infile metadata_file into pandas
	metadata_df = pd.read_csv(metadata_file)
	# Extract all gsms from the 'id' column
	allgsm = metadata_df['id'].values
	# Save only the unique values by converting the array to a set first, then to a list
	gsm_unique = list(set(allgsm))
	# Read expression data
	expression_dataframe = load_read_counts(archs_matrix, gsms = gsm_unique)
	# Set index name
	expression_dataframe.index.name = 'gene'
	# Save to outfile
	expression_dataframe.to_csv(outfile, sep="\t")

#############################################
########## 2. Normalize the gene expression table
#############################################
##### Use the method "count per million" to normalize gene expression
@transform(getCounts,
		   suffix('-counts.txt'),
		   '-cpm.txt')
def getCPM(infile, outfile):
	# Read the previously generated raw gene expression
	# index_col = False -> forces pandas to not use the 1st column as the index
	expression_dataframe = pd.read_table(infile, index_col='gene')
	# Perform counts per million normalization
	cpm_table = compute_CPMs(expression_dataframe)
	# Save the cpm normalized df to outfile
	cpm_table.to_csv(outfile, sep="\t")

#######################################################
#######################################################
########## S2. Generate signatures
#######################################################
#######################################################
# 

#############################################
########## 1. Label sample IDs with their status groups
#############################################
##### Make a dictionary where std_i is the key, and index of the df.
##### Create an array of give shape & type, filled with zeros, then replace with 1 or 2 depends on ctrl/perturb status.
##### Use the dictionary to create signatures of studies

@mkdir('s2-signatures.dir')

@merge([melanoma_masterfile, getCPM],
	   's2-signatures.dir/melanoma-idmatches.txt')

def getMatches(infiles, outfile):
	# Define infiles
	melanoma_masterfile, cpm = infiles
	# Read the previously generated infiles
	melanoma_masterdf = pd.read_csv(melanoma_masterfile,index_col='study_index')
	cpm_table = pd.read_table(cpm, index_col='gene')

	# create a dictionary called matches
	matches = {} 
	for std_i, row in melanoma_masterdf.iterrows(): # std_i is the key of the dictionary, index of the dataframe
		ctrl_gsms = row['ctrl_id'].split(' ')
		pert_gsms = row['perturb_id'].split(' ')
		match = np.zeros(cpm_table.shape[1], dtype=np.int32) # Return a new array of given shape and type, filled with zeros.
		match[np.in1d(cpm_table.columns, ctrl_gsms)] = 1
		match[np.in1d(cpm_table.columns, pert_gsms)] = 2
		matches[std_i] = match

	matches_df=pd.DataFrame(matches,index=cpm_table.columns)
	matches_df.index.name='sample'
	matches_df.to_csv(outfile, sep="\t")
#############################################
########## 2. Create signatures
#############################################
##### Describe the step
### Input: 
### Output: 
@transform(getMatches,
		   suffix('-idmatches.txt'),
		   add_inputs(getCPM),
		   '-signatures.txt')
def getSignatures(infiles,outfile):
	print (infiles, outfile)
	# # create an empty dataframe
	# cd_results = pd.DataFrame(index=cpm_table.index)
	# d_std_i_cd = {} # to top up/down genes

	# for std_i, match in matches.items():
	# 	cd_res = geode.chdir(cpm_table.values, match, cpm_table.index, 
	# 					gamma=.5, sort=False, calculate_sig=False)
	# 	cd_coefs = np.array(map(lambda x: x[0], cd_res))
	# 	cd_results[std_i] = cd_coefs
		
	# 	# sort CD in by absolute values in descending order
	# 	srt_idx = np.abs(cd_coefs).argsort()[::-1]
	# 	cd_coefs = cd_coefs[srt_idx][:600]
	# 	sorted_DEGs = cpm_table.index[srt_idx][:600]
		
	# 	# split up and down
	# 	up_genes = dict(zip(sorted_DEGs[cd_coefs > 0], cd_coefs[cd_coefs > 0]))
	# 	dn_genes = dict(zip(sorted_DEGs[cd_coefs < 0], cd_coefs[cd_coefs < 0]))
	# 	d_std_i_cd[std_i+'-up'] = up_genes
	# 	d_std_i_cd[std_i+'-dn'] = dn_genes

	# 	# Save cd_results as your signature table
	# 	cd_results.to_csv(outfile, sep='\t')
##################################################
##################################################
########## Run pipeline
##################################################
##################################################

#############################################
########## 2. Normalize the gene expression table
#############################################
##### Describe the step
### Input: 
### Output: 

##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1, forcedtorun_tasks=[])
print('Done!')