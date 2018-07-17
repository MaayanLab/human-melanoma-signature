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
import sys, os, json, operator #sklearn.metrics.pairwise
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import functools

##### 2. Custom modules #####
# Pipeline running (supporting codes)
sys.path.append('pipeline/scripts')
from Melanoma import *
# import geode, RNAseq
#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####
archs_matrix = 'rawdata/human_matrix.h5'
metadata_file = 'rawdata/mela_rna_3rep_short.csv'

#######################################################
#######################################################
########## S1. Find signatures of each melanoma study
#######################################################
#######################################################
# We process the ARCHS4 human gene expression data, to exact RNA sequencing data of the samples of interest

#############################################
########## 1. Read data from ARCHS4
#############################################

# Make a directory that will host all runned results
@mkdir('s1-expression_data.dir')
# Merge two starting files as a 1 list, and export this single output to the result folder
@merge([archs_matrix, metadata_file],
	   's1-expression_data.dir/melanoma-counts.txt')


def getCounts(infiles, outfile):
		# Read infiles, create a list of gsms to be analyzed
	archs_matrix, metadata_file = infiles
	# read metadata file in to pandas
	metadata_df = pd.read_csv(metadata_file)
	# Save all ctrl_ids to a new dataframe, as lists of lists
	ctrlgsm = np.array(metadata_df['ctrl_id'].str.strip().str.split())
	# Save all perturb_ids to a new dataframe, as lists of lists
	perturbgsm = np.array(metadata_df['perturb_id'].str.strip().str.split())
	gsmlist=np.concatenate((ctrlgsm, perturbgsm))
	# Reduce lists of lists to 1 long list of strings
	gsmlist = functools.reduce (operator.concat, gsmlist)
	# Read expression data
	expression_dataframe = load_read_counts(archs_matrix, gsms = gsmlist)
	# Save to outfile
	expression_dataframe.to_csv(outfile, sep='\t')
#
#############################################
########## 2. Normalize the gene expression table
#############################################
##### Use the method "count per million" to normalize gene expression
# @transform(expression_dataframe,suffix('-counts.txt'),'_cpm.csv')

# def getCPM (infile, ourfile):
# 	cpm_table = compute_CPMs (infile)
# 	cpm_table.to_csv(outfile, sep="\t")

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