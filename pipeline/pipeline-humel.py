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
import sys, os, json, operator, requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sklearn.metrics.pairwise as smp
from scipy.spatial.distance import cosine
import seaborn as sns; sns.set(color_codes=True)
from collections import Counter
from collections import defaultdict

##### 2. Custom modules #####
# Pipeline running (import supporting codes)
sys.path.append('pipeline/scripts')
from Melanoma import *
import geode_jupies
from RNAseq import *

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
	# Split & define infiles
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
	# Split & define infiles
	melanoma_masterfile, cpm = infiles
	# Read the previously generated infiles
	melanoma_masterdf = pd.read_csv(melanoma_masterfile,index_col='study_index')
	cpm_table = pd.read_table(cpm, index_col='gene')
	# Create a dictionary called matches
	matches = {} 
	for std_i, row in melanoma_masterdf.iterrows(): # std_i is the key of the dictionary, index of the dataframe
		ctrl_gsms = row['ctrl_id'].split(' ')
		pert_gsms = row['perturb_id'].split(' ')
		match = np.zeros(cpm_table.shape[1], dtype=np.int32) # Return a new array of given shape and type, filled with zeros.
		match[np.in1d(cpm_table.columns, ctrl_gsms)] = 1
		match[np.in1d(cpm_table.columns, pert_gsms)] = 2
		matches[std_i] = match
	# Convert a dict to a df, assign GSM IDs as index named "sample"
	matches_df=pd.DataFrame(matches,index=cpm_table.columns)
	matches_df.index.name='sample'
	matches_df.to_csv(outfile, sep="\t")

#############################################
########## 2. Create signatures
#############################################
##### Use the "matches" map to generate signatures based on study_index, ctrl/pert status, and genes 
@transform(getMatches,
		   suffix('-idmatches.txt'),
		   add_inputs(getCPM),
		   '-signatures.txt')

def getSignatures(infiles,outfile):
	# Split & define infiles
	matches, cpm = infiles
	# Read infiles into pandas
	matches_df=pd.read_table(matches, index_col='sample')
	cpm_table=pd.read_table(cpm, index_col='gene')
	# Initialize results
	results = []
	# Looooop
	for std_i, match in matches_df.items():
		print('Doing {}...'.format(std_i))
		cd_res = geode_jupies.chdir(cpm_table.values, match.values, cpm_table.index, 
						gamma=.5, sort=False, calculate_sig=False)		
		cd_dataframe = pd.DataFrame(cd_res).rename(columns={1: 'gene', 0: 'CD'})
		cd_dataframe['signature'] = std_i
		results.append(cd_dataframe)
	# Concatenate and merge
	result_dataframe = pd.concat(results)
	result_table = result_dataframe.pivot(index='gene', columns='signature', values='CD')
	# Write
	result_table.to_csv(outfile, sep='\t')

#############################################
########## 3. Subset up & down regulated genes
#############################################
##### Create signatures of up/down genes separately
@transform(getSignatures,
		   suffix('-signatures.txt'),
		   '-geneset-top.json')

def getGenesets(infile, outfile):
	 # Create empty dictionary
	dict_top_genes ={}
	sig_df = pd.read_table(infile).set_index('gene')
	# Number of studies to be processed
	n = len(sig_df.columns)
	i = 0
	# Loop through studies 
	for study in sig_df.columns[:]:
		i += 1 # Print computing status
		print('Doing study {study} ({i}/{n})...'.format(**locals()))
		# Extract column
		col = sig_df[study].sort_values(ascending=False) # Sort values in decreasing order
		genesets = {
			"top":col.index[:300].tolist(),
			"bottom":col.index[-300:].tolist()
			}
		# Extract the top/bottom 300 genes from the index
		dict_top_genes[study] = genesets	
	# open file
	f = open(outfile, 'w')
	# write to file
	f.write(json.dumps(dict_top_genes, ensure_ascii=False, indent=4))
	# close file
	f.close()

#######################################################
#######################################################
########## S3. Explore relationships between different studies
#######################################################
#######################################################
# Visualize 'similarities' of different studies via cosine distance heatmap

#############################################
########## 1. Calculate cosine distances of every study pair
#############################################
@mkdir('s3-distances.dir')
@files(getSignatures,'s3-distances.dir/melanoma-distance.txt')

def getCosdis(infile, outfile):
	# Import infile table, transpose df so index is n_samples, which is required for pair-wise distance func
	signature_df = pd.read_table(infile, index_col='gene').T
	# Create an empty dataframe with study_index as both index and column labels
	cos_distance_df = pd.DataFrame(index=signature_df.index, columns=signature_df.index)
	# Compute the distance matrix from a vector array X
	cos_distance = smp.pairwise_distances(signature_df, metric = 'cosine')
	# Note that cosine distance is defined as 1.0 minus the cosine similarity, hence a completely different function
	# Fill in the empty df earlier with values from cos_distance array
	cos_distance_df[:]=cos_distance
	# Export table out
	cos_distance_df.to_csv(outfile, sep='\t')

#############################################
########## 2. Visualization of the distances
#############################################
##### Create heatmap of cos distances between studies
@transform(getCosdis,
		   suffix('-distance.txt'),
		   '-3repstudies-distance-heatmap.png')

def getDismap (infile,outfile):
	# Import distance file, set 1st column as index
	cos_distance_df = pd.read_table(infile, index_col=0)
	# Visualize with the seaborn package
	g1 = sns.clustermap(cos_distance_df)
	# Export graph
	g1.savefig(outfile)

#######################################################
#######################################################
########## S4. Enrichment Analysis 
#######################################################
#######################################################
# Explore pathways activation via Enrichr
# Use Enrichr to compare up&down gene sets between studies
#############################################
########## 1. Gene set enrichment analysis
#############################################
##### Create enrichr results from up & down gene lists
@mkdir('s4-enrichment.dir')
@files(getGenesets,'s4-enrichment.dir/melanoma-enrichr-results.txt')

def getEnrichr (infile, outfile):
	# Read infile
	with open(infile) as openfile:
		dict_top_genes = json.load(openfile)
	# Defaultdict defaults a value if that key has not been set yet
	enr_res = []
	# Iterate through dictionary, create gene lists and save as 'genes'
	for signature in dict_top_genes:
		print('Doing '+signature)
		for direction in ('top', 'bottom'):
			genes = dict_top_genes[signature][direction] # value = dict [1st key] [2nd key]
			result_dataframe = run_enrichr(genes)
			result_dataframe['direction'] = direction
			result_dataframe['signature'] = signature
			enr_res.append(result_dataframe)

	# Concatenate
	results = pd.concat(enr_res)
	results.to_csv(outfile, sep='\t', index=False)
	# freeze default dict for read only
	# enr_res.default_factory = None
	# open file
	# f = open(outfile, 'w')
	# write to file
	# f.write(json.dumps(enr_res, ensure_ascii=False, indent=4))
	# close file
	# f.close()

#############################################
########## 2. Count
#############################################
# c = Counter()
# for i in [enr_res[signature]['top'] for signature in enr_res.keys()]:
#     c.update(i['term_name'].iloc[:5])
#############################################
########## 3. Use Clustergrammer to create heatmap of up & down genes
#############################################


# # POST the expression matrix to Clustergrammer and get the URL
# def getClustergrammer (infile, outfile):
# 	clustergrammer_url = 'http://amp.pharm.mssm.edu/clustergrammer/matrix_upload/'
# 	r = requests.post(clustergrammer_url, files={'file': open(dict_top_genes, 'rb')})
# 	link = r.text
# 	display_link(link)

#############################################
########## 2. Step Name
#############################################
##### Describe the step
### Input: 
### Output: 
#############################################
########## 2. Step Name
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