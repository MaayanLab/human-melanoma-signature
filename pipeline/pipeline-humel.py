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
# Differentially expressed genes are compiled as a signature/vector per study

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

#############################################
########## 4. Extract each study's signature individually
#############################################
##### Prepare each study's differentially expressed gene vector
@subdivide(getSignatures,
		   formatter(),
		   's2-signatures.dir/melanoma-signature-study-*.txt',
		   's2-signatures.dir/melanoma-signature-study-')
def getSinglesig(infile, outfiles, outfileRoot):	
	sig_df = pd.read_table(infile).set_index('gene')
	# Loop through studies 
	for study in sig_df.columns[:]:
		print('Doing '+ study)
		# Make each signature (column) to a np.array
		study_signature = sig_df[study].values
		outfile = '{outfileRoot}{study}.txt'.format(**locals())
		np.savetxt(outfile, study_signature)

#############################################
########## 5. Find cluster intersection/union
#############################################
##### Find intersections of genes between all 6 clusters		
@transform(getSignatures,
		   suffix('-signatures.txt'),
		   '-top50sig-union.txt')
def getUnion(infile, outfile):
	sig_df = pd.read_table(infile).set_index('gene')
	# Make each signature to a np.array
	s17a=sig_df['17_a'].values
	s17b=sig_df['17_b'].values
	s17c=sig_df['17_c'].values
	s2a=sig_df['2_a'].values
	s2b=sig_df['2_b'].values
	s2c=sig_df['2_c'].values
	s15a=sig_df['15_a'].values
	s15b=sig_df['15_b'].values
	s15c=sig_df['15_c'].values
	s18a=sig_df['18_a'].values
	s18b=sig_df['18_b'].values
	s18c=sig_df['18_c'].values
	s19a=sig_df['19_a'].values
	s19b=sig_df['19_b'].values
	s19c=sig_df['19_c'].values
	s19d=sig_df['19_d'].values
	s1=sig_df['1'].values
	s5a=sig_df['5_a'].values
	s5b=sig_df['5_b'].values
	s5c=sig_df['5_c'].values
	s11=sig_df['11'].values
	s3=sig_df['3'].values
	# Take vector sum per cluster as a representation of intersection
	cluster1_inter= s17a + s17b + s17c + s2a + s2b + s2c
	cluster2_inter= s15c + s19b + s19d
	cluster3_inter= s15a + s15b + s3 + s18b + s11 + s19a
	cluster4_inter= s18c + s5b +s5c
	cluster5_inter= s1 + s19c
	cluster6_inter= s18a + s5a
	# Take absolute values to include both up & down, sort the array in reserve order
	sort_idx_1 = np.argsort(np.abs(cluster1_inter))[::-1]
	sort_idx_2 = np.argsort(np.abs(cluster2_inter))[::-1]
	sort_idx_3 = np.argsort(np.abs(cluster3_inter))[::-1]
	sort_idx_4 = np.argsort(np.abs(cluster4_inter))[::-1]
	sort_idx_5 = np.argsort(np.abs(cluster5_inter))[::-1]
	sort_idx_6 = np.argsort(np.abs(cluster6_inter))[::-1]
	# Return an index of top 50 summed genes for this cluster
	sorted_genes_1 = sig_df.index[sort_idx_1][:50]
	sorted_genes_2 = sig_df.index[sort_idx_2][:50]
	sorted_genes_3 = sig_df.index[sort_idx_3][:50]
	sorted_genes_4 = sig_df.index[sort_idx_4][:50]
	sorted_genes_5 = sig_df.index[sort_idx_5][:50]
	sorted_genes_6 = sig_df.index[sort_idx_6][:50]
	# Union of all 6 clusters's top 50 differentially exp. genes
	export = set(sorted_genes_1) | set(sorted_genes_2) | set(sorted_genes_3) | set(sorted_genes_4) | set(sorted_genes_5) | set(sorted_genes_6)	export.to_csv(outfile, sep='\t', index=False)
	# Subset the original signature matrix with the union genes we extracted
	export_df=sig_df.loc[export]
	export_df.to_csv(outfile, sep='\t')
	# Intersection
	# center= set(sorted_genes_1) & set(sorted_genes_2) & set(sorted_genes_3) & set(sorted_genes_4) & set(sorted_genes_5) & set(sorted_genes_6)
	# top500inter=sig_df.loc[center]
	# top500inter.to_csv(outfile, sep='\t')

#############################################
########## 6. Find other general overlap
#############################################
##### Find top overlap of genes between studies, besides union and intersection
##### Use rank score to find top & bottom 250 genes across studies
@transform(getSignatures,
		   suffix('.txt'),
		   '-80percent.txt')
def getPercent(infile,outfile):
	sig_df=pd.read_table(infile, index_col='gene')
	sig_df_ranked=sig_df.rank()# Rank matrix
	like=np.zeros_like(sig_df_ranked) # Create an all zero array
	# define top 250 and bottom ranked 250 genes as 1
	like[sig_df_ranked<250]=1
	like[sig_df_ranked>(18992-250)]= 1
	ls=like.sum(axis=1) # Row sums of all 18992 row entries
	# subset 17/22=0.8, genes that are present top 250, 80% at the time
	subset80=sig_df[ls>17]
	subset80.to_csv(outfile, sep='\t')
	
######################################################
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
	enr_res = []
	# Iterate through dictionary, create gene lists and save as 'genes'
	for signature in dict_top_genes:
		print('Doing '+signature)
		for direction in ('top', 'bottom'):
			genes = dict_top_genes[signature][direction] # Value = dict [1st key] [2nd key]
			result_dataframe = run_enrichr(genes)
			# Create new columns with dictionary value from direction
			result_dataframe['direction'] = direction 
			result_dataframe['signature'] = signature
			enr_res.append(result_dataframe)
	# Concatenate & save
	results = pd.concat(enr_res)
	results.to_csv(outfile, sep='\t', index=False)
	# In case you want to write a dictionary instead
	# enr_res = defaultdict(dict)
	# # Iterate through dictionary, create gene lists and save as 'genes'
	# for signature in dict_top_genes:
	# 	print('Doing '+signature)
	# 	for direction in ('top', 'bottom'):
	# 		genes = dict_top_genes[signature][direction] # value = dict [1st key] [2nd key]
	# 		enr_res[signature][direction] = run_enrichr(genes)
	# # freeze default dict for read only
	# enr_res.default_factory = None
#############################################
########## 2. Analyze enrichr results cross-studies
#############################################
##### Find out the significantly overlapping upregulated pathways via FDR values
@subdivide(getEnrichr,
		   formatter(),
		   's4-enrichment.dir/enrichr-results-*.txt',
		   's4-enrichment.dir/enrichr-results-')

def getTopfdr (infile, outfiles, outfileRoot):
	enr_res_df = pd.read_table(infile)
	# Select rows whose column value equals a scalar 'top'

	for direction in ['top', 'bottom']:
		all_top = enr_res_df.loc[enr_res_df['direction'] == direction]
		
		# Create a table of pathways to study names, filtered by FDR value
		# Fill na value is set to 1 since that'd be the lowest for log transformation later
		toptable= all_top.pivot(index='term_name', columns = 'signature',values='FDR').fillna(1)
		# Log 10 transformation
		all_top_log = -np.log10(toptable)
		outfile = '{outfileRoot}{direction}.txt'.format(**locals())
		all_top_log.to_csv (outfile, sep='\t')

	# # Convert 2D array to a boolean array
	# mask = (all_top_log > 0.000001).sum(axis=1)
	# # Every 'true' is score 1, and use row sum to filter
	# fdrs_filtered = all_top_log[mask<10]
#############################################
########## 3. Visualize FDR of up-regulated pathways
#############################################
##### Use seaborn to visualize consensus pathways 
@transform(getTopfdr,
		   suffix('.txt'),
		   '-FDR-heatmap.png')
def getTopfdrheat (infile, outfile):
	all_top_log = pd.read_table(infile, index_col='term_name')
	# Unbiased variance over columns, take top 20 results
	top_terms = all_top_log.var(axis=1).sort_values(ascending=False).index[:20]
	n = 20
	plot_df=all_top_log.copy()
	# Visualize with the seaborn package
	g2 = sns.clustermap(plot_df.loc[top_terms])
	# Export graph
	g2.savefig(outfile)
#############################################
########## 4. Intersection of pathways 
#############################################
##### Closely examine the intersection of pathway between all 6 clusters
#############################################
########## ?. Use Clustergrammer to create heatmap of up & down genes
#############################################
# from IPython.display import HTML, display
# # to display hyperlink as <a> tag in output cells
# def display_link(url):
#     raw_html = '<a href="%s" target="_blank">%s</a>' % (url, url)
#     return display(HTML(raw_html))

# clustergrammer_url = 'http://amp.pharm.mssm.edu/clustergrammer/matrix_upload/'
# r = requests.post(clustergrammer_url, files={'file': open(infile2, 'rb')})
# link = r.text
# display_link(link)

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