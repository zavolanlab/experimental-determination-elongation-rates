#!/usr/bin/env python

import numpy as np
import pandas as pd
import os
import gc
import pysam
import json
import statistics
import subprocess
from scipy.optimize import curve_fit


#file paths
tsv_file="/scicore/home/zavolan/schlus0000/riboseq_pipeline/snakemake/prepare_annotation/results/homo_sapiens/transcript_id_gene_id_CDS.tsv"

time_files_dict = { 0 : ("all_data_timecourse3_HepG2_CHXHRT/transcripts.mapped.unique.sam","all_data_timecourse3_HepG2_CHXHRT/p_site_offsets/alignment_offset.json"), \
	80 : ("all_data_timecourse3_HepG2_HRT_1_3min/transcripts.mapped.unique.sam","all_data_timecourse3_HepG2_HRT_1_3min/p_site_offsets/alignment_offset.json"), \
	90 : ("all_data_timecourse3_HepG2_HRT_1_5min/transcripts.mapped.unique.sam","all_data_timecourse3_HepG2_HRT_1_5min/p_site_offsets/alignment_offset.json"), \
	100 : ("all_data_timecourse3_HepG2_HRT_1_6min/transcripts.mapped.unique.sam","all_data_timecourse3_HepG2_HRT_1_6min/p_site_offsets/alignment_offset.json"), \
	120 : ("all_data_timecourse3_HepG2_HRT_2min/transcripts.mapped.unique.sam","all_data_timecourse3_HepG2_HRT_2min/p_site_offsets/alignment_offset.json"), \
	150 : ("all_data_timecourse3_HepG2_HRT_2_5min/transcripts.mapped.unique.sam","all_data_timecourse3_HepG2_HRT_2_5min/p_site_offsets/alignment_offset.json"), \
	180 : ("all_data_timecourse3_HepG2_HRT_3min/transcripts.mapped.unique.sam","all_data_timecourse3_HepG2_HRT_3min/p_site_offsets/alignment_offset.json") }


#column names, choose from ['median_p_site','mean_p_site','SL_location']
column_names = ['SL_location']
error_column_names = ['err_SL_location']


#parameters
SL_binsize = 15
number_jackknife_bins=12
cutoff_number_reads_per_gene=100
cutoff_chisq_red=100.0


#define global dictionaries
CDS_start_ind = {}
CDS_stop_ind = {}
CDS_lengths = {}
transcript2gene = {}

CDS_mapped_p_sites = {}
CDS_mapped_reads = {}
TPM = {}



'''
method reads the .tsv file and initializes: 
 - a dictionary with start locations of the CDS
 - a dictionary with end locations of the CDS
 - a dictionary with CDS lengths
'''
def read_tsv(tsv_file):
	with open(tsv_file) as CDS_coordinates:
		for line in CDS_coordinates:
			sp_line = line.strip().split("\t")
			transcript_ID = sp_line[0]
			gene_ID = sp_line[1]
			CDS_start_ind[transcript_ID] = int(sp_line[2]) - 1
			CDS_stop_ind[transcript_ID] = int(sp_line[3])
			transcript2gene[transcript_ID] = gene_ID
			CDS_lengths[transcript_ID] = CDS_stop_ind[transcript_ID] - CDS_start_ind[transcript_ID]
	
	transcript2gene['all'] = 'all'



'''
method reads bam file
with possible file omissions from omit_start to omit_end
and saves result into a dictionary with p_site_locations on each transcript 
using the respective offset json file
need to call read_tsv first to set CDS coordinates
omits reads on UTRs and on the first 3 nts
returns dictionary
'''
def bam2dict(sam_filename,json_filename,jackknife_bin_ind,min_CDS_length=0):
	sam = pysam.AlignmentFile(sam_filename, 'rb')
	with open(json_filename, 'r') as json_file:
		 json_data = json_file.read()
	offset_obj = json.loads(json_data)
	
	ctr = 0
	all_p_sites = []
	CDS_mapped_p_sites = {}
	for read in sam.fetch():
		if read.is_reverse:
			continue
		
		transcript = read.reference_name
		read_length = len(read.seq)
		
		if ctr % number_jackknife_bins != jackknife_bin_ind and CDS_lengths[transcript] > min_CDS_length and str(read_length) in offset_obj:
			p_site_pos = (read.reference_start + offset_obj[str(read_length)] - CDS_start_ind[transcript])
			
			if p_site_pos in range(4,CDS_lengths[transcript]-2):
				if transcript not in CDS_mapped_p_sites:
					CDS_mapped_p_sites[transcript] = [p_site_pos]
				else:
					CDS_mapped_p_sites[transcript].append(p_site_pos)
					
				all_p_sites.append(p_site_pos)
		ctr += 1
	
	CDS_mapped_p_sites['all'] = all_p_sites
	
	return CDS_mapped_p_sites



'''
method computes TPM and reads per gene dictionary
needs input dictionary of p_sites per transcript
returns TPM and reads per gene dicts
CDS_mapped_reads omits reads
below the read per gene cutoff
appends length dict by global maximum of length in 'all'
'''
def compute_TPM_and_reads(CDS_mapped_p_sites):
	TPM = {} 
	CDS_mapped_reads = {}
	
	max_len = 0
	for transcript in CDS_mapped_p_sites:
		if transcript != 'all':
			if CDS_lengths[transcript] > max_len:
				max_len = CDS_lengths[transcript]
	CDS_lengths['all'] = max_len
	
	for transcript in CDS_mapped_p_sites:
		tmp_num = len(CDS_mapped_p_sites[transcript])
		TPM[transcript] = tmp_num / CDS_lengths[transcript]
		if tmp_num > cutoff_number_reads_per_gene:
			CDS_mapped_reads[transcript] = tmp_num
	
	total = 0.0
	for transcript in TPM:
		if transcript != 'all':
			total += TPM[transcript]
	
	for transcript in TPM:
		if total > 0.0:
			TPM[transcript] *= 1e6 / total
		else:
			TPM[transcript] = np.nan
	
	return TPM, CDS_mapped_reads



'''
method computes the mean ribosome location per transcript
needs input dictionary of p_sites per transcript
returns mean dict
'''
def mean_ribo_location(CDS_mapped_p_sites):
	mean_ribo_locs = {}
	
	for transcript in CDS_mapped_p_sites:
		mean_ribo_locs[transcript] = statistics.mean(CDS_mapped_p_sites[transcript])
	
	return mean_ribo_locs
	



'''
method computes the median ribosome location per transcript
needs input dictionary of p_sites per transcript
returns median dict
'''
def median_ribo_location(CDS_mapped_p_sites):
	median_ribo_locs = {}
	
	for transcript in CDS_mapped_p_sites:
		median_ribo_locs[transcript] = statistics.median(CDS_mapped_p_sites[transcript])
	
	return median_ribo_locs



'''
find maximum of the difference
'''
def find_max_difference(CDS_mapped_p_sites):
	max_difference = {}

	for transcript in CDS_mapped_p_sites:
		tlength = CDS_lengths[transcript]

		number_bins = tlength // SL_binsize + 1	#first codon is empty anyway

		data = np.array(CDS_mapped_p_sites[transcript][SL_binsize:])

		binned_data, _ = np.histogram(data, bins=number_bins)

		diff =  binned_data[1:] - binned_data[:-1]

		max_difference[transcript] = (np.argmax(diff) + 1) * SL_binsize

	return max_difference



'''
method computes the ribosome profiles
needs input dictionary of p_sites per transcript
returns profiles DataFrame
'''
def ribo_profile(CDS_mapped_p_sites,norm_win_lower,norm_win_upper):
	profiles = {}
	
	for transcript in CDS_mapped_p_sites:
		tmp_profile = np.zeros(norm_win_upper // SL_binsize)
		for p_site_pos in CDS_mapped_p_sites[transcript]:
			if p_site_pos in range(0,norm_win_upper):
				ind = p_site_pos // SL_binsize
				tmp_profile[ind] += 1.0
		profiles[transcript] = tmp_profile
	
	return profiles



'''
method computes and returns overall dataframe with raw data
specify omitted data
and minimum CDS length
'''
def compute_overall_df(jackknife_bin_index,min_CDS_length=3000,norm_win_lower=2400,norm_win_upper=3000):
	read_tsv(tsv_file)
	
	tmp_df_trans2gene = pd.DataFrame.from_dict(transcript2gene,orient='index',columns=['gene_ID'])
	
	overall_df = pd.DataFrame()
	for key, value in time_files_dict.items():
		CDS_mapped_p_sites = bam2dict(value[0],value[1],jackknife_bin_index,min_CDS_length)

		TPM, CDS_mapped_reads = compute_TPM_and_reads(CDS_mapped_p_sites)
		mean_p_sites = mean_ribo_location(CDS_mapped_p_sites)
		median_p_sites = median_ribo_location(CDS_mapped_p_sites)
		max_difference = find_max_difference(CDS_mapped_p_sites)
		p_site_profiles = ribo_profile(CDS_mapped_p_sites,norm_win_lower,norm_win_upper)
		
		tmp_df_jackknife_bin_index = pd.DataFrame([[transcript,jackknife_bin_index] for transcript in CDS_mapped_p_sites],columns=['transcript_ID','jackknife_bin_index'])
		tmp_df_jackknife_bin_index.set_index('transcript_ID',inplace=True)
		
		tmp_df_timediff = pd.DataFrame([[transcript,key] for transcript in CDS_mapped_p_sites],columns=['transcript_ID','timediff'])
		tmp_df_timediff.set_index('transcript_ID',inplace=True)
		
		#here because not all of them must be populated in all jackknife bins
		tmp_df_CDS_lengths = pd.DataFrame.from_dict(CDS_lengths,orient='index',columns=['CDS_length'])
		
		tmp_df_TPM = pd.DataFrame.from_dict(TPM,orient='index',columns=['TPM'])
		tmp_df_CDS_mapped_reads = pd.DataFrame.from_dict(CDS_mapped_reads,orient='index',columns=['reads'])
		
		tmp_df_median = pd.DataFrame.from_dict(median_p_sites,orient='index',columns=['median_p_site'])
		tmp_df_mean = pd.DataFrame.from_dict(mean_p_sites,orient='index',columns=['mean_p_site'])
		tmp_df_maxdiff = pd.DataFrame.from_dict(max_difference,orient='index',columns=['max_difference'])
		tmp_df_profiles = pd.DataFrame.from_dict(p_site_profiles,orient='index',columns=['profile_bin_'+str(i) for i in range(norm_win_upper // SL_binsize)])

		tmp_overall_df = pd.concat([tmp_df_trans2gene,tmp_df_CDS_lengths,tmp_df_jackknife_bin_index,tmp_df_timediff,tmp_df_TPM,tmp_df_CDS_mapped_reads,tmp_df_median,tmp_df_mean,tmp_df_maxdiff,tmp_df_profiles],axis=1,join='inner')
		tmp_overall_df.index.name = 'transcript_ID'
		tmp_overall_df.reset_index()
		
		overall_df = pd.concat([overall_df,tmp_overall_df],axis=0,join='outer')
	return overall_df



'''
function normalizes profile
to a given window (in nt, needs to be divided by SL_binsize)
returns normalized profile
'''
def normalize_profile(unnorm_profile,norm_win_lower,norm_win_upper):
	norm_win_length = norm_win_upper - norm_win_lower
	
	if len(unnorm_profile) < norm_win_upper // SL_binsize:
		return []
	
	norm = 0.0
	for index in range(norm_win_lower // SL_binsize,norm_win_upper // SL_binsize):
		norm += unnorm_profile[index]
	
	if norm > 0.0:
		normalized_profile = [ x / norm * norm_win_length / SL_binsize for x in unnorm_profile ]
		return normalized_profile
	else:
		return []



'''
method determines SL location (see Ingolia et. al.: Mammalian translation speed (2011))
which is the point where the ribosome saturationin the HRTringtonine sample 
#reaches predefined point (default 0.5) of the reference sample
needs reference and HRTringtonine sample as input as well as the saturation level
'''
def determine_saturation_point(ref_profile,profile,saturation_level=0.5,omitted_nt=120):
	for index in range(omitted_nt // SL_binsize,min(len(ref_profile),len(profile))):
		if profile[index-1] < saturation_level*ref_profile[index-1] and profile[index] >= saturation_level*ref_profile[index]:
			return index*SL_binsize
	
	return np.NAN



'''
convert profile columns to SL points
returns new dataframe
'''
def convert_profiles_to_SL(df,norm_win_lower=2400,norm_win_upper=3000):
	ret_df = df[['transcript_ID','gene_ID','jackknife_bin_index','timediff','CDS_length','TPM','reads',*column_names[:-1]]].copy()
	
	ret_df['SL_location'] = np.nan
	
	for transcript in df.transcript_ID.unique():
		for jackknife_bin_index in range(0,number_jackknife_bins+1):
			#reference profile
			unnorm_ref_profile = df.loc[(df['transcript_ID']==transcript) & (df['timediff']==0),df.columns[6:]].values.reshape(-1).tolist()
			
			norm_ref_profile = normalize_profile(unnorm_ref_profile,norm_win_lower,norm_win_upper)
			
			if len(norm_ref_profile) >= (norm_win_upper // SL_binsize):
				#iterate through the HRT samples
				for timediff in df.loc[(df['transcript_ID']==transcript) & (df['timediff'] > 0)].timediff.unique():
					unnorm_profile = df.loc[(df['transcript_ID']==transcript) & (df['timediff']==timediff),df.columns[6:]].values.reshape(-1).tolist()
			
					norm_profile = normalize_profile(unnorm_profile,norm_win_lower,norm_win_upper)
					
					if len(norm_profile) >= (norm_win_upper // SL_binsize):
						tmp = determine_saturation_point(norm_ref_profile,norm_profile)
						if np.isnan(tmp) == 0:
							ret_df.loc[(ret_df['transcript_ID']==transcript) & (ret_df['timediff']==timediff) & (ret_df['jackknife_bin_index']==jackknife_bin_index),['SL_location']] = tmp 
						
	return ret_df



'''
converts dataframe with jacknife bins
to one with best estimates and errors
for a list of columns
returns new dataframe
CAUTION: an overall timing error of 
9nt is added to every statistical error
'''
def convert_jackknife_to_errors(df,colnames):
	#fill return dataframe with best estimates
	ret_df = df[df['jackknife_bin_index']==number_jackknife_bins].filter(['transcript_ID','gene_ID','timediff','CDS_length','TPM','reads',*colnames]).copy()
	
	for colname in colnames:
		error_colname = "err_"+colname
		
		#add error column
		ret_df[error_colname] = np.nan
		
		#fill error columns
		for transcript in ret_df.transcript_ID.unique():
			for timediff in time_files_dict:
				tmp_list = df.loc[(df['transcript_ID']==transcript) & (df['timediff']==timediff),colname].to_list()
				
				if len(tmp_list) > 1:
					best_estimate = tmp_list[-1]
					tmp_list.pop()
				else:
					best_estimate = np.nan
					
				tmp = 0.0
				ctr = 0
				for jbr in tmp_list:
					if np.isnan(jbr)==0:
						tmp += (jbr - best_estimate) ** 2
						ctr += 1
				
				#modify back!
				if ctr > 1:
					tmp *= (ctr - 1) / ctr
					tmp = np.sqrt(tmp) + 6.0
				else:
					tmp = 6.0
				
				ret_df.loc[(df['transcript_ID']==transcript) & (df['timediff']==timediff),error_colname] = tmp
		
	return ret_df



'''
linear fit function
at point x with slope m*m and intercept b
CAUTION: m>=0 enforced!
'''
def linear(x, m, b):
	return m*m*x + b



'''
determines chi^2 reduced
at:
	- fit function callable fit_func
	- fit parameter vector p
	- location vector x
	- function value vector y
	- error vector dy
'''
def chisqred(fit_func, p, x, y, dy):
	ret = 0.0
	
	for xi, yi, dyi in zip(x,y,dy):
		ret += ( (fit_func(xi, *p) - yi) / dyi ) ** 2
	
	dof = len(x) - len(p)
	
	if dof > 0:
		ret /= dof
	else:
		ret = np.nan
	
	return ret



'''
function takes input array for fitting:
x-value array
y-value array
dy-value array
and removes nan values in y
and adds pseudocount of 0.1 to errors in dy
returns updated versions
'''
def prepare_fitting_data(x,y,dy):
	rx = []
	ry = []
	rdy = []
	for cx,cy,cdy in zip(x,y,dy):
		if np.isnan(cy)==0:
			rx.append(cx)
			ry.append(cy)
			
			if np.isnan(cdy)==0 and cdy > 0.0:
				rdy.append(cdy+0.1)
			else:
				rdy.append(0.1)
		
	return rx,ry,rdy



'''
converts dataframe with timepoints and errors
to one with fitted estimates of translation speed
(automatically adapts fitting range)
for a list of columns and their respective error columns
automatically adapts fitting range if ribosomes would have fallen off
returns new dataframe
'''
def fit_elong_speeds(df,colnames,error_colnames):
	#fill return dataframe with best estimates
	min_time = min(time_files_dict,key=time_files_dict.get)
	ret_df = df[df['timediff']==min_time].filter(['transcript_ID','gene_ID','CDS_length','TPM','reads']).copy()
	
	for colname, error_colname in zip(colnames,error_colnames):
		fit_colname = "elong_speed_"+colname
		fit_error_colname = "error_"+fit_colname
		fit_chisqred_colname = "chisq_red_"+fit_colname
		
		#add fit and error column
		ret_df[fit_colname] = np.nan
		ret_df[fit_error_colname] = np.nan
		ret_df[fit_chisqred_colname] = np.nan
		
		#fill error columns
		for transcript in ret_df.transcript_ID.unique():
			trnscrpt_CDS_len = ret_df.loc[ret_df['transcript_ID']==transcript,'CDS_length'].values[0]
			
			tx = df.loc[df['transcript_ID']==transcript,'timediff'].to_list()
			ty = df.loc[df['transcript_ID']==transcript,colname].to_list()
			tdy = df.loc[df['transcript_ID']==transcript,error_colname].to_list()
			x, y, dy = prepare_fitting_data(tx,ty,tdy)
			
			if len(x) > 2:
				#possibly automatically adapt fitting range
				popt, pcov = curve_fit(linear, x, y, sigma=dy)
				perr = np.sqrt(np.diag(pcov))
				tmp_chisqred = chisqred(linear, popt, x, y, dy)

				if (popt[0]**2 * x[-1] < trnscrpt_CDS_len and tmp_chisqred < cutoff_chisq_red) or transcript == 'all':
					ret_df.loc[df['transcript_ID']==transcript,fit_colname] = popt[0] ** 2
					ret_df.loc[df['transcript_ID']==transcript,fit_error_colname] = abs(2 * perr[0] * popt[0])
					ret_df.loc[df['transcript_ID']==transcript,fit_chisqred_colname] = tmp_chisqred
				
	return ret_df



'''
BEGIN SCRIPT
'''
if os.path.isfile("p_sites.tsv") == False:
	all_df = pd.DataFrame()
	for jackknife_ind in range(number_jackknife_bins+1):
		tmp_jackknife_df = compute_overall_df(jackknife_ind,min_CDS_length=3000)
		all_df = pd.concat([all_df,tmp_jackknife_df],axis=0,join='outer')
		print('Bin number '+str(jackknife_ind)+' processed.')


	#unset index
	all_df.reset_index(inplace=True)
	all_df.rename(columns={'index' : 'transcript_ID'},inplace=True)

	#save dataframe to csv file
	all_df.to_csv("p_sites.tsv",sep="\t",index=False)
	print("P-site positions calculated.")
else:
	all_df = pd.read_csv("p_sites.tsv",sep="\t")
	print("P-site positions read.")


if os.path.isfile("locations.tsv") == False:
	if 'SL_location' in column_names:
		loc_df = convert_profiles_to_SL(all_df)
		loc_df.to_csv("locations.tsv",sep="\t",index=False)
	else:
		loc_df = all_df.copy()
	print("SL-location calculated.")
else:
	loc_df = pd.read_csv("locations.tsv",sep="\t")
	print("SL-location read.")

del all_df
gc.collect()


if os.path.isfile("locations_with_errors.tsv") == False:
	df_w_err = convert_jackknife_to_errors(loc_df,column_names)
	df_w_err.to_csv("locations_with_errors.tsv",sep="\t",index=False)
	print("Jackknife errors estimated.")
else:
	df_w_err = pd.read_csv("locations_with_errors.tsv",sep="\t")
	print("Jackknife errors read.")

del loc_df
gc.collect()


fit_df = fit_elong_speeds(df_w_err,column_names,error_column_names)
fit_df.dropna(subset=["elong_speed_" + column_name for column_name in column_names],thresh=1,inplace=True)
fit_df.to_csv("elongation_speeds.tsv",sep="\t",index=False)


