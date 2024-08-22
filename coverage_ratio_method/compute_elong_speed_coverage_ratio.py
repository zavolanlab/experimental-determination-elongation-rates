#!/usr/bin/env python

import numpy as np
import pandas as pd
import os
import gc
import pysam
import json
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



#parameters
number_jackknife_bins=20
cutoff_chisq_red=100000.0


#define global dictionaries
CDS_start_ind = {}
CDS_stop_ind = {}
CDS_lengths = {}
transcript2gene = {}



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
			CDS_start_ind[transcript_ID] = int(sp_line[2])
			CDS_stop_ind[transcript_ID] = int(sp_line[3]) + 1
			transcript2gene[transcript_ID] = gene_ID
			CDS_lengths[transcript_ID] = CDS_stop_ind[transcript_ID] - CDS_start_ind[transcript_ID]
	
	transcript2gene['all'] = 'all'
	CDS_lengths['all'] = max(CDS_lengths.values())



'''
method reads bam file
with possible file omissions from omit_start to omit_end
returns reads on CDS wo start and stop
normalized by CDS reads on stop dictionary
'''
def bam2coverage_ratio(sam_filename,json_filename,jackknife_bin_ind,min_CDS_length=0):
	sam = pysam.AlignmentFile(sam_filename, 'rb')
	with open(json_filename, 'r') as json_file:
		 json_data = json_file.read()
	offset_obj = json.loads(json_data)

	ctr = 0
	readcount_first_half = {}
	readcount_second_half = {}
	for read in sam.fetch():
		if read.is_reverse:
			continue

		transcript = read.reference_name
		read_length = len(read.seq)

		if (ctr % number_jackknife_bins) != jackknife_bin_ind and CDS_lengths[transcript] > min_CDS_length and str(read_length) in offset_obj:
			p_site_pos = (read.reference_start + offset_obj[str(read_length)] - CDS_start_ind[transcript])	# 0-based

			if p_site_pos in range(45,CDS_lengths[transcript]//2):
				if transcript not in readcount_first_half:
					readcount_first_half[transcript] = 1
				else:
					readcount_first_half[transcript] += 1
			elif p_site_pos in range(CDS_lengths[transcript]//2,CDS_lengths[transcript]-45):
				if transcript not in readcount_second_half:
					readcount_second_half[transcript] = 1.0
				else:
					readcount_second_half[transcript] += 1.0

		ctr += 1

	readcount_first_half['all'] = sum(readcount_first_half.values())

	norm_readcount_first_half = {}
	for t in readcount_second_half:
		if t in readcount_first_half:
			norm_readcount_first_half[t] = (CDS_lengths[transcript]//2 - 45) * readcount_first_half[t] / readcount_second_half[t]
		else:
			norm_readcount_first_half[t] = 0.0

	total_reads = {}
	for t in CDS_lengths:
		total_reads[t] = 0.0

		if t in readcount_first_half:
			total_reads[t] += readcount_first_half[t]

		if t in readcount_second_half:
			total_reads[t] += readcount_second_half[t]

	return norm_readcount_first_half, total_reads



'''
method computes and returns overall dataframe with raw data
specify omitted data
and minimum CDS length
'''
def compute_overall_df(jackknife_bin_index,min_CDS_length):
	read_tsv(tsv_file)
	
	tmp_df_trans2gene = pd.DataFrame.from_dict(transcript2gene,orient='index',columns=['gene_ID'])
	
	overall_df = pd.DataFrame()
	for key, value in time_files_dict.items():
		coverage_ratio, total_readcounts = bam2coverage_ratio(value[0],value[1],jackknife_bin_index,min_CDS_length)
		
		tmp_df_jackknife_bin_index = pd.DataFrame([[transcript,jackknife_bin_index] for transcript in coverage_ratio],columns=['transcript_ID','jackknife_bin_index'])
		tmp_df_jackknife_bin_index.set_index('transcript_ID',inplace=True)
		
		tmp_df_timediff = pd.DataFrame([[transcript,key] for transcript in coverage_ratio],columns=['transcript_ID','timediff'])
		tmp_df_timediff.set_index('transcript_ID',inplace=True)
		
		#here because not all of them must be populated in all jackknife bins
		tmp_df_CDS_lengths = pd.DataFrame.from_dict(CDS_lengths,orient='index',columns=['CDS_length'])
		
		tmp_df_coverage_ratio = pd.DataFrame.from_dict(coverage_ratio,orient='index',columns=['coverage_ratio'])

		tmp_df_middle_reads = pd.DataFrame.from_dict(total_readcounts,orient='index',columns=['middle_reads'])

		tmp_overall_df = pd.concat([tmp_df_trans2gene,tmp_df_CDS_lengths,tmp_df_jackknife_bin_index,tmp_df_timediff,tmp_df_coverage_ratio,tmp_df_middle_reads],axis=1,join='inner')
		tmp_overall_df.index.name = 'transcript_ID'
		tmp_overall_df.reset_index()
		
		overall_df = pd.concat([overall_df,tmp_overall_df],axis=0,join='outer')
	return overall_df



'''
converts dataframe with jacknife bins
to one with best estimates and errors
for a list of columns
returns new dataframe
CAUTION: an overall timing error of 
9nt is added to every statistical error
'''
def convert_jackknife_to_errors(df):
	#fill return dataframe with best estimates
	ret_df = df[df['jackknife_bin_index']==number_jackknife_bins].filter(['transcript_ID','gene_ID','timediff','CDS_length','coverage_ratio','middle_reads']).copy()

	error_colname = "err_coverage_ratio"

	#add error column
	ret_df[error_colname] = np.nan

	#fill error columns
	for transcript in ret_df.transcript_ID.unique():
		for timediff in time_files_dict:
			tmp_list = df.loc[(df['transcript_ID']==transcript) & (df['timediff']==timediff),'coverage_ratio'].to_list()

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

			if ctr > 1:
				tmp *= (ctr - 1) / ctr
				tmp = np.sqrt(tmp)

			ret_df.loc[(df['transcript_ID']==transcript) & (df['timediff']==timediff),error_colname] = tmp

	ret_df = ret_df[(ret_df['coverage_ratio'] > 0.0) & (ret_df[error_colname] > 0.0)]

	return ret_df



'''
linear fit function
at point x with slope m*m and intercept b
CAUTION: m>=0 enforced!
'''
def linear(x, m, b):
	return -m*m*x + b



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
			if np.isnan(cdy)==0 and cdy > 0.0:
				rx.append(cx)
				ry.append(cy)
				rdy.append(cdy)
		
	return rx,ry,rdy



'''
converts dataframe with timepoints and errors
to one with fitted estimates of translation speed
(automatically adapts fitting range)
for a list of columns and their respective error columns
automatically adapts fitting range if ribosomes would have fallen off
returns new dataframe
'''
def fit_elong_speeds(df):
	#fill return dataframe with best estimates
	min_time = min(time_files_dict,key=time_files_dict.get)
	ret_df = df[df['timediff']==min_time].filter(['transcript_ID','gene_ID','CDS_length']).copy()

	fit_colname = "elong_speed_coverage_ratio"
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
		ty = df.loc[df['transcript_ID']==transcript,'coverage_ratio'].to_list()
		tdy = df.loc[df['transcript_ID']==transcript,'err_coverage_ratio'].to_list()
		x, y, dy = prepare_fitting_data(tx,ty,tdy)

		if len(x) > 2:
			#possibly automatically adapt fitting range
			popt, pcov = curve_fit(linear, x, y, sigma=dy, maxfev=2000)
			perr = np.sqrt(np.diag(pcov))
			tmp_chisqred = chisqred(linear, popt, x, y, dy)

			if tmp_chisqred < cutoff_chisq_red or transcript == 'all':
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
		tmp_jackknife_df = compute_overall_df(jackknife_ind,min_CDS_length=90)
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


if os.path.isfile("locations_with_errors.tsv") == False:
	df_w_err = convert_jackknife_to_errors(all_df)
	df_w_err.to_csv("locations_with_errors.tsv",sep="\t",index=False)
	print("Jackknife errors estimated.")
else:
	df_w_err = pd.read_csv("locations_with_errors.tsv",sep="\t")
	print("Jackknife errors read.")

del all_df
gc.collect()


fit_df = fit_elong_speeds(df_w_err)
fit_df.dropna(subset=["elong_speed_coverage_ratio"],thresh=1,inplace=True)
fit_df.to_csv("elongation_speeds.tsv",sep="\t",index=False)


