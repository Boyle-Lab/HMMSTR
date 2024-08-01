#I dont need all of these, remove all unused at the end, many will be encapsulated by other files
#from collections import Counter, defaultdict
import argparse
import os, glob
import gzip
#import numpy as np
import pandas as pd
from time import perf_counter, sleep
import multiprocessing as mp
from multiprocessing import Pool, Queue, Process, Manager, Semaphore
import pickle as pkl
import sys
#custom imports
from profile_HMM.profile_HMM import ProfileHMM
from HMMSTR_utils.HMMSTR_utils import *
from process_read.process_read import Process_Read
from GMM_stats.GMM_stats import GMMStats
from KDE_stats.KDE_stats import KDE_cluster
#imports for motif plotting
from motif_comp.motif_comp_code import *
import pysam
import mappy
import queue  # imported for using queue.Empty exception


class LoadFromFile (argparse.Action):
    #for reading from input file instead of command line
    def __call__ (self, parser, namespace, values, option_string = None):
        with values as f:
            # parse arguments in the file and store them in the target namespace
            parser.parse_args(f.read().split(), namespace)

def read_my_file(argv):
    '''
    Custom parser for file input
    '''
    new_argv = []
    with open(argv,"r") as params:
        for a in params: #read each file
            for b in a.rstrip().split(" "):
                new_argv.append(b)
    return new_argv
def write_input_command(args):
    '''
    Function to write out args to file compatible with read_my_file for future runs
    '''
    with open(args.out + "_run_input.txt","w") as outfile:
        outfile.write(vars(args)["subcommand"] + "\n")
        #check subcommand 
        if vars(args)["subcommand"] == "coordinates":
            #write out coordinates specific arguments
            outfile.write(vars(args)["coords_file"] + "\n")
            outfile.write(vars(args)["chrom_sizes_file"] + "\n")
            outfile.write(vars(args)["ref"] + "\n")
            outfile.write("--input_flank_length"+" "+str(vars(args)["input_flank_length"]) + "\n")
        else:
            #targets
            outfile.write(vars(args)["targets"] + "\n")
        outfile.write(vars(args)["out"] + "\n")
        outfile.write(vars(args)["inFile"] + "\n")

        #TODO supress bootstrapping options if bootstrapping isn't selected
        for key,val in vars(args).items():
            if key in ["subcommand","out","inFile","coords_file","targets","ref","chrom_sizes_file","input_flank_length"]: #these are positional fields, put these where they go first
                continue
            if val in [None, False]: #check if we want this parameter in our input
                continue
            if val == True: #stores true
                outfile.write("--"+key +"\n")
                continue
            outfile.write("--"+key + " "+ str(val)+"\n")
    return
def convert_to_fasta(tsv,out):
    '''
    Function to write "read" entry for a given target for both prefix and suffix
    fasta files to be used as queries in blastn

    Args:
        tsv (Pandas Dataframe): input tsv containing name, prefix, repeat, suffix
        out (str): prefix/path for outputted files

    Returns:
        nothing, writes prefix and suffix fasta files to out_prefix.fa and out_suffix.fa respectively

    '''
    for row in tsv.itertuples():
        curr_out = open(out+"_prefix.fa",'a')
        curr_out.write(">"+row.name+"\n")
        curr_out.write(row.prefix + "\n")
        curr_out.close()
        curr_out = open(out +"_suffix.fa",'a')
        curr_out.write(">"+row.name+"\n")
        curr_out.write(row.suffix + "\n")
        curr_out.close()
    return
def build_all(curr_target, out, background,alphabet,transitions,mismatch_probs, repeat_probs,flanking_size):
    '''
    Function to initialize a ProfileHMM object for a row in input targets tsv

    Parameters
    ----------------------------------------------------------------------------
    curr_target: pandas Series. Row from input tsv corresponding to a single target
    out: str. Output prefix
    background: dictionary. Dictionary of background probabilities to use as background frequencies
    alphabet: list. List of characters used to define alphabet of HMM
    transitions: dictionary. Dictionary of transition probabilites to encode in models
    mismatch_probs: dictionary. Dictionary of mismatch probabilites to encode as emission probabilites in models
    repeat_probs: dictionary. Dictionary of custom repeat match probabilities
    flanking_size: int. Number of bases to encode in both the prefix and suffix of each model
    Returns
    ------------------------------------------------------------------------------
    None, writes ProfileHMM in format needed for modified umdhmm-v1.02
    '''
    name = curr_target[0]
    prefix = curr_target[1].rstrip()
    repeat = curr_target[2].rstrip()
    suffix = curr_target[3].rstrip()
    #initialize ProfileHMM object
    curr_hmm = ProfileHMM(name=name,alphabet=alphabet, prefix=prefix, repeat=repeat, suffix=suffix, background=background, out=out, transitions=transitions,mismatch_probs=mismatch_probs, repeat_probs=repeat_probs,flanking_size=flanking_size)
    curr_hmm.build_profile("forward") #builds and writes forward ProfileHMM object
    curr_hmm.build_profile("reverse") #builds and writes reverse ProfileHMM object
    return

    
def fetch_reads(bam_file, regions, mapq_cutoff, *args):
    '''
    Function to fetch reads from a BAM file given a list of regions
    '''
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for region in regions:
            for read in bam.fetch(region[0],region[1], region[2]):
                if read.mapping_quality >= mapq_cutoff:
                    seq = mappy.revcomp(read.query_sequence) if read.is_reverse else read.query_sequence
                    read_id = read.query_name
                    #check if name for region is given, if it is, region is name given in regions tuple (4th element), if no 4th element, region is the string of the region
                    if len(region) > 3:
                        region_name = region[3]
                    else:
                        region_name = f'{region[0]}:{region[1]}-{region[2]}'
                    yield (read_id, seq, region_name) + args

def worker(args):
    # print('Worker function called')  # Add this line
    try:
        read_id, seq, region, hmm_file, rev_hmm_file, hidden_states, hidden_states_rev, out, targets, build_pre, mode, mapq_cutoff, k, w, use_full_read, flanking_size, output_labelled_seqs, coords_file = args
        # print('Calling process_read')  # Add this line
        process_read(read_id, seq, hmm_file, rev_hmm_file, hidden_states, hidden_states_rev, out, targets, build_pre, mode, mapq_cutoff, k, w, use_full_read, flanking_size, output_labelled_seqs, region, coords_file)
        return True
    except Exception as e:
        print(f'Exception in worker function: {e}')  # Add this line
def process_read(header,seq,hmm_file,rev_hmm_file,hidden_states,hidden_states_rev,out,targets, build_pre, mode, cutoff, k, w, use_full_seq, flanking_size,output_labelled_seqs,region=None, coord_file=None):
    '''
    Wrapper function to call all Process_Read methods. This will initialize all operations
    for a given read including (1) alignment (2) target identification and filtering 
    (3) repeat identification with viterbi (4) repeat counting and (5) recording relevant metrics

    Parameters
    ----------------------------------------------------------------------------------------------
    header: str. Header for read from fasta file
    seq: str. Raw read sequence
    hmm_file: str. Path to file containing written profile hmm from build_all
    rev_hmm_file: str. Path to file continaing the reverse strand profile hmm from build_all
    hidden_states: str. Path to hidden states file 
    hidden_states_rev: str. Path to hidden states file for reverse strand
    out: str. Output prefix for run including full or relative path
    targets: pandas DataFrame. DataFrame containing input information for each target to search for in reads
    build_pre: str. Path and prefix for hmm files, used in case targets have already been used and user doesn't want duplicate files
    mode: str. Sequencing technology used, map-ont for Nanopore, pb for PacBio
    cutoff: int. Mapq cutoff for mappy aligner to determine read maping location
    k: int. k-mer length parameter to be passed to mappy
    w: int. Minimizer window size to be passed to mappy
    use_full_seq: bool. Override to supress subsetting read sequence before Viterbi
    output_labelled_seqs: bool. Output labelled sequence output by Viterbi to a text file.
    region: str. Region this read mapped to, if known from bam file
    coord_file: str. Path to file containing coordinates for reads to be processed. Applicable in case of bam input.
    
    Returns
    ------------------------------------------------------------------------------------------------
    None, all outputs written to out
    '''
    #check if we are reading input from a bam file or a fastx file
    if coord_file is not None:
        read_id = header
        target = region
    elif coord_file is None:
        read_id = header.split(" ")[0][1:]
        target= None
    curr_read = Process_Read(read_id, seq, cutoff, mode, out, k,w, use_full_seq)
    curr_read.align_mappy(target) #sets prefix_df and suffix_df fields for read

    curr_read.assign_targets(targets) #sets target_info for current read, this is a dictionary of keys=target name and values=info, all information needed for the target
    #The check for no assigned targets is in run_vit
    targets_found = curr_read.run_viterbi(hmm_file=hmm_file,rev_hmm_file=rev_hmm_file,hidden_states=hidden_states, rev_states=hidden_states_rev, out=out,build_pre=build_pre, prefix_idx = flanking_size,output_labelled_seqs=output_labelled_seqs)
    return

def decide_method(target,out,out_count_name):
    '''
    Function to automatically choose a peakcalling method based on data IQR and outliers

    Parameters
    ----------------------------------------------------------------------------------------------
    target: pandas Series. Row from targets dataframe to generate method for
    out: str. Output prefix
    out_count_name: str. Path to file containing raw counts data for this target
    
    Returns
    ------------------------------------------------------------------------------------------------
    decision: str. Method to use for peakcalling
    '''
    #read in counts for current row
    name = target[0]
    out_count_file = out + "_" + name + out_count_name
    if os.path.exists(out_count_file):
        counts_file= pd.read_csv(out_count_file, sep=" ",header=None)
        counts_file.columns = ["read_id","strand","align_score","neg_log_likelihood","subset_likelihood","repeat_likelihood","repeat_start","repeat_end","align_start", "align_end","counts"]
        counts=counts_file.counts
    else:
        return None
    
    QR_range, IQR_range = get_IQR(counts,0.10)

    #number of outliers according to IQR
    lower_counts = len(counts[counts <= IQR_range[0]])
    upper_counts = len(counts[counts >= IQR_range[1]])
    #number of outliers according to quartiles -- not currently used, can be if tuning is necessary in the future
    lower_Q1_counts = len(counts[counts <= QR_range[0]])
    upper_Q3_counts = len(counts[counts <= QR_range[1]])

    data_width = max(counts) - min(counts)
    QR_width = QR_range[1]-QR_range[0]
    IQR_width = IQR_range[1] - IQR_range[0]

    #3 cases: (1) KDE--throw-outliers (2) KDE (3) GMM

    #check if we should use KDE, no outliers found -- can play with data_width definition here if need be
    if data_width < 4 or (QR_width <= 3 and data_width < 2*QR_width): #not including checking for outliers here rn, might need to
        decision = "kde"
    elif QR_width <= 4 and (lower_counts > 0 or upper_counts > 0): #might need to play with this
        decision = "kde_throw_outliers"
    else:
        decision = "gmm"
    return decision

def call_peaks(row, out, out_count_name, plot_hists, max_peaks, filter_outliers=False, filter_quantile=0.25,bootstrap=False, CI_width=0.95, resample_size=100,allele_specific_plots=False,allele_specif_CIs=False, bandwidth='scott',kernel="gaussian",flanking_like_filter=False):
    '''
    This method is the wrapper function for calling both KDE and GMM classes based on the decision per target, it also insures the outputs are uniform across all methods

    Parameters
    ----------------------------------------------------------------------------------------------
    row: pandas Series. Row from targets dataframe to call peaks for
    out: str. Output prefix
    out_count_name: str. Path to file containing raw counts data for this target
    plot_hists: bool. Boolean designating if plots should be saved.
    max_peaks: int. Maximum number of peaks to call.
    filter_outliers: bool. Boolean designating if outlier counts should be filtered before peak calling
    filter_quantile: float. If filter_outliers is set, determine outler counts based on this quantile.
    bootstrap: bool. Boolean designating if bootstrapping of allele calls should be performed.
    CI_width: float. If bootstrapping is set, calculate the confidence intervals at this width.
    resample_size: int. Number of iterations to resample for bootstraping.
    allele_specific_plots: bool. If true, output plots of allele count distributions in addition to overall read histograms and peak calls.
    allele_specif_CIs: bool. Output confidence intervals for each allele independently after clusters are assigned in addition to overall CIs.
    bandwidth: str or float. KDE specific parameter, either a string designating estimation method or float designating bandwidth to use.
    kernel: str. KDE specific parameter, KDE kernel to use in peakcalling
    flanking_like_filter: bool. If set, reads with exceedingly low likelihood at their flanking sequence will be filtered before clustering.
    
    Returns
    ------------------------------------------------------------------------------------------------
    final: pandas Series. Genotype calls and supporting information for final output.
    '''
    #make cluster object based on decision
    decision = row[4]

    if decision is not None and "kde" in decision:
        if decision == "kde_throw_outliers":
            filter_outliers = True #override when this option is chosen, will be the same as kde if user overrides
        curr_kde = KDE_cluster(row,out,out_count_name, filter_outliers,flanking_like_filter) #discard outliers will still be included in all of these regardless of if KDE-throw-outliers is chosen, that is the only option where it is the default
        if curr_kde.data is None: 
            allele_calls = {}
            allele_calls = {"name":curr_kde.name}
            allele_calls['bandwidth'] = -1
            allele_calls["peak_calling_method"] = -1
            allele_calls["num_supporting_reads"] = 0
            for i in range(max_peaks):
                allele_calls["A"+str(i + 1)+":median"] = 0
                allele_calls["A"+str(i + 1)+":mode"] = 0
                allele_calls["A"+str(i + 1)+":supporting_reads"] = 0
                allele_calls["A"+str(i+1)+":SD"] = 0
                if bootstrap:
                    allele_calls["A"+str(i)+":median_CI"] = (0,0)
                if allele_specif_CIs:
                    allele_calls["A"+ str(i+1)+":median_CI_allele_specific"] = (0,0)
            allele_calls_series = pd.Series(allele_calls)
            cols = allele_calls_series.index.tolist()
            cols.sort()
            final = allele_calls_series[cols]
            return final
        #get read info tsv
        curr_kde.data = curr_kde.get_stats(plot_hists, filter_outliers, filter_quantile, flanking_like_filter)
        
        clusters, allele_calls, outliers, flanking_outliers = curr_kde.call_clusters(kernel=kernel, bandwidth=bandwidth,max_k=max_peaks, output_plots =  plot_hists, filter_quantile=filter_quantile)
        assignments = curr_kde.assign_clusters(clusters,outliers, flanking_outliers)

        #plot allele_specific plots if applicable
        if allele_specific_plots:
            for assignment in assignments.cluster_assignments.unique():
                curr_kde.call_clusters(kernel=kernel, bandwidth=bandwidth,max_k=1, output_plots =  True, subset=assignments[assignments.cluster_assignments == assignment], allele_specific=True, allele=assignment, filter_quantile=filter_quantile)
    
        #write out read assignments
        assignments["name"] = curr_kde.name
        assignments["peak_calling_method"] = decision
        
        if flanking_like_filter:
            assignments[['name','read_id','strand', 'align_score','neg_log_likelihood', 'subset_likelihood', 'repeat_likelihood','repeat_start', 'repeat_end', 'align_start', 'align_end','counts', 'freq','cluster_assignments',"outlier","flanking_outlier","peak_calling_method"]].to_csv(out + "_read_assignments.tsv", sep="\t", header=None,index=False, mode="a")
        else:
            assignments[['name','read_id','strand', 'align_score','neg_log_likelihood', 'subset_likelihood', 'repeat_likelihood','repeat_start', 'repeat_end', 'align_start', 'align_end','counts', 'freq','cluster_assignments',"outlier","peak_calling_method"]].to_csv(out + "_read_assignments.tsv", sep="\t", index=False, mode="a", header=None)

        if len(assignments[assignments.outlier == False][assignments.flanking_outlier == False].cluster_assignments.unique()) < max_peaks:
            for i in range(max_peaks - len(assignments[assignments.outlier == False][assignments.flanking_outlier == False].cluster_assignments.unique())):
                j = len(assignments[assignments.outlier == False][assignments.flanking_outlier == False].cluster_assignments.unique())+1
                allele_calls["A"+str(i+j)+":median"] = 0
                allele_calls["A"+str(i+j)+":mode"] = 0
                allele_calls["A"+str(i+j)+":supporting_reads"] = 0
                allele_calls["A"+str(i+j)+":SD"] = 0
                if allele_specif_CIs:
                    allele_calls["A"+str(i+j)+":median_CI_allele_specific"] = (0,0)
        allele_calls["name"] = curr_kde.name
        allele_calls_series = pd.Series(allele_calls)
        #bootstrap
        assignments = pd.merge(left=assignments, right=curr_kde.data[["read_id","freq"]], how="left")
        if bootstrap:
            if clusters != -1:
                print("Starting bootstrap...")
                median_CIs = curr_kde.bootstrap_KDE(assignments[assignments.outlier == False][assignments.flanking_outlier == False], resample_size, CI_width, max_peaks, out)
                allele_calls_series = pd.concat([allele_calls_series,pd.Series(median_CIs)])
        if allele_specif_CIs: 
            for assignment in assignments.cluster_assignments.unique():
                if clusters == -1:
                    break
                if pd.isna(assignment):
                    continue
                print("Starting allele-specific bootstrap...")
                allele_calls_series["A"+ str(assignment)+":median_CI_allele_specific"] = curr_kde.bootstrap_KDE_allele_specific(assignments[assignments.cluster_assignments == assignment], resample_size, CI_width, out)
        
        #get total number of supporting reads
        allele_calls_series["num_supporting_reads"] = 0
        for k in range(max_peaks):
            allele_calls_series["num_supporting_reads"] = allele_calls_series["num_supporting_reads"] + allele_calls["A"+str(k+1) + ":supporting_reads"]
        allele_calls_series["peak_calling_method"] = decision
        cols = allele_calls_series.index.tolist()
        cols.sort()
        final = allele_calls_series[cols]
        return final
        
    else:#GMM is the only other option
        #check existence of count file, may not exist if no reads identified for a given target
        name = row[0]
        out_count_file = out + "_" + name + out_count_name
        #print(out_count_file)
        if os.path.exists(out_count_file) == False:
            print(out_count_file + " does not exist, writing null row for", name,"...")
            # returning null at the beginning of the run may be causing issues with apply, try returning a null row of the right dimensions
            curr_dict = {"name":name}
            curr_dict['bandwidth'] = -1
            curr_dict["peak_calling_method"] = -1
            for i in range(1,max_peaks+1):
                curr_median = "A"+str(i)+":median"
                curr_mode = "A"+str(i)+":mode"
                curr_sd = "A" + str(i) + ":SD"
                curr_support = "A" + str(i) + ":supporting_reads"
                curr_allele_spec = "A"+ str(i)+":median_CI_allele_specific"
                curr_boot = "A"+str(i)+":median_CI"
                if curr_median not in curr_dict.keys():
                    curr_dict[curr_median] = 0
                if curr_mode not in curr_dict.keys():
                    curr_dict[curr_mode] = 0
                if curr_sd not in curr_dict.keys():
                    curr_dict[curr_sd] = 0
                if curr_support not in curr_dict.keys():
                    curr_dict[curr_support] = 0
                if bootstrap and curr_boot not in curr_dict.keys():
                    curr_dict[curr_boot] = (0,0)
                if allele_specif_CIs and curr_allele_spec not in curr_dict.keys():
                    curr_dict[curr_allele_spec] = (0,0)
                curr_dict["num_supporting_reads"] = 0
            
            return pd.Series(curr_dict)
        #initialize GMMStats object
        gmm_stats = GMMStats(target_row=row) #contains all target attributes as well as E and A dictionaries
        final_data = gmm_stats.get_stats(out_count_file, out,plot_hists, filter_outliers, filter_quantile=filter_quantile,flanking_like_filter=flanking_like_filter)
        
        #peak calling
        final_data2 = final_data[final_data.outlier == False][final_data.flanking_outlier == False][final_data.counts != 0].copy()
        curr_row, final_data2['cluster_assignments'] = gmm_stats.call_peaks(final_data2, out, max_peaks, plot=plot_hists, save_allele_plots = allele_specific_plots)
        if curr_row is None:
            print("current row doesnt exist")
            return
        
        if bootstrap:
            print("Starting bootstrap...")
            curr_row = gmm_stats.bootstrap_gmm(curr_row,final_data2, resample_size, CI_width, max_peaks, out)
        if allele_specif_CIs:
            for assignment in final_data2['cluster_assignments'].unique():
                print("Starting allele-specific bootstrap...")
                median_CIs = gmm_stats.bootstrap_gmm_allele_specific(final_data2[final_data2.cluster_assignments == assignment],resample_size,CI_width,out)
                curr_row["A"+ str(assignment+1)+":median_CI_allele_specific"] = median_CIs
            #check if we have a value for every column in max_peaks
            if len(final_data2['cluster_assignments'].unique()) < max_peaks:
                for i in range(max_peaks - len(final_data2['cluster_assignments'].unique())):
                    j = len(final_data2['cluster_assignments'].unique())+1
                    curr_row["A"+str(i+j)+":median_CI_allele_specific"] = (0,0)

        #write out cluster assignments to file
        assignments = pd.merge(left = final_data2[['read_id','counts','cluster_assignments']],right=final_data, on = ["read_id","counts"], how="right")
        assignments["name"] = gmm_stats.name
        assignments["peak_calling_method"] = decision
        assignments["cluster_assignments"] = assignments["cluster_assignments"]+1

        if flanking_like_filter:
            assignments[['name','read_id','strand', 'align_score','neg_log_likelihood', 'subset_likelihood', 'repeat_likelihood','repeat_start', 'repeat_end', 'align_start', 'align_end','counts', 'freq','cluster_assignments',"outlier","flanking_outlier","peak_calling_method"]].to_csv(out + "_read_assignments.tsv", sep="\t", header=None,index=False, mode="a")
        else:
            assignments[['name','read_id','strand', 'align_score','neg_log_likelihood', 'subset_likelihood', 'repeat_likelihood','repeat_start', 'repeat_end', 'align_start', 'align_end','counts', 'freq','cluster_assignments',"outlier","peak_calling_method"]].to_csv(out + "_read_assignments.tsv", sep="\t", index=False, mode="a", header=None)
    #add bandwidth column for kde calls
    curr_row['bandwidth'] = -1
    curr_row["peak_calling_method"] = decision
    return curr_row #return genotype
    
def call_peaks_stranded(row, out, out_count_name, plot_hists, max_peaks, filter_outliers=False, filter_quantile=0.25,bootstrap=False, CI_width=0.95, resample_size=100,allele_specific_plots=False,allele_specif_CIs=False, bandwidth='scott',kernel="gaussian",flanking_like_filter=False, strand=None):
    '''
    This method is the wrapper function for calling both KDE and GMM classes based on the decision per target, it also insures the outputs are uniform across all methods. This version runs only on the given strand.

    Parameters
    ----------------------------------------------------------------------------------------------
    row: pandas Series. Row from targets dataframe to call peaks for
    out: str. Output prefix
    out_count_name: str. Path to file containing raw counts data for this target
    plot_hists: bool. Boolean designating if plots should be saved.
    max_peaks: int. Maximum number of peaks to call.
    filter_outliers: bool. Boolean designating if outlier counts should be filtered before peak calling
    filter_quantile: float. If filter_outliers is set, determine outler counts based on this quantile.
    bootstrap: bool. Boolean designating if bootstrapping of allele calls should be performed.
    CI_width: float. If bootstrapping is set, calculate the confidence intervals at this width.
    resample_size: int. Number of iterations to resample for bootstraping.
    allele_specific_plots: bool. If true, output plots of allele count distributions in addition to overall read histograms and peak calls.
    allele_specif_CIs: bool. Output confidence intervals for each allele independently after clusters are assigned in addition to overall CIs.
    bandwidth: str or float. KDE specific parameter, either a string designating estimation method or float designating bandwidth to use.
    kernel: str. KDE specific parameter, KDE kernel to use in peakcalling
    flanking_like_filter: bool. If set, reads with exceedingly low likelihood at their flanking sequence will be filtered before clustering.
    strand: str. Current strand to call peaks for.
    
    Returns
    ------------------------------------------------------------------------------------------------
    final: pandas Series. Genotype calls and supporting information for final output.
    '''
    #make cluster object based on decision
    decision = row[4]
    if decision is not None and "kde" in decision:
        if decision == "kde_throw_outliers":
            filter_outliers = True #override when this option is chosen, will be the same as kde if user overrides
        curr_kde = KDE_cluster(row,out,out_count_name, filter_outliers,flanking_like_filter,strand) #discard outliers will still be included in all of these regardless of if KDE-throw-outliers is chosen, that is the only option where it is the default
        if curr_kde.data is None:
            allele_calls = {}
            allele_calls = {"name":curr_kde.name}
            allele_calls['bandwidth'] = -1
            allele_calls["peak_calling_method"] = -1
            for i in range(max_peaks):
                allele_calls["A"+str(i + 1)+":median"] = 0
                allele_calls["A"+str(i + 1)+":mode"] = 0
                allele_calls["A"+str(i + 1)+":supporting_reads"] = 0
                allele_calls["A"+str(i+1)+":SD"] = 0
                if bootstrap:
                    allele_calls["A"+str(i)+":median_CI"] = (0,0)
                if allele_specif_CIs:
                    allele_calls["A"+ str(i+1)+":median_CI_allele_specific"] = (0,0)
            allele_calls["strand"] = strand
            allele_calls_series = pd.Series(allele_calls)
            cols = allele_calls_series.index.tolist()
            cols.sort()
            final = allele_calls_series[cols]
            return final #TODO figure out how to return in a stranded manner
        #get read info tsv
        curr_kde.data = curr_kde.get_stats(plot_hists, filter_outliers, filter_quantile, flanking_like_filter, strand=strand)
        
        clusters, allele_calls, outliers, flanking_outliers = curr_kde.call_clusters(kernel=kernel, bandwidth=bandwidth,max_k=max_peaks, output_plots =  plot_hists, filter_quantile=filter_quantile,strand=strand)
        assignments = curr_kde.assign_clusters(clusters,outliers, flanking_outliers)

        #plot allele_specific plots if applicable
        if allele_specific_plots:
            for assignment in assignments.cluster_assignments.unique():
                curr_kde.call_clusters(kernel=kernel, bandwidth=bandwidth,max_k=1, output_plots =  True, subset=assignments[assignments.cluster_assignments == assignment], allele_specific=True, allele=assignment, filter_quantile=filter_quantile, strand=strand)
    
        #write out read assignments
        assignments["name"] = curr_kde.name
        assignments["peak_calling_method"] = decision
        
        if flanking_like_filter:
            assignments[['name','read_id','strand', 'align_score','neg_log_likelihood', 'subset_likelihood', 'repeat_likelihood','repeat_start', 'repeat_end', 'align_start', 'align_end','counts', 'freq','cluster_assignments',"outlier","flanking_outlier","peak_calling_method"]].to_csv(out + "_"+strand+"_read_assignments.tsv", sep="\t", header=None,index=False, mode="a")
        else:
            assignments[['name','read_id','strand', 'align_score','neg_log_likelihood', 'subset_likelihood', 'repeat_likelihood','repeat_start', 'repeat_end', 'align_start', 'align_end','counts', 'freq','cluster_assignments',"outlier","peak_calling_method"]].to_csv(out + "_"+strand+"_read_assignments.tsv", sep="\t", index=False, mode="a", header=None)
        #account for differing column number for less than max_k peaks
        if len(assignments[assignments.outlier == False].cluster_assignments.unique()) < max_peaks:
            for i in range(max_peaks - len(assignments[assignments.outlier == False].cluster_assignments.unique())):
                allele_calls["A"+str(i+len(assignments[assignments.outlier == False].cluster_assignments.unique())+1)+":median"] = 0
                allele_calls["A"+str(i+len(assignments[assignments.outlier == False].cluster_assignments.unique())+1)+":mode"] = 0
                allele_calls["A"+str(i+len(assignments[assignments.outlier == False].cluster_assignments.unique())+1)+":supporting_reads"] = 0
                allele_calls["A"+str(i+len(assignments[assignments.outlier == False].cluster_assignments.unique())+1)+":SD"] = 0
                if allele_specif_CIs:
                    allele_calls["A"+str(i+len(assignments[assignments.outlier == False].cluster_assignments.unique())+1)+":median_CI_allele_specific"] = (0,0)
        allele_calls["name"] = curr_kde.name
        allele_calls["strand"] = strand
        allele_calls_series = pd.Series(allele_calls)

        #bootstrap
        assignments = pd.merge(left=assignments, right=curr_kde.data[["read_id","freq"]], how="left")
        if bootstrap:
            if clusters != -1:
                print("Starting bootstrap...")
                median_CIs = curr_kde.bootstrap_KDE(assignments[assignments.outlier == False], resample_size, CI_width, max_peaks, out)
                allele_calls_series = pd.concat([allele_calls_series,pd.Series(median_CIs)])
        if allele_specif_CIs: 
            for assignment in assignments.cluster_assignments.unique():
                if clusters == -1:
                    break
                if pd.isna(assignment):
                    continue
                print("Starting allele-specific bootstrap...")
                allele_calls_series["A"+ str(assignment)+":median_CI_allele_specific"] = curr_kde.bootstrap_KDE_allele_specific(assignments[assignments.cluster_assignments == assignment], resample_size, CI_width, out)
        
        #get total number of supporting reads
        allele_calls_series["num_supporting_reads"] = 0
        for k in range(max_peaks):
            allele_calls_series["num_supporting_reads"] = allele_calls_series["num_supporting_reads"] + allele_calls["A"+str(k+1) + ":supporting_reads"]
        allele_calls_series["peak_calling_method"] = decision
        cols = allele_calls_series.index.tolist()
        cols.sort()
        final = allele_calls_series[cols]
        return final 
    
    else:#GMM is the only other option
        #check existence of count file, may not exist if no reads identified for a given target
        name = row[0]
        out_count_file = out + "_" + name + out_count_name
        if os.path.exists(out_count_file) == False:
            print(out_count_file + " does not exist, writing null row for", name,"...")
            # returning null at the beginning of the run may be causing issues with apply, try returning a null row of the right dimensions
            curr_dict = {"name":name}
            curr_dict['bandwidth'] = -1
            curr_dict["peak_calling_method"] = -1
            for i in range(1,max_peaks+1):
                curr_median = "A"+str(i)+":median"
                curr_mode = "A"+str(i)+":mode"
                curr_sd = "A" + str(i) + ":SD"
                curr_support = "A" + str(i) + ":supporting_reads"
                #TODO add column names for alleles here, right now all of these are being set to (0,0) when they shouldn't be, somehow they are being reset to nan
                curr_allele_spec = "A"+ str(i)+":median_CI_allele_specific"
                curr_boot = "A"+str(i)+":median_CI"
                if curr_median not in curr_dict.keys():
                    curr_dict[curr_median] = 0
                if curr_mode not in curr_dict.keys():
                    curr_dict[curr_mode] = 0
                if curr_sd not in curr_dict.keys():
                    curr_dict[curr_sd] = 0
                if curr_support not in curr_dict.keys():
                    curr_dict[curr_support] = 0
                if bootstrap and curr_boot not in curr_dict.keys():
                    curr_dict[curr_boot] = (0,0)
                if allele_specif_CIs and curr_allele_spec not in curr_dict.keys():
                    curr_dict[curr_allele_spec] = (0,0)
                curr_dict["num_supporting_reads"] = 0
            curr_dict["strand"] = strand
            return pd.Series(curr_dict)
        #initialize GMMStats object
        gmm_stats = GMMStats(target_row=row) #contains all target attributes as well as E and A dictionaries
        final_data = gmm_stats.get_stats(out_count_file, out,plot_hists, filter_outliers, filter_quantile=filter_quantile,flanking_like_filter=flanking_like_filter, curr_strand=strand)
    
        final_data.to_csv(out +"_"+ gmm_stats.name +"_"+strand+"_final_out.tsv", index = False, sep="\t")
        
        #peak calling
        final_data2 = final_data[final_data.outlier == False][final_data.flanking_outlier == False][final_data.counts != 0].copy()
        curr_row, final_data2['cluster_assignments'] = gmm_stats.call_peaks(final_data2, out, max_peaks, plot=plot_hists, save_allele_plots = allele_specific_plots, strand=strand)
        if curr_row is None:
            print("current row doesnt exist")
            return
        
        if bootstrap:
            print("Starting bootstrap...")
            curr_row = gmm_stats.bootstrap_gmm(curr_row,final_data2, resample_size, CI_width, max_peaks, out)
        if allele_specif_CIs:
            for assignment in final_data2['cluster_assignments'].unique():
                median_CIs = gmm_stats.bootstrap_gmm_allele_specific(final_data2[final_data2.cluster_assignments == assignment],resample_size,CI_width,out)
                curr_row["A"+ str(assignment+1)+":median_CI_allele_specific"] = median_CIs
            if len(final_data2['cluster_assignments'].unique()) < max_peaks:
                for i in range(max_peaks - len(final_data2['cluster_assignments'].unique())):
                    j = len(final_data2['cluster_assignments'].unique())+1
                    curr_row["A"+str(i+j)+":median_CI_allele_specific"] = (0,0)
    
        #write out cluster assignments to file
        assignments = pd.merge(left = final_data2[['read_id','counts','cluster_assignments']],right=final_data, on = ["read_id","counts"], how="right")
        assignments["name"] = gmm_stats.name
        assignments["peak_calling_method"] = decision
        assignments["cluster_assignments"] = assignments["cluster_assignments"]+1
        if flanking_like_filter:
            assignments[['name','read_id','strand', 'align_score','neg_log_likelihood', 'subset_likelihood', 'repeat_likelihood','repeat_start', 'repeat_end', 'align_start', 'align_end','counts', 'freq','cluster_assignments',"outlier","flanking_outlier","peak_calling_method"]].to_csv(out + "_"+strand+"_read_assignments.tsv", sep="\t", header=None,index=False, mode="a")
        else:
            assignments[['name','read_id','strand', 'align_score','neg_log_likelihood', 'subset_likelihood', 'repeat_likelihood','repeat_start', 'repeat_end', 'align_start', 'align_end','counts', 'freq','cluster_assignments',"outlier","peak_calling_method"]].to_csv(out + "_"+strand+"_read_assignments.tsv", sep="\t", index=False, mode="a", header=None)
    #add bandwidth column for kde calls
    curr_row['bandwidth'] = -1
    curr_row["peak_calling_method"] = decision
    curr_row["strand"] = strand
    return curr_row


def main():
    parser = argparse.ArgumentParser()
    #Required inputs
    subparsers = parser.add_subparsers(dest='subcommand')
    subparsers.required = True
    #Target input tsv parser
    parser_infile = subparsers.add_parser('targets_tsv')
    parser_infile.add_argument(dest ="targets",type=str, help='TSV with name, prefix, repeat, suffix, header required. Recommended at least 100bp with default alignment parameters')
    #coordinates mode parser
    parser_coords = subparsers.add_parser('coordinates')
    parser_coords.add_argument(dest = 'coords_file',type=str,help='Path to input custom bed file with either 4 (chr,start,end,repeat motif) or 5 (chr,start,end,repeat motif, name) columns')
    parser_coords.add_argument(dest='chrom_sizes_file',type=str,help = 'Path to chromosome sizes file')
    parser_coords.add_argument(dest="ref",type=str,help="Path to reference genome to get flanking sequence from")
    parser_coords.add_argument("--input_flank_length",type=int,help="Length of prefix and suffix to get from reference, must be longer than 30bp (Default) or flanking_size parameter (default: %(default)s)",default=200)

    #same as above but for targets_tsv
    parser.add_argument(dest ="out",type=str, help='Output prefix including directory path')
    parser.add_argument(dest = "inFile", type=str, help= 'Sequence to search and annotate in fasta or fastq, gzipped accepted. Bam files accepted in coordinate mode.')

    #Optional output options
    parser.add_argument("--output_plots",help="Output plots, default outputs supporting read histograms and model of best fit plots per target", action='store_true')
    parser.add_argument("--output_labelled_seqs",help="Output the model path through prefix, repeat, and suffix identified per read as context_labelled.txt per target", action='store_true')
    parser.add_argument("--max_peaks", type=int,help="Maximum number of peaks to calculate AIC and BIC for in peak calling (default: %(default)s)",default=2)
    parser.add_argument("--cpus", type=int,help="Number of cpus to use. If none given, half of available cpus used (default: %(default)s)",default=int(mp.cpu_count()/2))
    parser.add_argument("--flanking_size", type=int,help="Number of basepairs flanking repeat to encode in model, longer recommended for highly repeatitive regions (default: %(default)s)",default=100)

    #alignment parameters
    parser.add_argument("--mode", type=str,help="map-ont (Nanopore), pb (PacBio), or sr (short accurate reads, use for accurate short flanking sequence input) (default: %(default)s)",default="map-ont")
    parser.add_argument("--mapq_cutoff", type=int,help="MapQ cutoff for prefix and suffix alignment (default: %(default)s)",default=30)
    parser.add_argument("--k", type=int,help="Kmer parameter to be passed to mappy",default=None)
    parser.add_argument("--w", type=int,help="Window parameter to be passed to mappy",default=None)
    parser.add_argument("--use_full_read", help="Flag to not subset the read by the alignment, optimal for reads with short prefix and suffix",action='store_true')

    #Peakcalling options
    parser.add_argument("--peakcalling_method", type=str, help="If set, this peak calling method will be used instead of chosen by decision process. Options: gmm, kde, or kde_throw_outliers (default: %(default)s)",default="auto")
    parser.add_argument("--bandwidth", type=int, help="Bandwidth to use for KDE (if chosen). Will default to scott method if not provided", default=None)
    parser.add_argument("--kernel", type=str, help="Kernel to use for KDE (if chosen) (default: %(default)s)", default="gaussian")

    #optional boostrapping parameters
    parser.add_argument("--bootstrap", help="Boolean designating to output boostrapped confidence intervals with genotype calls",action='store_true')
    parser.add_argument("--call_width", type=float, help="Decimal percentage designating width of confidence interval for allele calls (default: %(default)s)", default=0.95)
    parser.add_argument("--resample_size", type=int, help="Number of times to resample in bootstrapping (default: %(default)s)", default=100)

    #allele specific stats options
    parser.add_argument("--allele_specific_CIs",  help="Output allele-specific bootsrapped confidence intervals", action='store_true')
    parser.add_argument("--allele_specific_plots", help="Output allele-specific histograms with model of best fit", action='store_true')
    parser.add_argument("--discard_outliers", help="Discard outliers based on quantile", action='store_true')
    parser.add_argument("--filter_quantile", type=float, help="Float designating quantile of count frequency to discard when filtering outliers (default: %(default)s)", default=0.25)

    #post read processing filters/options
    parser.add_argument("--flanking_like_filter", help="If set, filter reads per target that are designated as outliers based on flanking sequence likelihood",action='store_true')
    parser.add_argument("--stranded_report", help="If set, genotypes are called for each strand separately as well as together and strand bias is reported if found",action='store_true')

    #for debugging
    parser.add_argument("--cluster_only", help="Run clustering on output corresponding to required arguments, can only be run after a run with --save_intermediates", action='store_true')
    parser.add_argument("--save_intermediates", help="Flag designating to save intermediate files such as model inputs, raw count files, and state sequence files", action='store_true')

    #optional, if input, models and accompanying files already produced and use the provided prefix
    parser.add_argument("--background", type=str, help="Pickle of custom background dictionary, must include empty string for deletion character and total must add up to 1", default=None)
    parser.add_argument("--E_probs", type=str, help= 'TSV of custom emissions probability matrix, see format specifications', default=None)
    parser.add_argument("--A_probs", type=str, help= 'TSV of custom transition probabilities, see format specifications', default=None)
    parser.add_argument("--custom_RM", type=str, help= 'TSV of custom repeat match state probability matrix, see format specifications, used for motif mosacism, single target only', default=None)
    parser.add_argument("--hmm_pre", type=str,help="Prefix for files produced by build function, use if running the same targets across multiple input files. Only compatible with --save_intermediates option from previous run")


    # New arguments for motif composition analysis
    parser.add_argument("--motif_comp", action='store_true', help="Flag to run motif composition analysis")
    parser.add_argument("--motif_targets", type=str, nargs='+', help="List of targets for motif composition")

    #check input to see if file or command line was used
    if len(sys.argv) != 2:
        args = parser.parse_args()
    elif len(sys.argv) == 2:
        if sys.argv[1] == "--help" or sys.argv[1] == "-h":
            args = parser.parse_args()
            sys.exit()
        args = parser.parse_args(read_my_file(sys.argv[1]))

    print("Parameters for this HMMSTR run: ")
    with open(args.out + "_run_parameters.txt","w") as outfile:
        for key,val in vars(args).items():
            print(key,":",val)
            outfile.write(key + ":"+ str(val)+"\n")
    write_input_command(args)

    if args.output_plots:
        directory = args.out + "_plots"
        if os.path.exists(directory) == False:
            os.mkdir(directory) 
            print("Directory '% s' created" % directory)
        else:
            print("Directory '% s' already exists! Plots will be output to existing directory..." % directory)
    if args.output_labelled_seqs:
        directory = args.out + "_labelled_seqs"
        if os.path.exists(directory) == False:
            os.mkdir(directory) 
            print("Directory '% s' created" % directory)
        else:
            print("Directory '% s' already exists! Labelled sequences will be output to existing directory..." % directory) #if running on the same targets as a previous run the labelled sequences will be appended to the current files
    #check input mode and check if inputs compatible for coordinates
    if args.subcommand == 'coordinates':
        #check inputs
        if args.input_flank_length < args.flanking_size:
            print("Prefix and suffix length (input_flank_length,%d) are incompatible with model flanking length (flanking_size, %d), exiting..."%(args.input_flank_length,args.flanking_size))
            return
        #generate infile from coordinates given
        print("Coordinates input selected, getting %dbp flanking sequence from reference at %s..."%(args.input_flank_length,args.ref))
        targets = generate_input_sheet(args.coords_file, args.chrom_sizes_file, args.ref, args.input_flank_length)
        print("Inputs generated from coordinates! Saving target inputs to %s..."%(args.out+"_inputs.tsv"))
        targets.to_csv(args.out + "_inputs.tsv", sep="\t",index=False)


    if os.path.exists(args.inFile) == False:
        print("Infile file does not exist, exiting...")
        return
    if args.subcommand == 'targets_tsv':
        #read in targets
        print("Target tsv input selected! Reading input from  %s..."%(args.targets))
        targets = pd.read_csv(args.targets, sep="\t")
    else:
        print("Targets generated from coordinates. Reading in target inputs from %s..."%(args.out + "_inputs.tsv"))
    #convert custom inputs
    if args.E_probs is not None:
        E_probs = read_model_params(args.E_probs,"emissions")
    else:
        E_probs = None
    if args.A_probs is not None:
        A_probs = read_model_params(args.A_probs,"transitions")
    else:
        A_probs = None
    if args.background is not None:
        background = read_model_params(args.background,"background")
    else:
        background = None
    if args.custom_RM is not None:
        custom_RM = read_model_params(args.custom_RM,"repeat")
        #print(custom_RM)
    else:
        custom_RM = None
    
    if args.bandwidth is None:
        args.bandwidth = 'scott'

    #default alphabet that is compatible with all methods
    alphabet = ['A','T','C','G','']
    if args.cluster_only == False:
        if args.hmm_pre is None:
            pool_start = perf_counter()
            #build all models
            targets.apply(build_all,axis=1,args=(args.out, background, alphabet,A_probs, E_probs, custom_RM, args.flanking_size))
            pool_end = perf_counter()
            print("All models built! finished .... time was: ", str(pool_end-pool_start))
            build_pre = args.out 
        else:
            print("Using previously written input files with prefix: " + args.hmm_pre)
            build_pre = args.hmm_pre
        convert_to_fasta(targets,args.out) #write out all target flanking sequences to fasta files to be used by mappy

        #file prefixes and suffixes to use
        hmm_file = ".hmm" #this is a suffix, does not include individual names
        rev_hmm_file = "_revcomp" + ".hmm"
        hidden_states = ".hidden_states.txt"
        hidden_states_rev = "_revcomp.hidden_states.txt"

        #start multiprocess of read processing
        pool_start = perf_counter()
        print("Starting multiprocess, cpu's in use: " + str(args.cpus))
        pool = mp.Pool(processes=args.cpus)

        #check input file type and pass to appropriate parser
        #bam file input
        if args.subcommand == 'coordinates' and args.inFile.endswith('.bam'): # check if bam infile
            with open(args.coords_file, 'r') as f:
                num_columns = len(f.readline().split())
                if num_columns == 4:
                    regions = [(line.split()[0], int(line.split()[1]), int(line.split()[2])) for line in f] #FIXME this current implementation uses the coordinates from the bed file to get regions and ignores names (if provided), this will break in process read since we write our outputs in terms of names, thus currently only compatible with no names given
                elif num_columns == 5:
                    regions = [(line.split()[0], int(line.split()[1]), int(line.split()[2]),line.split()[4]) for line in f] #names are given and these will correspond to the outputs and model names
            pool.imap(worker, fetch_reads(args.inFile, regions, args.mapq_cutoff, hmm_file, rev_hmm_file, hidden_states, hidden_states_rev, args.out, targets, build_pre, args.mode, args.mapq_cutoff, args.k, args.w, args.use_full_read, args.flanking_size, args.output_labelled_seqs, args.coords_file), chunksize=10) #10 for now, ran faster than args.cpus
        elif args.subcommand == 'targets_tsv' and args.inFile.endswith('bam'): #not currently allowed
            print("Bam file input not allowed with targets_tsv input, please run with coordinates input or fast(a/q) infile. exiting...")
            return
        #check if compressed
        else:
            if args.inFile.endswith('gz'):
                #gzipped file
                if args.inFile.endswith("fasta.gz") or args.inFile.endswith("fa.gz"): #fasta file
                    [pool.apply_async(process_read, args=(header,seq,hmm_file,rev_hmm_file,hidden_states,hidden_states_rev,args.out,targets, build_pre, args.mode, args.mapq_cutoff, args.k, args.w, args.use_full_read, args.flanking_size, args.output_labelled_seqs)) for header, seq, bool in read_fasta(gzip.open(args.inFile,'rt'))]
                elif args.inFile.endswith("fastq.gz") or args.inFile.endswith("fq.gz"): #fastq
                    [pool.apply_async(process_read, args=(header,seq,hmm_file,rev_hmm_file,hidden_states,hidden_states_rev,args.out,targets, build_pre, args.mode, args.mapq_cutoff, args.k, args.w, args.use_full_read, args.flanking_size,args.output_labelled_seqs)) for header, seq, bool in read_fastq(gzip.open(args.inFile,'rt'))]
            else:
                if args.inFile.endswith('a'): #fasta file
                    [pool.apply_async(process_read, args=(header,seq,hmm_file,rev_hmm_file,hidden_states,hidden_states_rev,args.out,targets, build_pre, args.mode, args.mapq_cutoff, args.k, args.w, args.use_full_read,args.flanking_size,args.output_labelled_seqs)) for header, seq, bool in read_fasta(open(args.inFile))]
                elif args.inFile.endswith('q'): #fastq
                    [pool.apply_async(process_read, args=(header,seq,hmm_file,rev_hmm_file,hidden_states,hidden_states_rev,args.out,targets, build_pre, args.mode, args.mapq_cutoff, args.k, args.w, args.use_full_read,args.flanking_size,args.output_labelled_seqs)) for header, seq, bool in read_fastq(open(args.inFile))]
        pool.close()
        pool.join()
        pool_end = perf_counter()
        print("All reads processed, the pooled run took: ", pool_end-pool_start)

    #stats runs and allele calls
    if args.hmm_pre is None:
        build_pre = args.out #need to make sure this is compatible with new build all implentation
    else:
        build_pre = args.hmm_pre
    out_count_name = "_counts.txt"

    pool_start = perf_counter()

    #Choose a peakcalling method for allele calls
    if args.peakcalling_method == "auto":
        targets["peak_call_method"] = targets.apply(decide_method,args=(args.out, out_count_name), axis=1)
    elif args.peakcalling_method in ["kde","kde_throw_outliers","gmm"]:
        targets["peak_call_method"] = args.peakcalling_method
    else:
        print("Invalid method inputed for peak calling, using auto...")
        targets["peak_call_method"] = targets.apply(decide_method,args=(args.out, out_count_name), axis=1)

    #default output
    if args.stranded_report == False:
        #Note: this will override the current read_assignments file if it exists
        if args.flanking_like_filter:
            pd.DataFrame(['name','read_id','strand', 'align_score','neg_log_likelihood', 'subset_likelihood', 'repeat_likelihood','repeat_start', 'repeat_end', 'align_start', 'align_end','counts', 'freq','cluster_assignments',"outlier","flanking_outlier","peak_calling_method"]).T.to_csv(args.out + "_read_assignments.tsv", sep="\t", index=False, header=None)
        else:
            pd.DataFrame(['name','read_id','strand', 'align_score','neg_log_likelihood', 'subset_likelihood', 'repeat_likelihood','repeat_start', 'repeat_end', 'align_start', 'align_end', 'counts','freq','cluster_assignments',"outlier","peak_calling_method"]).T.to_csv(args.out + "_read_assignments.tsv", sep="\t", index=False,header=None)

        #call genotypes with assigned or given peakcalling method
        geno_df = targets.apply(call_peaks, args=(args.out, out_count_name, args.output_plots,args.max_peaks,args.discard_outliers,args.filter_quantile,args.bootstrap, args.call_width, args.resample_size,args.allele_specific_plots,args.allele_specific_CIs, args.bandwidth,args.kernel,args.flanking_like_filter), axis=1)
        #read resulting read_assignments for sorting
        read_assignments = pd.read_csv(args.out + "_read_assignments.tsv", sep="\t")
        #check if valid genotype output
        if isinstance(geno_df, pd.DataFrame) and isinstance(read_assignments, pd.DataFrame):
            #sort outputs
            geno_df_final = geno_df.apply(sort_outputs, args=(args.max_peaks,args.bootstrap,args.allele_specific_CIs), axis=1)
            read_assignments_final = sort_read_assign_alleles(read_assignments, geno_df_final, args.max_peaks)
            read_assignments_final.to_csv(args.out + "_read_assignments.tsv", sep="\t", index=False)
            geno_df_final.to_csv(args.out + "_genotype_calls.tsv",sep="\t",index=False)
            pool_end = perf_counter()
            print("Genotyping run done! Took: ", pool_end-pool_start)
        else:
            print(geno_df)
            print("Results are not a DataFrame! Something went wrong...")
    else: #stranded output
        if args.flanking_like_filter:
            pd.DataFrame(['name','read_id','strand', 'align_score','neg_log_likelihood', 'subset_likelihood', 'repeat_likelihood','repeat_start', 'repeat_end', 'align_start', 'align_end','counts', 'freq','cluster_assignments',"outlier","flanking_outlier","peak_calling_method"]).T.to_csv(args.out + "_forward_read_assignments.tsv", sep="\t", index=False,header=None)
            pd.DataFrame(['name','read_id','strand', 'align_score','neg_log_likelihood', 'subset_likelihood', 'repeat_likelihood','repeat_start', 'repeat_end', 'align_start', 'align_end','counts', 'freq','cluster_assignments',"outlier","flanking_outlier","peak_calling_method"]).T.to_csv(args.out + "_reverse_read_assignments.tsv", sep="\t", index=False,header=None)

        else:
            pd.DataFrame(['name','read_id','strand', 'align_score','neg_log_likelihood', 'subset_likelihood', 'repeat_likelihood','repeat_start', 'repeat_end', 'align_start', 'align_end', 'counts','freq','cluster_assignments',"outlier","peak_calling_method"]).T.to_csv(args.out + "_forward_read_assignments.tsv", sep="\t", index=False,header=None)
            pd.DataFrame(['name','read_id','strand', 'align_score','neg_log_likelihood', 'subset_likelihood', 'repeat_likelihood','repeat_start', 'repeat_end', 'align_start', 'align_end', 'counts','freq','cluster_assignments',"outlier","peak_calling_method"]).T.to_csv(args.out + "_reverse_read_assignments.tsv", sep="\t", index=False,header=None)

        #genotype each strand independently
        geno_forward_df = targets.apply(call_peaks_stranded, args=(args.out, out_count_name, args.output_plots,args.max_peaks,args.discard_outliers,args.filter_quantile,args.bootstrap, args.call_width, args.resample_size,args.allele_specific_plots,args.allele_specific_CIs, args.bandwidth,args.kernel,args.flanking_like_filter,"forward"), axis=1)
        geno_reverse_df = targets.apply(call_peaks_stranded, args=(args.out, out_count_name, args.output_plots,args.max_peaks,args.discard_outliers,args.filter_quantile,args.bootstrap, args.call_width, args.resample_size,args.allele_specific_plots,args.allele_specific_CIs, args.bandwidth,args.kernel,args.flanking_like_filter, "reverse"), axis=1)
        
        #concat stranded results
        geno_df = pd.concat([geno_forward_df,geno_reverse_df]).sort_values(by="name")
        if isinstance(geno_df, pd.DataFrame):
            geno_df_final = geno_df.apply(sort_outputs, args=(args.max_peaks,args.bootstrap,args.allele_specific_CIs,True), axis=1)
            for i in range(1,args.max_peaks+1):
                geno_df_final["A"+str(i)+":strand_std"] = geno_df_final.groupby("name")["A"+str(i)+":median"].transform(np.std)
            geno_df_final.to_csv(args.out + "_genotype_calls.tsv",sep="\t",index=False)
            pool_end = perf_counter()
            print("Genotyping done! Took: ", pool_end-pool_start)
        else:
            print(geno_df)
            print("Results are not a DataFrame! Something went wrong...")
    #clean up intermediates
    if not args.save_intermediates:
        for file_subset in glob.glob(args.out+"*hidden_states.txt"):
            if os.path.exists(file_subset):
                os.remove(file_subset)
        for file_subset in glob.glob(args.out+"*.hmm"):
            if os.path.exists(file_subset):
                os.remove(file_subset)
        for file_subset in glob.glob(args.out+"*labeled_seqs.txt"):
            if os.path.exists(file_subset):
                os.remove(file_subset)
        #FIXME counts file is required for cluster only
        for file_subset in glob.glob(args.out+"*counts.txt"):
            os.remove(file_subset)
        if os.path.exists(args.out+"_prefix.fa"):
            os.remove(args.out+"_prefix.fa")
        if os.path.exists(args.out+"_suffix.fa"):
            os.remove(args.out+"_suffix.fa")

    # Code block to run motif composition pipeline

    if args.motif_comp:
        motif_info = motif_comp.motif_info.tsv
        read_assignments = os.path.join(args.out, f"{args.out}_read_assignments.tsv")
        sample_file = args.inFile
        output_folder = os.path.join(args.out, "motif_comp")
        os.makedirs(output_folder, exist_ok=True)
        motif_targets = args.motif_targets if args.motif_targets else "all"

        # Generate consensus sequences
        generate_consensus(read_assignments, sample_file, output_folder, targets=motif_targets)

        # Define input and output for motifscope
        input_fasta = os.path.join(output_folder, "consensus_sequence_file.fa")
        motif_data = args.targets  

        # Run motifscope processing
        os.chdir(output_folder)
        process_and_run_motifscope(input_fasta, output_folder, motif_data, targets=motif_targets)

        # Generate plots
        for file in os.listdir(output_folder):
            if "compressed_representation" in file:
                target = [file.split(".fa")[0]]
                input_file = os.path.join(output_folder, file)
                graphing(input_file, targets=target)
                os.remove(file)  # Clean intermediates for motif comp
            elif "motif_in_hex" in file:
                os.remove(file)
            elif "motif_to_hex" in file:
                os.remove(file)

if __name__ == "__main__":
  main()