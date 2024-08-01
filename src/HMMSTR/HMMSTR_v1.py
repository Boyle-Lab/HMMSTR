#I dont need all of these, remove all unused at the end, many will be encapsulated by other files
#from collections import Counter, defaultdict
import argparse
import os, glob
import gzip
#import numpy as np
import pandas as pd
from time import perf_counter
import multiprocessing as mp
import pickle as pkl

#custom imports
from profile_HMM.profile_HMM import ProfileHMM
from HMMSTR_utils.HMMSTR_utils import *
from process_read.process_read import Process_Read
from GMM_stats.GMM_stats import GMMStats

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
    curr_target: pandas Series
        Row from input tsv corresponding to a single target
    alphabet: list
        list of characters used to define alphabet of HMM
    background: dictionary
        Dictionary of background probabilities to use as background frequencies
    out: str
        output prefix

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

    #write out the transition and emission matrices to pickle files for original scaling protocol CAN BE REMOVED IN FINAL
    # with open(out + "_" + name + "_A.pkl", "wb") as outfile:
    #     pkl.dump(curr_hmm.A, outfile)
    # with open(out + "_" + name + "_E.pkl", "wb") as outfile:
    #     pkl.dump(curr_hmm.E, outfile)
    return
def process_read(header,seq,hmm_file,rev_hmm_file,hidden_states,hidden_states_rev,out,targets, build_pre, mode, cutoff, k, w, use_full_seq, flanking_size):
    '''
    Wrapper function to call all Process_Read methods. This will initialize all operations
    for a given read including (1) alignment (2) target identification and filtering 
    (3) repeat identification with viterbi (4) repeat counting and (5) recording relevant metrics

    Parameters
    ----------------------------------------------------------------------------------------------
    header: str
        header for read from fasta file
    seq: str
        raw read sequence
    hmm_file: str
        path to file containing written profile hmm from build_all
    rev_hmm_file: str
        path to file continaing the reverse strand profile hmm from build_all
    hidden_states: str
        path to hidden states file 
    hidden_states_rev: str
        path to hidden states file for reverse strand
    out: str
        Output prefix for run including full or relative path
    targets: pandas DataFrame
        DataFrame containing input information for each target to search for in reads
    build_pre: str
        Path and prefix for hmm files, used in case targets have already been used and user doesn't want duplicate files
    mode: str
        Sequencing technology used, map-ont for Nanopore, pb for PacBio
    cutoff: int
        Mapq cutoff for mappy aligner to determine read maping location
    
    Returns
    ------------------------------------------------------------------------------------------------
    None, all outputs written to out
    '''
    curr_read = Process_Read(header, seq, cutoff, mode, out, k,w, use_full_seq)
    curr_read.align_mappy() #sets prefix_df and suffix_df fields for read
    # if len(curr_read.suffix_df.name) > 0:
    #     print(curr_read.prefix_df)
    #     print(curr_read.suffix_df) 
    curr_read.assign_targets(targets) #sets target_info for current read, this is a dictionary of keys=target name and values=info, all information needed for the target
    #The check for no assigned targets is in run_vit
    targets_found = curr_read.run_viterbi(hmm_file=hmm_file,rev_hmm_file=rev_hmm_file,hidden_states=hidden_states, rev_states=hidden_states_rev, out=out,build_pre=build_pre, prefix_idx = flanking_size)
    return
def ratio_gmm_stats(row, out, out_count_name,plot_hists,max_peaks, filter_outliers, bootstrap, CI_width, resample_size, allele_specific_plots,allele_specif_CIs): 
    '''
    This is a wrapper function to run stats methods on all targets. This version corresponds to our original stats methods
    which includes normalizing by a perfect match path through the hmm followed by GMM clustering and allele calls

    Args:
        target(pandas series): row of target dataframe
        out(str): output suffix
        out_count_name(str): output suffix for final result tsv

    Returns:
        None, results written out to tsv per target

    '''
    #check existence of count file, may not exist if no reads identified for a given target
    name = row[0]
    out_count_file = out + "_" + name + out_count_name
    #print(out_count_file)
    if os.path.exists(out_count_file) == False:
        print(out_count_file + " does not exist, writing null row for", name,"...")
        # returning null at the beginning of the run may be causing issues with apply, try returning a null row of the right dimensions
        curr_dict = {"name":name}
        for i in range(1,max_peaks+1):
            #curr_mean = "H"+str(i)+":mean"
            curr_median = "H"+str(i)+":median"
            curr_mode = "H"+str(i)+":mode"
            curr_sd = "H" + str(i) + ":SD"
            curr_support = "H" + str(i) + ":supporting_reads"
            #if curr_mean not in curr_dict.keys():
             #   curr_dict[curr_mean] = 0
            if curr_median not in curr_dict.keys():
                curr_dict[curr_median] = 0
            if curr_mode not in curr_dict.keys():
                curr_dict[curr_mode] = 0
            if curr_sd not in curr_dict.keys():
                curr_dict[curr_sd] = 0
            if curr_support not in curr_dict.keys():
                curr_dict[curr_support] = 0
            #curr_dict["H"+str(i)] = -1
            curr_dict["num_supporting_reads"] = 0

            #check if we were suppose to bootstrap for this run so we can return the correct number of columns
            if bootstrap:
                curr_dict["H"+str(i)+":median_CI"] = None
            if allele_specif_CIs:
                curr_dict["H"+ str(i)+":median_CI_allele_specific"] = None
        return pd.Series(curr_dict)
    #initialize GMMStats object
    gmm_stats = GMMStats(target_row=row) #contains all target attributes as well as E and A dictionaries
    final_data = gmm_stats.get_stats(out_count_file, out,plot_hists, filter_outliers)

    final_data.to_csv(out +"_"+ gmm_stats.name +"_read_info.tsv", index = False, sep="\t") #previously final_out.tsv

    #peak calling
    final_data = final_data[final_data.outlier == False].copy()
    curr_row, final_data['cluster_assignments'] = gmm_stats.call_peaks(final_data, out, max_peaks, plot=plot_hists, save_allele_plots = allele_specific_plots)
    if curr_row is None:
        print("current row doesnt exist")
        return
    
    if bootstrap:
        curr_row, cluster_assignments_bootstrap = gmm_stats.bootstrap_gmm(curr_row,final_data[final_data.outlier == False], resample_size, CI_width, max_peaks, out)
    if allele_specif_CIs:
        for assignment in final_data['cluster_assignments'].unique():
            curr_row["H"+ str(assignment+1)+":median_CI_allele_specific"] = gmm_stats.bootstrap_gmm_allele_specific(final_data[final_data.cluster_assignments == assignment],resample_size,CI_width,out)

    #write out cluster assignments to file
    assignments = final_data[['read_id','counts']]
    assignments["Name"] = gmm_stats.name
    assignments['assignment'] = final_data['cluster_assignments']
    
    #check if this is the first entry, if it is, save the header, if its not, append to existing file
    if os.path.exists(out + "read_assignments_gmm.tsv"):
        assignments[['Name','read_id','counts','assignment']].to_csv(out + "read_assignments_gmm.tsv",header=None, sep="\t", index=False, mode="a")
    else:
        assignments[['Name','read_id','counts','assignment']].to_csv(out + "read_assignments_gmm.tsv", sep="\t", index=False, mode="a")

    return curr_row #return genotype


def main():
    print("HMMSTR script started")
    parser = argparse.ArgumentParser()
    #Required inputs, TODO subcommands for input types, test this
    subparsers = parser.add_subparsers(dest='subcommand')
    subparsers.required = True
        #  subparser for dump
    parser_infile = subparsers.add_parser('targets_tsv')
    # add a required argument
    parser_infile.add_argument(dest ="targets",type=str, help='TSV with name, prefix, repeat, suffix, header required. Recommended at least 100bp with default alignment parameters')
    #coordinates mode parser
    parser_coords = subparsers.add_parser('coordinates')
    parser_coords.add_argument(dest = 'coords_file',type=str,help='Path to input bed file with either 4 (chr,start,end,repeat motif) or 5 (chr,start,end,repeat motif, name) columns')
    parser_coords.add_argument(dest='chrom_sizes_file',type=str,help = 'Path to chromosome sizes file')
    parser_coords.add_argument(dest="ref",type=str,help="Path to reference genome to get flanking sequence from")
    parser_coords.add_argument("--input_flank_length",type=int,help="Length of prefix and suffix to get from reference, must be longer than 30bp (Default) or flanking_size parameter (default: %(default)s)",default=200)
    #parser.add_argument(dest ="targets",type=str, help='TSV with name, prefix, repeat, suffix, header required. Recommended at least 100bp with default alignment parameters')
    parser.add_argument(dest ="out",type=str, help='Output prefix including directory path')
    parser.add_argument(dest = "inFile", type=str, help= 'Sequence to search and annotate in fasta or fastq, gzipped accepted')

    #coordinate input options ***these inputs are dependent on other inputs and I need to figure out how to deal with that
    #TODO figure out how to make arguments conditional ie how do I make it so you can pass either coordinates or target sheet? do I need a bunch of conditionals to control this or is there a
    #streamlined way to do this?

    #optional, if input, models and accompanying files already produced and use the provided prefix
    parser.add_argument("--background", type=str, help="Pickle of custom background dictionary, must include empty string for deletion character and total must add up to 1", default=None)
    parser.add_argument("--E_probs", type=str, help= 'TSV of custom emissions probability matrix, see format specifications', default=None)
    parser.add_argument("--A_probs", type=str, help= 'TSV of custom transition probabilities, see format specifications', default=None)
    parser.add_argument("--custom_RM", type=str, help= 'TSV of custom repeat match state probability matrix, see format specifications, used for motif mosacism, currently supports single target', default=None)
    parser.add_argument("--hmm_pre", type=str,help="Prefix for files produced by build function, use if running the same targets across multiple input files") #if I clean up the directory after I will get rid of this potentially
    #defaults subject to change
    parser.add_argument("--output_hist",help="Output supporting read histogram", action='store_true')
    #parser.add_argument("--exclude_sd", type=int,help="Number of standard deviations from the mean of tallest peaks to exclude neighboring peaks", default=2)
    parser.add_argument("--max_peaks", type=int,help="Maximum number of peaks to calculate AIC and BIC for in peak calling (default: %(default)s)",default=2)
    parser.add_argument("--cpus", type=int,help="Number of cpus to use. If none given, half of available cpus used (default: %(default)s)",default=int(mp.cpu_count()/2))
    parser.add_argument("--flanking_size", type=int,help="Number of basepairs flanking repeat to encode in model, longer recommended for highly repeatitive regions (default: %(default)s)",default=30)
    #alignment parameters
    parser.add_argument("--mode", type=str,help="map-ont (Nanopore), pb (PacBio), or sr (short accurate reads, use for accurate short flanking sequence input) (default: %(default)s)",default="map-ont")
    parser.add_argument("--cutoff", type=int,help="MapQ cutoff for prefix and suffix alignment (default: %(default)s)",default=30)
    parser.add_argument("--k", type=int,help="Kmer parameter to be passed to mappy",default=None)
    parser.add_argument("--w", type=int,help="Window parameter to be passed to mappy",default=None)
    parser.add_argument("--use_full_read", help="Flag to not subset the read by the alignment, optimal for reads with short prefix and suffix",action='store_true')

    #optional boostrapping parameters
    parser.add_argument("--bootstrap", help="Boolean designating to output boostrapped confidence intervals with genotype calls",action='store_true')
    parser.add_argument("--call_width", type=float, help="Decimal percentage designating width of confidence interval for allele calls (default: %(default)s)", default=0.95)
    parser.add_argument("--resample_size", type=int, help="Number of times to resample in bootstrapping (default: %(default)s)", default=500)

    #allele specific stats options
    parser.add_argument("--allele_specific_CIs",  help="Output allele-specific bootsrapped confidence intervals", action='store_true')
    parser.add_argument("--allele_specific_plots", help="Output allele-specific histograms with model of best fit", action='store_true')
    parser.add_argument("--discard_outliers", help="Discard outliers based on quantile", action='store_true')
    parser.add_argument("--filter_quantile", type=float, help="Float designating quantile of count frequency to discard when filtering outliers (default: %(default)s)", default=0.25)
    #for debugging
    parser.add_argument("--cluster_only", help="Run clustering on output corresponding to required arguments, mostly for debugging past runs", action='store_true')
    parser.add_argument("--save_intermediates", help="Flag designating to save intermediate files such as model inputs, raw count files, and state sequence files", action='store_true')

    # New arguments for motif composition analysis
    parser.add_argument("--motif_comp", action='store_true', help="Flag to run motif composition analysis")
    parser.add_argument("--targets", type=str, nargs='+', help="List of targets for motif composition")

    args = parser.parse_args()

    #check input mode and check if inputs compatible for coordinates
    if args.subcommand == 'coordinates':
        #check inputs
        if args.input_flank_length < args.flanking_size:
            print("Prefix and suffix length (input_flank_length,%d) are incompatible with model flanking length (flanking_size, %d), exiting..."%(args.input_flank_length,args.flanking_size))
            return
        #generate infile from coordinates given
        print("Coordinates input selected, getting %dbp flanking sequence from reference at %s..."%(args.input_flank_length,args.ref))
        targets = generate_input_sheet(args.coords_file, args.chrom_sizes_file, args.ref, args.input_flank_length)
        #TODO decide if this should be an option
        print("Inputs generated from coordinates! Saving target inputs to %s..."%(args.out+"_inputs.tsv"))
        targets.to_csv(args.out + "_inputs.tsv", sep="\t",index=False)
    else:
        #read in targets
        targets = pd.read_csv(args.targets, sep="\t")
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
        print(custom_RM)
    else:
        custom_RM = None

    #default alphabet that is compatible with all methods, subject to change if users want to have additional inputs such as methylated bases ect in the future
    alphabet = ['A','T','C','G','']
    if args.cluster_only == False:
        if args.hmm_pre is None:
            pool_start = perf_counter()
            targets.apply(build_all,axis=1,args=(args.out, background, alphabet,A_probs, E_probs, custom_RM, args.flanking_size)) #added flanking size
            pool_end = perf_counter()
            print("build_profile finished .... time was: ", str(pool_end-pool_start))
            build_pre = args.out #need to make sure this is compatible with new build all implentation
        else:
            print("Using previously written input files with prefix: " + args.hmm_pre)
            build_pre = args.hmm_pre
        convert_to_fasta(targets,args.out) #write out all target flanking sequences to fasta files to be used by mappy

        #file prefixes and suffixes to use
        hmm_file = ".hmm" #this is a suffix, does not include individual names****
        rev_hmm_file = "_revcomp" + ".hmm"
        hidden_states = ".hidden_states.txt"
        hidden_states_rev = "_revcomp.hidden_states.txt"

        #start multiprocess of read processing
        pool_start = perf_counter()
        print("cpu's in use: " + str(args.cpus))
        pool = mp.Pool(processes=args.cpus)

        #check input file type and pass to appropriate parser
        #check if compressed
        if args.inFile.endswith('gz'):
            #gzipped file
            if args.inFile.endswith("fasta.gz") or args.inFile.endswith("fa.gz"): #fasta file
                [pool.apply_async(process_read, args=(header,seq,hmm_file,rev_hmm_file,hidden_states,hidden_states_rev,args.out,targets, build_pre, args.mode, args.cutoff, args.k, args.w, args.use_full_read, args.flanking_size)) for header, seq, bool in read_fasta(gzip.open(args.inFile,'rt'))]
            elif args.inFile.endswith("fastq.gz") or args.inFile.endswith("fq.gz"): #fastq
                [pool.apply_async(process_read, args=(header,seq,hmm_file,rev_hmm_file,hidden_states,hidden_states_rev,args.out,targets, build_pre, args.mode, args.cutoff, args.k, args.w, args.use_full_read, args.flanking_size)) for header, seq, bool in read_fastq(gzip.open(args.inFile,'rt'))]
        else:
            if args.inFile.endswith('a'): #fasta file
                [pool.apply_async(process_read, args=(header,seq,hmm_file,rev_hmm_file,hidden_states,hidden_states_rev,args.out,targets, build_pre, args.mode, args.cutoff, args.k, args.w, args.use_full_read,args.flanking_size)) for header, seq, bool in read_fasta(open(args.inFile))]
            elif args.inFile.endswith('q'): #fastq
                [pool.apply_async(process_read, args=(header,seq,hmm_file,rev_hmm_file,hidden_states,hidden_states_rev,args.out,targets, build_pre, args.mode, args.cutoff, args.k, args.w, args.use_full_read,args.flanking_size)) for header, seq, bool in read_fastq(open(args.inFile))]
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
    #call original allele call procedure
    geno_df = targets.apply(ratio_gmm_stats,axis=1,args=(args.out,out_count_name, args.output_hist,args.max_peaks,args.discard_outliers, args.bootstrap, args.call_width, args.resample_size,args.allele_specific_plots,args.allele_specific_CIs))
    #geno_df.dropna(axis="rows",how="any", inplace=True)
    if isinstance(geno_df, pd.DataFrame):
        geno_df.to_csv(args.out + "genotype_calls_gmm.tsv",sep="\t",index=False)
        pool_end = perf_counter()
        print("GMM run done! Took: ", pool_end-pool_start)
    else:
        print(geno_df)
        print("Results are not a DataFrame! Something went wrong...")

    #clean up intermediates
    #TODO make this an option
    #all files produced for models
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

    from motif_comp_code import load_tandem_repeats_info, extract_tandem_repeats, update_repeat_positions, find_matching_key, concensus_gen, generate_consensus, process_and_run_motifscope, graphing

    if args.motif_comp:
        read_assignments = os.path.join(args.out, f"{args.name}_read_assignments.tsv")
        sample_file = args.inFile
        output_folder = os.path.join(args.out, "motif_comp")
        os.makedirs(output_folder, exist_ok=True)
        motif_targets = args.targets if args.targets else ["all"]

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