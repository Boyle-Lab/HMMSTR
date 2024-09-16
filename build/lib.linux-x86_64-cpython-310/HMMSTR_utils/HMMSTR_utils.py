import pandas as pd
import mappy
import numpy as np
import pysam
def rev_comp(strand):
    '''
    Get the reverse complement of a given sequence
    '''
    comp = ""
    for i in range(0, len(strand)):
        next_base = strand[len(strand) - i - 1]
        #check what the base is to switch to complement
        if next_base == "A":
            comp = comp + "T"
        if next_base == "T":
            comp = comp + "A"
        if next_base == "G":
            comp = comp + "C"
        if next_base == "C":
            comp = comp + "G"
        if next_base == "N":
            comp = comp + "N"
    return comp

def seq2int(seq):
  '''
  Convert DNA sequence to integer representation for Viterbi algorithm
  '''
  d = {'A': '1', 'T': '2', 'C': '3', 'G': '4'}
  seqInt = ' '.join(d[s] if s in d else s for s in seq.upper())
  return seqInt
  
def read_fasta(fileObject):
  '''
  Generator function to read in reads from fasta
  '''
  header = ''
  seq = ''
  # skip any useless leading information
  for line in fileObject:
    if line.startswith('>'):
      header = line.strip()
      break
  for line in fileObject:
    if line.startswith('>'):
      if 'N' in seq.upper():
        yield header, seq, False
      else:
        yield header, seq, True
      header = line.strip()
      seq = ''
    else:
      seq += line.strip()
  if header:
    if 'N' in seq.upper():
      yield header, seq, False
    else:
      yield header, seq, True

def read_fastq(fileObject):
  '''
  Generator function to read in reads from fastq
  '''
  header = ''
  seq = ''
  # skip any useless leading information
  for line in fileObject:
    #print(line)
    if line.startswith('@'):
      header = line.strip()
      break
  for line in fileObject:
    if line.startswith('@'):
      if 'N' in seq.upper():
        yield header, seq, False
      else:
        yield header, seq, True
      header = line.strip()
      seq = ''
    #check if optional or quality field
    elif line.startswith('+'):
       #assume quailty is one line and is after optional line
       next(fileObject)
       continue
    else:
      seq += line.strip()
  if header:
    if 'N' in seq.upper():
      yield header, seq, False
    else:
      yield header, seq, True

def remove_outlier_IQR(df,quantile=0.25):
    '''
    Removes outlier repeat copy numbers from distribution
    '''
    Q1=df.quantile(quantile)
    Q3=df.quantile(1-quantile)
    IQR=Q3-Q1
    df_final=df[~((df<(Q1-1.5*IQR)) | (df>(Q3+1.5*IQR)))]
    outliers= df[~(df.isin(df_final))]
    return df_final, outliers

def throw_low_cov(df):
    '''
    Filters any copy number with less than one supporting read
    '''
    filtered = df[df.freq > 1]
    return filtered

def get_IQR(df,quantile=0.25):
    '''
    Calculates the inter quantile range of the supporting read distribution
    '''
    Q1=df.quantile(quantile)
    Q3=df.quantile(1-quantile)
    IQR=Q3-Q1
    return (Q1, Q3), (Q1-1.5*IQR, Q3+1.5*IQR)

#input adapters for probabilities
def read_model_params(input_file, param_type):
   '''
   Converts user tsv input of custom model paramaters to dictionaries compatible with the
   original implementation
   Args:
      input_file (str): file name of input parameter tsv
      param_type (str): string designating which parameter dictionary to convert to, options: transitions, emissions, repeats, background
    Returns:
      param_dict (dictionary): dictionary of designated parameters
      
   '''
   #convert to respective dictionary
   if param_type in ["transitions","background"]: #input_df is a series here
    #convert to transition dictionary
    input_df = pd.read_csv(input_file, sep="\t")
    param_dict = input_df.iloc[0].to_dict()
    if param_type == 'background':
      param_dict[''] = 0
   if param_type in ["emissions","repeat"]: #input_df is a dataframe here
    #convert to transition dictionary
    input_df = pd.read_csv(input_file, sep="\t",index_col=0)
    if param_type == 'repeat':
       #add blank row
      input_df.loc[len( pd.DataFrame(input_df).index)] = 0
      input_df.index = ['A','G','C','T','']
    param_dict = input_df.to_dict()
   return param_dict

def generate_input_sheet(coords_file, chrom_sizes_file, ref,flanking_length=200):
  '''
  Function to generate input sheet with the following info from 4 or 5 column bedfile with the following columns: chromosome, start, end, repeat motif, name
  (name is an optional column with default name being the coordinated in chr:start-end format)
  Args:
    coords_file (str): bed file with 4 or 5 columns in chromosome, start, end, repeat motif, name format
    chrom_sizes_file (str): chrom.sizes file containing the sizes of each chromosome in the reference, used to ensure prefix/suffix does not exceed the length of the chromosome
    ref (str): path to reference genome to use
    flanking_length (int): length of desired prefix and suffix, must be greater than the length of flank parameter (if given) or 30 (if default), recommended >100bp
  Returns:
    targets (DataFrame): pandas dataframe of targets with columns name, prefix, repeat, suffix
  '''
  #load reference genome
  ref_aligner = mappy.Aligner(ref)
  #load chromosome sizes
  chrom_sizes = pd.read_csv(chrom_sizes_file, sep="\t",names=["chr","chr_length"])
  #load bedfile
  coords=pd.read_csv(coords_file,sep="\t",header=None)
  #check if name column is included
  if len(coords.columns) == 5:
    coords.columns =  ["chr","start","end","repeat","name"]
  elif len(coords.columns) == 4:
    #add name column
    coords.columns = ["chr","start","end","repeat"]
    coords["name"] = coords.chr + ":" +coords.start.astype(str) + "-" + coords.end.astype(str)
  else:
    print("Bed file has unexpected number of columns, please ensure file has 4 (repeat motif) or 5 (motif and name) columns. Exiting...")
    return
  #get prefix and suffix coordinates
  coords['prefix_start'] = coords.start - flanking_length
  coords['suffix_end'] = coords.end + flanking_length
  #replace invalid flanking starts and ends
  coords_final = pd.merge(left=coords,right=chrom_sizes,how="left",on="chr")
  coords_final['prefix_start'] = np.where(coords_final['prefix_start'] < 0,0,coords_final['prefix_start'])
  coords_final['suffix_end'] = np.where(coords_final['suffix_end'] > coords_final['chr_length'], coords_final['chr_length'], coords_final['suffix_end'])

  #set up dictionary
  targets_dict = {"name":[],"prefix":[],"repeat":[],"suffix":[]}
  for target in coords_final.itertuples():
    targets_dict['name'].append(target.name)
    targets_dict["prefix"].append(ref_aligner.seq(name=target.chr,start=target.prefix_start,end=target.start))
    targets_dict['repeat'].append(target.repeat)
    targets_dict['suffix'].append(ref_aligner.seq(name=target.chr,start=target.end,end=target.suffix_end))
  return pd.DataFrame(targets_dict)

def remove_flanking_outliers(counts_data):
  '''
  filter outlier reads based on their flanking sequence likelihood.
  this is useful for filtering off target reads and reads with large insertions
  in their repeat sequence
  '''
  counts_data["flanking_likelihood"] = counts_data["subset_likelihood"] - counts_data["repeat_likelihood"]
  thresh = counts_data["flanking_likelihood"].mean() - counts_data["flanking_likelihood"].std()
  df_final=counts_data[counts_data.flanking_likelihood > thresh]
  outliers= counts_data[~(counts_data.read_id.isin(df_final.read_id))]
  return df_final, outliers.read_id

def sort_alleles(curr_row, new_row, bootstrap, allele_specif_CIs): #this is only for diploid, need to translate to any allele number
    #establish an order of alleles based on median
    allele_med_cols = curr_row[curr_row.index.str.endswith("median")]
    allele_med_cols_order = allele_med_cols.sort_values(ascending=False) #largest to smallest, can change
    allele_med_cols_order = allele_med_cols.sort_values() 
    for i,allele in enumerate(allele_med_cols_order):
        #get current corresponding original row
        curr_allele = allele_med_cols_order.index[i].split(":")[0]
        curr_allele_cols = curr_row[curr_row.index.str.contains(curr_allele)]
        new_row["A"+str(i+1)+":median"] = curr_allele_cols[curr_allele + ":median"]
        new_row["A"+str(i+1)+":mode"] = curr_allele_cols[curr_allele + ":mode"]
        new_row["A"+str(i+1)+":SD"] = curr_allele_cols[curr_allele + ":SD"]
        new_row["A"+str(i+1)+":supporting_reads"] = curr_allele_cols[curr_allele + ":supporting_reads"]
        if bootstrap:
          new_row["A"+str(i+1)+":median_CI"] = curr_allele_cols[curr_allele + ":median_CI"] 
        if allele_specif_CIs:
          curr_CI = curr_allele_cols[curr_allele + ":median_CI_allele_specific"][0] #tuple, take lower bound and find closest allele assignment
          closest_allele_bootstrap = allele_med_cols_order[(allele_med_cols_order-curr_CI).abs().argsort()].index[0].split(":")[0]
          new_row["A"+ str(i+1)+":median_CI_allele_specific"] = curr_row[closest_allele_bootstrap + ":median_CI_allele_specific"]
    return new_row

def sort_outputs(curr_row, max_peaks,bootstrap,allele_specif_CIs, stranded=False):
  '''
  Use this to sort the column names to be in the same order and order alleles from smallest to largest 
  '''
  column_names_in_order = []
  column_names_in_order.append("name")
  if stranded:
      column_names_in_order.append("strand")
  for i in range(1,max_peaks+1):
      column_names_in_order.append("A"+str(i)+":median")
      column_names_in_order.append("A"+str(i)+":mode")
      column_names_in_order.append("A"+str(i)+":SD")
      column_names_in_order.append("A"+str(i)+":supporting_reads")
      if bootstrap:
        column_names_in_order.append("A"+str(i)+":median_CI")
      if allele_specif_CIs:
        column_names_in_order.append("A"+ str(i)+":median_CI_allele_specific")
  column_names_in_order.append("num_supporting_reads")
  column_names_in_order.append("bandwidth")
  column_names_in_order.append("peak_calling_method")
  new_row = curr_row[column_names_in_order].copy()
  new_row_final = sort_alleles(curr_row, new_row,bootstrap, allele_specif_CIs)
  return new_row_final

def sort_read_assign_alleles(read_assign, geno, max_k):
    '''
    Sorts the alleles in the read_assign file to match the geno file first.
    This function now takes max_k as an argument to dynamically handle the number of alleles.
    '''
    # Generate the list of allele median columns dynamically based on max_k
    allele_columns = [f"A{i}:median" for i in range(1, max_k + 1)]
    
    # Calculate medians for read_assign
    read_assign_allele_medians = read_assign.groupby(["name", "cluster_assignments"]).counts.median().reset_index()
    
    # Melt the geno dataframe to long format for easier comparison, using the dynamically generated allele columns
    geno_allele_medians = geno.melt(id_vars=["name"], value_vars=allele_columns, var_name="allele", value_name="median")
    geno_allele_medians["allele"] = geno_allele_medians["allele"].str.extract('A(\d+):median').astype(int)
    
    # Copy read_assign to avoid modifying the original dataframe
    new_read_assign = read_assign.copy()
    
    # Iterate over unique names
    unique_names = read_assign["name"].unique()
    for name in unique_names:
        geno_allele_medians_subset = geno_allele_medians[geno_allele_medians["name"] == name]
        for allele in read_assign["cluster_assignments"].dropna().unique():
            read_median = read_assign[(read_assign["name"] == name) & (read_assign["cluster_assignments"] == allele)].counts.median()
            closest_allele = geno_allele_medians_subset.iloc[(geno_allele_medians_subset["median"]-read_median).abs().argsort()[:1]]["allele"].to_list()[0]
            new_read_assign["cluster_assignments"] = np.where((read_assign["name"] == name) & (read_assign["cluster_assignments"] == allele), closest_allele, new_read_assign["cluster_assignments"])
    
    return new_read_assign

def get_read_ids(bam_file, chrom,start,end,region_name):
    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Get all reads that align to the specified region
        reads = bam.fetch(chrom, start,end)
        # Get the IDs of the reads
        read_ids = {read.query_name:region_name for read in reads}
    return read_ids
def bam_assign(bam_file,target_coordinates):
    '''
    This function assigns reads to targets based on their alignment to the reference genome
    this mode will be optional and fo now its an additional rather than an alternative input
    will be used to pre-assign reads from global alignment to targets to avoid misassignemnt during WGS runs
    or runs with a large number of targets with non-unique direct flanking sequence
    This mode will also only allow for one-to-one assignment while local alignment only allows multiple
    which is suboptimal for a loarge number of targets in this iteration of the software
    '''
    #for every target region, get all reads that align to it and add them to a dictionary of read_ids to target names
    targets = pd.read_csv(target_coordinates,sep="\t", header=None)
    #check if there is a 5th column in the bed file
    if len(targets.columns) == 5:
        targets.columns = ["chr","start","end","repeat","name"]
    else:
        targets.columns = ["chr","start","end","repeat"]
        targets["name"] = targets.chr + ":" +targets.start.astype(str) + "-" + targets.end.astype(str)
    #for every target in targets, get the reads that align to it
    read_ids = {}
    for target in targets.itertuples():
        read_ids.update(get_read_ids(bam_file, target.chr ,target.start, target.end,target.name))
    return read_ids

