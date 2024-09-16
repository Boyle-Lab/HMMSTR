import pandas as pd
import mappy
import numpy as np
def rev_comp(strand): #can refactor this to use dictionaries
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
  d = {'A': '1', 'T': '2', 'C': '3', 'G': '4'}
  seqInt = ' '.join(d[s] if s in d else s for s in seq.upper())
  #print(seqInt)
  return seqInt
  
def read_fasta(fileObject):
  '''
  Generator function to read in reads
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
  Generator function to read in reads
  '''
  header = ''
  seq = ''
  # skip any useless leading information
  for line in fileObject:
    print(line)
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

def remove_outlier_IQR(df):
    Q1=df.quantile(0.25)
    Q3=df.quantile(0.75)
    IQR=Q3-Q1
    df_final=df[~((df<(Q1-1.5*IQR)) | (df>(Q3+1.5*IQR)))]
    outliers= df[~(df.isin(df_final))]
    return df_final, outliers

def throw_low_cov(df):
    filtered = df[df.freq > 1]
    return filtered

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
      print(input_df)
      input_df.loc[len( pd.DataFrame(input_df).index)] = 0
      input_df.index = ['A','G','C','T','']
    param_dict = input_df.to_dict()
    print(param_dict)
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