import pandas as pd
import mappy
from os.path import exists
from HMMSTR_utils.HMMSTR_utils import seq2int
from subprocess import run, PIPE
from importlib_resources import files

#custom imports
from process_vit_res.process_vit_res import *
class Process_Read:

    def __init__(self, header, seq, cutoff=30, mode = "map-ont", out = ".", k = None, w = None, use_full_seq = False):
        # self.read_id = header.split(" ")[0][1:]
        self.read_id = header #changed so we can just pass the read_id from different sources
        self.seq = seq
        self.cutoff = cutoff
        self.mode = mode
        self.use_full_seq = use_full_seq

        self.k = k
        self.w = w
        #check for if fastas exist?
        if exists(out + "_prefix.fa") and exists(out + "_suffix.fa"):
            self.prefix_fa = out + "_prefix.fa"
            self.suffix_fa = out+"_suffix.fa"
        elif exists(out + "_prefix.fa") == False:
            print("Prefix fasta file does not exist! Please check output path")
        else:
            print("Suffix fasta file does not exist! Please check output path")
        #call alignment and such to set other fields? <-- you can set fields in class methods without initializing them

    def align_mappy(self,target=None):
        '''
        This function aligns the prefixes and suffixes to the current read to assign it to the most likely target.
        '''
        if self.k  is None and self.w is None: # default
            aligner = mappy.Aligner(seq=self.seq,preset=self.mode,best_n = 1) #only return the best alignment per target
        #adjusted k
        elif self.k is not None and self.w is None:
            aligner = mappy.Aligner(seq=self.seq,preset=self.mode,best_n = 1, k =  self.k)
        elif self.k is None and self.w is not None:
            aligner = mappy.Aligner(seq=self.seq,preset=self.mode,best_n = 1, w = self.w)
        else:
            aligner = mappy.Aligner(seq=self.seq,preset=self.mode,best_n = 1, k =self.k , w= self.w) #custom parameters

        #align prefixes
        prefix_dict = {"name":[],"prefix_start":[],"prefix_end":[],"prefix_mapq":[],"strand":[],"alignment_length":[]}
        suffix_dict = {"name":[],"suffix_start":[],"suffix_end":[],"suffix_mapq":[],"strand":[],"alignment_length":[]}

        for name, seq,qual in mappy.fastx_read(self.prefix_fa):
            for hit in aligner.map(seq):
                if hit.mapq < self.cutoff:
                    continue
                #check if we are using bam file assignments, if we are only accept alignments to that target prefix 
                if target is not None:
                    if name != target:
                        continue
                prefix_dict["name"].append(name)
                prefix_dict["prefix_start"].append(hit.r_st)
                prefix_dict["prefix_end"].append(hit.r_en)
                prefix_dict["prefix_mapq"].append(hit.mapq)
                prefix_dict["strand"].append(hit.strand)
                prefix_dict["alignment_length"].append(hit.blen)

        #align_suffixes
        for name, seq,qual in mappy.fastx_read(self.suffix_fa):
            for hit in aligner.map(seq):
                if hit.mapq < self.cutoff:
                    continue
                if target is not None:
                    if name != target:
                        continue
                suffix_dict["name"].append(name)
                suffix_dict["suffix_start"].append(hit.r_st)
                suffix_dict["suffix_end"].append(hit.r_en)
                suffix_dict["suffix_mapq"].append(hit.mapq)
                suffix_dict["strand"].append(hit.strand)
                suffix_dict["alignment_length"].append(hit.blen)

        #record if no valid alignment found
        if len(prefix_dict['name']) < 1:
            self.prefix_df = False
        else:
            self.prefix_df = pd.DataFrame(prefix_dict)
        if len(suffix_dict['name']) < 1:
            self.suffix_df = False
        else:
            self.suffix_df = pd.DataFrame(suffix_dict)
        return
    
    def keep_region(self, prefix_info, suffix_info):
        '''
        This function determines if a read's alignments are eligible prefix and suffix alignments

        Parameters
        --------------------------------------------------------------------------------------------------
        prefix_info: pandas Series. Series contianing alignment info for the prefix
        suffix_info: pandas Series. Series containing alignment info for the suffix
        
        Returns
        --------------------------------------------------------------------------------------------------
        valid: boolean
            True if alignments suggest read is on target, False if alignments in the wrong orientation
        '''
        #check that both contain alignemnts
        if len(suffix_info.index) == 0 or len(prefix_info.index) == 0:
            return False
        #check alignments are on the same strand
        if prefix_info.strand[0] == suffix_info.strand[0]: 
            return True
        else:
            return False
        
    def get_align_info(self, row, prefix_info, suffix_info):
        '''
        This function gets attributes of a given read given that a target has been identified

        Parameters
        ----------------------------------------------------------------------------------------------------
        row: pandas Series. Row for current target.
        prefix_info: pandas Series. Series contianing alignment info for the prefix
        suffix_info: pandas Series. Series containing alignment info for the suffix

        Returns
        ----------------------------------------------------------------------------------------------------
        info: dictionary. Dictionary of alignment and subset information for the current read
        '''
        #dictionary of info for current target
        info = {}
        #may or may not want these in the dictionary, seems like a good idea to keep them together
        info["repeat"] = row.repeat.rstrip()
        info["prefix"] = row.prefix.rstrip()
        info["suffix"] = row.suffix.rstrip() 
        info["prefix_align_length"] = prefix_info.alignment_length[0]
        info["suffix_align_length"] = suffix_info.alignment_length[0]
        #get strand and start and end coordinates
        if prefix_info.strand[0] == 1 and suffix_info.strand[0] == 1:
            info["strand"] = "forward"
            info["align_start"] = prefix_info.prefix_start[0]
            info["align_end"] = suffix_info.suffix_end[0]
            info["end_length"] = info["suffix_align_length"]
            info["start_length"] = info["prefix_align_length"]
        else: #reverse
            info["strand"] = "reverse"
            info["align_start"] = suffix_info.suffix_start[0] #switched back cause I was wrong about which coordinates mappy returns for the reverse case
            info["align_end"] = prefix_info.prefix_end[0]
            info["end_length"] = info["prefix_align_length"]
            info["start_length"] = info["suffix_align_length"]
        #record mapqs
        info["prefix_mapq"] = prefix_info.prefix_mapq[0]
        info["suffix_mapq"] = suffix_info.suffix_mapq[0]

        if self.use_full_seq:
            info["subset"] = self.seq
            #ADDED record subset coords so we can calculate repeat coordinates
            info["subset_start"] = 0
            info["subset_end"] = len(self.seq) -1
            return info

        if info["align_start"] + info["start_length"] - 400 < 0 and info["align_end"]-info["end_length"]+400 < len(self.seq):
            info["subset"] = self.seq[info["align_start"] + info["start_length"] - 50: info["align_end"]-info["end_length"]+400]
            info["subset_start"] = info["align_start"] + info["start_length"] - 50
            info["subset_end"] = info["align_end"]-info["end_length"]+400 -1

        elif info["align_start"] + info["start_length"] - 400 > 0 and info["align_end"]-info["end_length"]+400 > len(self.seq):
            info["subset"] = self.seq[info["align_start"] + info["start_length"] - 400: info["align_end"]-info["end_length"]+50]
            info["subset_start"] = info["align_start"] + info["start_length"] - 400
            info["subset_end"] = info["align_end"]-info["end_length"]+50 -1
        elif info["align_start"] + info["start_length"] - 400 < 0 and info["align_end"]-info["end_length"]+400 > len(self.seq):
            if info["align_start"] + info["start_length"] - 50 < 0 and info["align_end"]-info["end_length"]+50 > len(self.seq):
                info["subset"] = self.seq[info["align_start"] : info["align_end"]]
                info["subset_start"] = info["align_start"]
                info["subset_end"] = info["align_end"]-1
            else:
                info["subset"] = self.seq[info["align_start"] + info["start_length"] - 50: info["align_end"]-info["end_length"]+50] #switch back to 50 if this doesnt help
                info["subset_start"] = info["align_start"] + info["start_length"] - 50
                info["subset_end"] = info["align_end"]-info["end_length"]+50 -1
        else:
            info["subset"] = self.seq[info["align_start"] + info["start_length"] - 400: info["align_end"]-info["end_length"]+400]
            info["subset_start"] = info["align_start"] + info["start_length"] - 400
            info["subset_end"] = info["align_end"]-info["end_length"]+400 -1
        if len(info["subset"]) < 1:
            return #return nothing if there is no repeat in this sequence, spurious alignment
        return info

    def assign_targets(self, targets_df):#TODO this would need to change to accomodate only having either suffix or prefix info, I currently filter reads on this
        '''
        This function takes alignment results for current process read object and determines which target(s)
        are in results

        Parameters
        -----------------------------------------------------------------------------------------------------
        targets_df: pandas DataFrame. DataFrame of tandem repeat target loci to compare to alignment results
        '''
        #check if any alignemnts returned
        if isinstance(self.prefix_df, (bool)) or isinstance(self.suffix_df, (bool)): #no alignments
            self.target_info = {}
            return False #need to decide on final returns for this function still
        #subset to only get targets that aligned according to mappy
        candidate_targets = targets_df[targets_df.name.isin(self.prefix_df.name)] #previously sub_targ
        #get the best alignments per target identified (previously sub_prefixes and sub_suffixes)
        candidate_prefix_aligns = self.prefix_df.groupby('name').head(1).reset_index()
        candidate_suffix_aligns = self.suffix_df.groupby('name').head(1).reset_index()
        self.target_info = {}
        #filter candidates that aren't in a compatible orientation and save final candidates
        for row in candidate_targets.itertuples():
            prefix_info =  candidate_prefix_aligns[candidate_prefix_aligns.name == row.name].reset_index()
            suffix_info = candidate_suffix_aligns[candidate_suffix_aligns.name == row.name].reset_index()
            #save valid regions' attributes
            if self.keep_region(prefix_info, suffix_info):
                self.target_info[row.name] = self.get_align_info(row, prefix_info, suffix_info)

    def run_viterbi(self,hmm_file,rev_hmm_file,hidden_states,rev_states,out,build_pre, prefix_idx,output_labelled_seqs):
        '''
        This function runs viterbi on the current read across all identified targets.

        Parameters
        ----------------------------------------------------------------------------------------------------
        hmm_file: str. Suffix for model file for Viterbi.
        rev_hmm_file: str. Suffix for reverse model file for Viterbi.
        hidden_states: str. Suffix for hidden_states file.
        rev_states: str. Suffix reverse hidden_states file.
        out: str. Output prefix
        build_pre: str. Prefix for all model files
        prefix_idx: int. Offset for prefix states
        output_labelled_seqs: bool. Output labelled read sequence

        Returns
        ----------------------------------------------------------------------------------------------------
        bool for if the run was successful
        '''
        #loop across all identified targets
        #if no targets, return
        if self.target_info == {} or self.target_info is None:
            return False
        for name in self.target_info.keys():
            #choose hmm to use based on strand of target
            curr_hmm_file = build_pre + "_" + name + hmm_file
            curr_rev_file = build_pre + "_" + name + rev_hmm_file
            if self.target_info[name]["strand"] == "forward":
                curr_hmm = curr_hmm_file
                curr_hidden_states_file = open(build_pre + "_" + name + hidden_states,'r')
                curr_states = curr_hidden_states_file.readline().split(".")
                curr_hidden_states_file.close()
            else:
                curr_hmm = curr_rev_file
                curr_hidden_states_rev_file = open(build_pre + "_" + name + rev_states,'r')
                curr_states = curr_hidden_states_rev_file.readline().split(".")
                curr_hidden_states_rev_file.close()
            #convert seqeunce to numeric so it is compatible with the C code
            numeric_seq = str(seq2int(self.target_info[name]["subset"].upper()))
            T = str(len(self.target_info[name]["subset"]))
            repeat_len = len(self.target_info[name]["repeat"])
            repeat_len_str = str(repeat_len)
            prefix_idx_str = str(3*prefix_idx + 1)
            test_hmm_cython_path = files('c_files').joinpath('test_hmm_cython.py')
            command = ['python',test_hmm_cython_path,curr_hmm,T,numeric_seq, repeat_len_str, prefix_idx_str]
            #run viterbi on current read
            result = run(command, universal_newlines=True,capture_output=True, text=True)
            vit_out = result.stdout

            #methods are either static in this class or will be imported by name
            labeled_seq, pointers,MLE = label_states(vit_out, curr_states)
            likelihood, sub_labels,repeats,context, final_repeat_like, repeat_start, repeat_end = calc_likelihood(vit_out, pointers,labeled_seq, curr_states, self.target_info[name]["subset"], self.target_info[name]["subset_start"],self.target_info[name]["subset_end"])
            #save state labels for KMeans method, if time we can figure out how to do this without saving a file
            label_file = open(out+"_"+ name + "_labeled_seqs.txt","a")
            label_file.write(self.read_id + "\t" +".".join(labeled_seq)+"\n")
            label_file.close()
            count = count_repeats(labeled_seq,pointers,repeat_len,self.target_info[name]["subset"])
            score = self.target_info[name]["prefix_mapq"] + self.target_info[name]["suffix_mapq"]

            out_file = open(out+"_"+name+"_counts.txt","a")
            out_file.write(self.read_id + " " + self.target_info[name]["strand"] + " "+ str(score) + " " + str(MLE) + " " + str(likelihood)+ " " + str(final_repeat_like) + " " + str(repeat_start) + " "+ str(repeat_end) + " "+ str(self.target_info[name]["align_start"]) + " "+str(self.target_info[name]["align_end"])+ " " + str(count) + "\n")
            out_file.close()
            #output labelled sequence to context file for given target if output_labelled_seqs is set
            if output_labelled_seqs:
                print_labelled(self.read_id,self.target_info[name]["strand"],sub_labels,context,pointers,out + "_labelled_seqs/"+name+"_context_labeled.txt")
        return True

        

            






        