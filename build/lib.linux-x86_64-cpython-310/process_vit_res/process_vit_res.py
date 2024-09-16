#helper functions to process viterbi outputs
from colorama import Fore, Back, Style
import numpy as np
def label_states(vit_out, hidden_states):
    '''
    Function to convert output from testvit.c to a single sequence of labeled hidden states

    Args:
        vit_out (text file): text file output of testvit.c
        hidden_states_in (text file): text file of "." delimited states in correct order to label with

    Returns:
        labeled_seq (list): a list corresponding to the input sequence in terms of hidden_states
        pointers (dict): a dictionary of "pointers" to the beginning of a state set (P,R, or S) so we can jump there later
        MLE (float): Maximum likelihood estimate of the viterbi path

    '''
    labeled_seq = []
    #pointers to record indeces of the start of each part of our model in the read sequence
    pointers = {"P":False, "R":False, "S":False, "G2":False, "D":[], "I":[]}
    #get likelihoods from viterbi output
    MLE = float(vit_out.split("\n")[0].rstrip().split(" ")[6])
    seq = vit_out.split("\n")[3]
    seq_list = seq.rstrip().split(" ")
    seq_list = [int(i) for i in seq_list]

    #build dictionary to convert from numeric to states
    conversion = {}
    for i,state in enumerate(hidden_states):
        conversion[i+1] = state #viterbi labels are one indexed
    #convert labeled_seq
    for i,obs in enumerate(seq_list):
        labeled_seq.append(conversion[seq_list[i]])
        #check if we need to add a "pointer"
        if  "P" in conversion[seq_list[i]] and pointers["P"] == False:
            pointers["P"] = i
        if "R" in conversion[seq_list[i]] and pointers["R"] == False:
            pointers["R"] = i
        if  "S" in conversion[seq_list[i]] and pointers["S"] == False:
            pointers["S"] = i
        if conversion[seq_list[i]] == "G2" and pointers["G2"] == False:
            pointers["G2"] = i
        if "D" in conversion[seq_list[i]]: #get indeces of all deletions
            pointers["D"].append(i)
        if "I" in conversion[seq_list[i]]: #get indeces of all deletions
            pointers["I"].append(i)
    return labeled_seq, pointers,MLE


def calc_likelihood(vit_out, pointers, labeled_seq,hidden_states, read,subset_start, subset_end):
    '''
    Function to calculate likelihoods corresponding to subseqeunces to be used later

    Args:
        vit_out (text file): text file output of testvit.c
        pointers (dict): dictionary of indeces in labeled_seq corresponding to beginning of different sections of the model
        labeled_seq (str): labeled state seqeunce corresponding to viterbi output
        hidden_states (text file): text file of "." delimited states in correct order to label with
        read (str): raw read sequence
        subset_start (int): integer indicating subset start position in original read
        subset_end (int): integer indeicating subset end position in original read

    Returns:
        final (float): final likelihood corresponding to prefix --> suffix viterbi path
        final_labels(str): seqeunce labels for prefix --> suffix, "-" delimited
        repeats (str): sequence including deletion labels corresponding to only the target repeat in the seqeunce
        context (str): sequence including deletion labels corresponding to the path through prefix-->suffix
        final_repeat_like (float): likelihood of the tandem repeat sequence
        repeat_start (int): repeat start position relative to the original read
        repeat_end (int): repeat end position relative to the original read

    '''
    #we can calculate the likelihood using the labeled sequence
    likelihoods_int = vit_out.split("\n")[4]
    likelihoods_dec = vit_out.split("\n")[5]
    #process likelihood outputs, combine leading integers with decimals
    likelihood_int_list = likelihoods_int.rstrip().split(" ")
    likelihood_dec_list = likelihoods_dec.rstrip().split(" ")
    likelihood_int_list = [int(i) for i in likelihood_int_list]
    likelihood_dec_list = [abs(int(i)) for i in likelihood_dec_list]
    likelihood_list = [float('.'.join(str(i) for i in x)) for x in list(zip(likelihood_int_list,likelihood_dec_list))]
    #get the last G1 and S likelihoods
    Sn = pointers["G2"] - 1
    Gn = pointers["P"] - 1
    final = likelihood_list[Sn] - likelihood_list[Gn]
    final_repeat_like = likelihood_list[pointers["S"]-1] - likelihood_list[pointers["R"]-1] #likelihood of just the repeat

    #get repeat coordinates per read for downstream motif analysis
    #FIXME this doesn't account for the fact that deletions are included in pointers and NOT subset start should be as follows (update install)
    repeat_start = pointers["R"] + subset_start - sum(np.where(np.array(pointers["D"]) < pointers["R"], True, False))
    repeat_end = pointers["S"] + subset_start - sum(np.where(np.array(pointers["D"]) < pointers["S"], True, False))
    # repeat_end =pointers["S"]-1 + subset_start
    # repeat_start = pointers["R"] + subset_start

    final_labels = "-".join(labeled_seq[pointers["P"]:pointers["G2"]])
    num_pd = final_labels.count('PD')
    sub_read = read[pointers["P"]+num_pd:] #need to adjust for deletions but I think it maybe postion dependent

    #get indeces in labelled with deletions
    deletions = [i - pointers["P"] for i in pointers["D"]] #convert to be relative to prefix start position

    #for each deletion, add a "-" in the raw sequence (this may not be the most efficient way to do this)
    context = ""
    offset = 0
    #print(final_labels)
    for i in range(len(labeled_seq)):
        if 'G' in labeled_seq[i]:
            continue
        if 'D' in labeled_seq[i]:
            context= context + "-"
            offset+=1
        else:
            context = context + read[i-offset]

    #now our read is the same length as the labelled seq so we can index with the pointers
    R_start = pointers["R"] - pointers["P"]
    R_end = pointers["S"] - pointers["P"]
    repeats = context[R_start:R_end] #need to see if this includes last index of repeat
    return final, final_labels, repeats, context, final_repeat_like, repeat_start, repeat_end
def count_repeats(labeled_seq, pointers,repeat_len,seq):
    '''
    Function to count the number of repeats in a sequence given a sequence labeled by its hidden states

    Args:
        labeled_seq (list): hidden state annotation of query sequence
        pointers (dict): dictionary of "pointers" to the first match to prefix, repeat, and suffix
        repeat_len (int): length of target repeat
        seq (str): read sequence

    Returns:
        count (int): the number of repeats identified in the query sequence

    '''
    repeat_region = labeled_seq[pointers["R"]:pointers["S"]]
    adjusted_length = len(repeat_region)
    #count insertions to account for them in count
    for state in repeat_region:
        if state[1] == "I":
            adjusted_length -= 1
    #calculate counts from adjusted length, deletions already accounted for in labeled sequence
    count = adjusted_length/repeat_len
    return count

def print_labelled(read_id,strand,sub_labels,context,pointers,out):
    '''
    Function to output color-coded context seqeunce

    Args:
        read_id(str): read id for input sequence
        strand(str): strand of read to output
        sub_labels(str): "-" delimited labels for prefix-->suffix
        context(str): context sequence including deletion labels
        pointers(dict): dictionary of pointers corresponding to indeces of insertions, deletions, and start positions
        out(str): output file suffix
    Returns:
    None, outputs new context sequence string to context file for given target
    '''
    # FIXME there is currently an edge case where if there is a deletion at the end of the repeat the sequence will continue to be labelled in white
    R_start = pointers["R"] - pointers["P"]
    R_end = pointers["S"] - pointers["P"]
    context_list = list(context)
    I = [i - pointers["P"] for i in pointers["I"]]
    D = [i - pointers["P"] for i in pointers["D"]]
    for i in I:
        if i < R_start or i > R_end:
            context_list[i] = '\x1b[5;37;42m' + context_list[i] + Style.RESET_ALL + '\x1b[1;30;40m'
        else:
            context_list[i] = '\x1b[5;37;42m' + context_list[i] +  Style.RESET_ALL + '\x1b[1;37;40m'
    for i in D:
        if i < R_start or i > R_end:
            context_list[i] = '\x1b[5;31;41m'+" " + Style.RESET_ALL + '\x1b[1;30;40m'
        else:
            context_list[i] = '\x1b[5;31;41m' + " "+  Style.RESET_ALL + '\x1b[1;37;40m'
    #replace all inserted bases with color + base to be printed
    file = open(out,"a")
    file.write(read_id + " " + strand + " " + '\x1b[1;30;40m' + "".join(context_list[:R_start]) + Style.RESET_ALL + '\x1b[1;37;40m' + "".join(context_list[R_start:R_end]) + Style.RESET_ALL +'\x1b[1;30;40m'+ "".join(context_list[R_end:])+ Style.RESET_ALL + "\n")
    file.close()
    return
    