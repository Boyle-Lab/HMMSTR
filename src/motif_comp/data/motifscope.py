#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
import argparse
import sys

from Bio import SeqIO
from Bio import Align
from multiprocess import Pool
from itertools import combinations_with_replacement, groupby

import umap
import copy

import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle

from random import random
import scipy.cluster.hierarchy as shc
from scipy.spatial.distance import squareform

import Levenshtein
import binascii
import pylibsais


def parse_fasta(filename):
    """Parse a FASTA file and return a dictionary of sequences"""
    sequences = {}
    for record in SeqIO.parse(filename, 'fasta'):
        if record.seq != "":
            sequences[record.description] = str(record.seq)
    return sequences

def prepare_suffix_string(sequence_dict):
    """Prepare a sequence dictionary in format required by pylibsais."""
    keys = []
    values = []
    index = []
    pos = 0
    for k,v in sequence_dict.items():
        keys.append(k)
        values.append(v)
        pos += len(v) + 1
        index.append(pos)
    index = np.array(index)
    seq = '$'.join(values) + '$'

    return (seq, index) 

def get_positions(suffix_ar, kmer_idx, kmer_cnt):
    """Get all (sorted) positions of a kmer in a sequence"""
    positions = suffix_ar[kmer_idx:kmer_idx+kmer_cnt]
    return np.sort(positions)

def get_kmer_sequence(seq, suffix_ar, kmer_idx, kmer_len):
    """Get the sequence of a kmer"""
    #get location of one of the kmer copies in original sequencedd
    xloc = suffix_ar[kmer_idx]
    #get the kmer sequence
    return seq[xloc:xloc+kmer_len]

def pairwise_alignment(seq1, seq2):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 0
    aligner.mismatch_score = -2
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -5
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0
    alignments = aligner.align(seq1, seq2)
    scores = []
    for alignment in alignments:
        scores.append(alignment.score)
    return max(scores)

def select_best_perm_kmer(new_kmer, kmer_list, previous_kmer):
    perm_kmer = [new_kmer[x:] + new_kmer[:x] for x in range(len(new_kmer))]
    all_kmer_list = [kmer['kmer'] for kmer in kmer_list]
    perm_kmer_exist = []
    for kmer in perm_kmer:
        if kmer in all_kmer_list:
            perm_kmer_exist += [kmer]
    perm_kmer_score = list(map(lambda i: sum([pairwise_alignment(a, i) for a in previous_kmer]), perm_kmer_exist))
    kmer_idx = perm_kmer_score.index(max(perm_kmer_score))
    return perm_kmer_exist[kmer_idx]

def select_best_kmers(k_min, k_max, seq, index, used_kmer, round, min_count=2, min_indiv_count=2, min_consecutive_count=2, min_consecutive_bp=6):
    """Select k-mers based on the amount of sequence masked.

    :param k_min: the minimum k-mer length
    :param k_max: the maximum k-mer length
    :param seq: the sequence to search, as prepared by prepare_suffix_string
    :param index: the index of the end of each sequence in seq, as prepared by prepare_suffix_string
    :param min_count: the minimum number of times a k-mer should occur in the full combined sequence (including overlaps)
    :param min_indiv_count: the minimum number of times a k-mer should occur in a single sequence (excluding overlaps)
    :param min_consecutive_count: the minimum number of consecutive times a k-mer should occur in a single sequence
    :param min_consecutive_bp: the minimum number of consecutive bases that need to be covered by the k-mer
    """

    #create suffix and LCP array
    suffix_ar, lcp_ar = pylibsais.sais(seq)
    #determine maximum length of valid suffixes at each position (should stop at $ and # symbol)
    mkmer_ar = pylibsais.max_suffix(seq)
    

    #get all kmers with min_count or more copies
    #returns no repetitive k-mers, but does return overlapping k-mers
    kmers = list(pylibsais.kmer_count(seq, suffix_ar, lcp_ar, mkmer_ar, k_min, k_max, min_count))
    kmers.sort(key=lambda x: (x[0] * x[2]), reverse=True) #sort on length * count, so that kmers that mask longer sequences are first

    kmers_with_min = []
    res = []
    max_continuous_masked_bp = 0
    evaluated = 0
    #walk across possible kmers
    for kmer_len, kmer_idx, kmer_cnt in kmers:
        #stop if we cannot improve on the current best
        if (kmer_cnt * kmer_len) < max_continuous_masked_bp:
            break

        
        #determine how much of the sequence is masked by this kmer
        total_masked, max_indiv_seq_count, max_consecutive_count = pylibsais.kmer_mask_potential(suffix_ar, mkmer_ar, index, kmer_len, kmer_idx, kmer_cnt)
        evaluated += 1

        #do not report kmer if it is worse than the current best
        #apply filter constraints (see function parameters)
        if max_consecutive_count * kmer_len < max_continuous_masked_bp or \
            max_indiv_seq_count < min_indiv_count or \
            max_consecutive_count < min_consecutive_count or \
            max_consecutive_count * kmer_len < min_consecutive_bp:
            continue

        kmers_with_min += [{'kmer_len': kmer_len, 'kmer_idx': kmer_idx, 'kmer_cnt':kmer_cnt, 'kmer':get_kmer_sequence(seq, suffix_ar, kmer_idx, kmer_len), 'min_kmer': pylibsais.min_string(get_kmer_sequence(seq, suffix_ar, kmer_idx, kmer_len))}]

        max_continuous_masked_bp = max_consecutive_count * kmer_len

        #get the kmer sequence
        kmer_s = get_kmer_sequence(seq, suffix_ar, kmer_idx, kmer_len)
        min_kmer = pylibsais.min_string(kmer_s)
        
        #get all positions of the kmer in 'seq' (can be overlapping)
        positions = get_positions(suffix_ar, kmer_idx, kmer_cnt)

        res.append({'kmer':kmer_s, 'min_kmer': min_kmer, 'suffix_cnt': kmer_cnt, 'total_masked': total_masked, 
                        'max_indiv_seq_count':max_indiv_seq_count, 'max_consecutive_masked':max_consecutive_count * kmer_len, 'pos':positions, 'idx':kmer_idx})
    
    res.sort(key=lambda x: (x['max_consecutive_masked'], x['max_indiv_seq_count'], x['total_masked'], len(x['kmer']), x['kmer']), reverse=True)

    
    print(f'KMER EVALUATED: {evaluated}')
    print(f'KMER SELECTED: {len(res)}')

    if len(res) == 0:
        return ({}, suffix_ar, mkmer_ar, used_kmer)

    else:

        if round == 0:
            #sort kmers on priority: max continuous masked, then max count in individual sequence, then total masked, then length, then alphabetically
            candidate_kmer = res[0]


        else:
            selected_kmer_object = {}
            kmer_set_selected = res[0]['min_kmer']
            min_kmer = res[0]['min_kmer']
            kmer_selected = select_best_perm_kmer(kmer_set_selected, kmers_with_min, used_kmer)
            if seq.index(kmer_selected*2):
                kmer_set = []
                for kmer in kmers_with_min:
                    if kmer['min_kmer'] == kmer_set_selected:
                        kmer_set += [kmer]
                
                for kmer in kmer_set:
                    if kmer['kmer'] == kmer_selected:
                        selected_kmer_object = kmer

                positions = get_positions(suffix_ar, selected_kmer_object['kmer_idx'], selected_kmer_object['kmer_cnt'])
                total_masked, max_indiv_seq_count, max_consecutive_count = pylibsais.kmer_mask_potential(suffix_ar, mkmer_ar, index, selected_kmer_object['kmer_len'], selected_kmer_object['kmer_idx'], selected_kmer_object['kmer_cnt'])
                candidate_kmer = {'kmer':kmer_selected, 'min_kmer': min_kmer, 'suffix_cnt': selected_kmer_object['kmer_cnt'], 'total_masked': total_masked, 
                                'max_indiv_seq_count':max_indiv_seq_count, 'max_consecutive_masked':max_consecutive_count * selected_kmer_object['kmer_len'], 'pos':positions, 'idx':selected_kmer_object['kmer_idx']}
            else:
                candidate_kmer = res[0]
        used_kmer += [candidate_kmer['kmer']]

        return (candidate_kmer, suffix_ar, mkmer_ar, used_kmer)





def select_all_kmer(seq, index, mink, maxk):

    #kmers that are selected
    selected_kmers = []

    #positions that are masked (list of tuples of position and kmer)
    marked_positions = []

    used_kmer = []
    n = 0
    #repeat until no (consequtive) kmers are found
    while True:    
        res, sa, mask, used_kmer = select_best_kmers(mink, maxk, seq, index, used_kmer, n)
        
        if res == {}:
            break
        
        selected = res
        selected_kmers.append(selected) 

        print(f"SELECT KMER: {selected['kmer']}")
        for k,v in selected.items():
            if k != 'kmer':
                print(f"- {k}: {v}")
        
        print('MASKED:')
        rseq, rmarked_pos = pylibsais.kmer_mask(seq, sa, mask, len(selected['kmer']), selected['idx'], selected['suffix_cnt'], 2, '.')
        print(rseq)
        print('\n' * 2)
        if(rseq.count('.') == 0):
            kmer = selected['kmer'] * 2
            #kmer = selected['kmer']
            idx = rseq.index(kmer)
            raise RuntimeError('No masked positions found')
        #mask sequence with # symbol. The '2' indicates that only stretches of at least 2 consecutive kmers are masked.
        seq, marked_pos = pylibsais.kmer_mask(seq, sa, mask, len(selected['kmer']), selected['idx'], selected['suffix_cnt'], 2, '#')
        marked_positions.extend([(e, selected['kmer']) for e in marked_pos])
        n += 1

    
    for selected in selected_kmers:
        seq, marked_pos = pylibsais.kmer_mask_simple(seq, selected['kmer'], '#')
        marked_positions.extend([(e, selected['kmer']) for e in marked_pos])

    return selected_kmers, marked_positions, seq


def select_best_kmers_motif_guided(k_min, k_max, seq, index, used_kmer, round, ref_motifs_dict_r, min_count=2, min_indiv_count=2, min_consecutive_count=2, min_consecutive_bp=6):
    """Select k-mers based on the amount of sequence masked.

    :param k_min: the minimum k-mer length
    :param k_max: the maximum k-mer length
    :param seq: the sequence to search, as prepared by prepare_suffix_string
    :param index: the index of the end of each sequence in seq, as prepared by prepare_suffix_string
    :param min_count: the minimum number of times a k-mer should occur in the full combined sequence (including overlaps)
    :param min_indiv_count: the minimum number of times a k-mer should occur in a single sequence (excluding overlaps)
    :param min_consecutive_count: the minimum number of consecutive times a k-mer should occur in a single sequence
    :param min_consecutive_bp: the minimum number of consecutive bases that need to be covered by the k-mer
    """

    #create suffix and LCP array
    suffix_ar, lcp_ar = pylibsais.sais(seq)
    #determine maximum length of valid suffixes at each position (should stop at $ and # symbol)
    mkmer_ar = pylibsais.max_suffix(seq)
    

    #get all kmers with min_count or more copies
    #returns no repetitive k-mers, but does return overlapping k-mers
    kmers = list(pylibsais.kmer_count(seq, suffix_ar, lcp_ar, mkmer_ar, k_min, k_max, min_count))
    kmers.sort(key=lambda x: (x[0] * x[2]), reverse=True) #sort on length * count, so that kmers that mask longer sequences are first


    kmers_with_min = []
    res = []
    max_continuous_masked_bp = 0
    evaluated = 0
    #walk across possible kmers
    for kmer_len, kmer_idx, kmer_cnt in kmers:
        kmers_with_min += [{'kmer_len': kmer_len, 'kmer_idx': kmer_idx, 'kmer_cnt':kmer_cnt, 'kmer':get_kmer_sequence(seq, suffix_ar, kmer_idx, kmer_len), 'min_kmer': pylibsais.min_string(get_kmer_sequence(seq, suffix_ar, kmer_idx, kmer_len))}]

    for kmer_len, kmer_idx, kmer_cnt in kmers:
        #stop if we cannot improve on the current best
        if (kmer_cnt * kmer_len) < max_continuous_masked_bp:
            break

        
        #determine how much of the sequence is masked by this kmer
        total_masked, max_indiv_seq_count, max_consecutive_count = pylibsais.kmer_mask_potential(suffix_ar, mkmer_ar, index, kmer_len, kmer_idx, kmer_cnt)
        evaluated += 1

        #do not report kmer if it is worse than the current best
        #apply filter constraints (see function parameters)
        if max_consecutive_count * kmer_len < max_continuous_masked_bp or \
            max_indiv_seq_count < min_indiv_count or \
            max_consecutive_count < min_consecutive_count or \
            max_consecutive_count * kmer_len < min_consecutive_bp:
            continue

        max_continuous_masked_bp = max_consecutive_count * kmer_len

        #get the kmer sequence
        kmer_s = get_kmer_sequence(seq, suffix_ar, kmer_idx, kmer_len)
        min_kmer = pylibsais.min_string(kmer_s)
        
        #get all positions of the kmer in 'seq' (can be overlapping)
        positions = get_positions(suffix_ar, kmer_idx, kmer_cnt)

        res.append({'kmer':kmer_s, 'min_kmer': min_kmer, 'suffix_cnt': kmer_cnt, 'total_masked': total_masked, 
                        'max_indiv_seq_count':max_indiv_seq_count, 'max_consecutive_masked':max_consecutive_count * kmer_len, 'pos':positions, 'idx':kmer_idx})
    
    res.sort(key=lambda x: (x['max_consecutive_masked'], x['max_indiv_seq_count'], x['total_masked'], len(x['kmer']), x['kmer']), reverse=True)
    
    print(f'KMER EVALUATED: {evaluated}')
    print(f'KMER SELECTED: {len(res)}')

    if len(res) == 0:
        return ({}, suffix_ar, mkmer_ar, used_kmer)

    else:
        selected_kmer_object = {}
        kmer_set_selected = res[0]['min_kmer']
        min_kmer = res[0]['min_kmer']
        kmer_set = []
        if round == 0:
            if kmer_set_selected in ref_motifs_dict_r:
                for kmer in kmers_with_min: 
                    if kmer["kmer"] == ref_motifs_dict_r[kmer_set_selected]:
                        selected_kmer_object = kmer
                positions = get_positions(suffix_ar, selected_kmer_object['kmer_idx'], selected_kmer_object['kmer_cnt'])
                total_masked, max_indiv_seq_count, max_consecutive_count = pylibsais.kmer_mask_potential(suffix_ar, mkmer_ar, index, selected_kmer_object['kmer_len'], selected_kmer_object['kmer_idx'], selected_kmer_object['kmer_cnt'])
                candidate_kmer = {'kmer':ref_motifs_dict_r[kmer_set_selected], 'min_kmer': min_kmer, 'suffix_cnt': selected_kmer_object['kmer_cnt'], 'total_masked': total_masked, 
                                'max_indiv_seq_count':max_indiv_seq_count, 'max_consecutive_masked':max_consecutive_count * selected_kmer_object['kmer_len'], 'pos':positions, 'idx':selected_kmer_object['kmer_idx']}
                
            else:
                candidate_kmer = res[0]


        else:
            if kmer_set_selected in ref_motifs_dict_r:
                kmer_selected = ref_motifs_dict_r[kmer_set_selected]
                for kmer in kmers_with_min:
                    if kmer['min_kmer'] == kmer_set_selected:
                        kmer_set += [kmer]
                
                for kmer in kmer_set:
                    if kmer['kmer'] == kmer_selected:
                        selected_kmer_object = kmer
                
                positions = get_positions(suffix_ar, selected_kmer_object['kmer_idx'], selected_kmer_object['kmer_cnt'])
                total_masked, max_indiv_seq_count, max_consecutive_count = pylibsais.kmer_mask_potential(suffix_ar, mkmer_ar, index, selected_kmer_object['kmer_len'], selected_kmer_object['kmer_idx'], selected_kmer_object['kmer_cnt'])
                candidate_kmer = {'kmer':kmer_selected, 'min_kmer': min_kmer, 'suffix_cnt': selected_kmer_object['kmer_cnt'], 'total_masked': total_masked, 
                                'max_indiv_seq_count':max_indiv_seq_count, 'max_consecutive_masked':max_consecutive_count * selected_kmer_object['kmer_len'], 'pos':positions, 'idx':selected_kmer_object['kmer_idx']}
            
            else:
                kmer_selected = select_best_perm_kmer(kmer_set_selected, kmers_with_min, used_kmer)

                if kmer_selected*2 in seq:
                    kmer_set = []
                    for kmer in kmers_with_min:
                        if kmer['min_kmer'] == kmer_set_selected:
                            kmer_set += [kmer]
                    
                    for kmer in kmer_set:
                        if kmer['kmer'] == kmer_selected:
                            selected_kmer_object = kmer

                    positions = get_positions(suffix_ar, selected_kmer_object['kmer_idx'], selected_kmer_object['kmer_cnt'])
                    total_masked, max_indiv_seq_count, max_consecutive_count = pylibsais.kmer_mask_potential(suffix_ar, mkmer_ar, index, selected_kmer_object['kmer_len'], selected_kmer_object['kmer_idx'], selected_kmer_object['kmer_cnt'])
                    candidate_kmer = {'kmer':kmer_selected, 'min_kmer': min_kmer, 'suffix_cnt': selected_kmer_object['kmer_cnt'], 'total_masked': total_masked, 
                                    'max_indiv_seq_count':max_indiv_seq_count, 'max_consecutive_masked':max_consecutive_count * selected_kmer_object['kmer_len'], 'pos':positions, 'idx':selected_kmer_object['kmer_idx']}
                else:
                    candidate_kmer = res[0]

        used_kmer += [candidate_kmer['kmer']]

        return (candidate_kmer, suffix_ar, mkmer_ar, used_kmer)


def select_all_kmer_motif_guided(seq, index, mink, maxk, ref_motifs_list):
    ref_motifs_dict = {}
    for motif in ref_motifs_list:
        ref_motifs_dict[motif] = min([motif[x:] + motif[:x] for x in range(len(motif))])
    ref_motifs_dict_r = {value: key for key, value in ref_motifs_dict.items()}

    #kmers that are selected
    selected_kmers = []

    #positions that are masked (list of tuples of position and kmer)
    marked_positions = []

    used_kmer = []
    n = 0
    #repeat until no (consequtive) kmers are found
    while True:    
        res, sa, mask, used_kmer = select_best_kmers_motif_guided(mink, maxk, seq, index, used_kmer, n, ref_motifs_dict_r)
        
        if res == {}:
            break
        
        selected = res
        selected_kmers.append(selected) 

        print(f"SELECT KMER: {selected['kmer']}")
        for k,v in selected.items():
            if k != 'kmer':
                print(f"- {k}: {v}")
        
        print('MASKED:')
        rseq, rmarked_pos = pylibsais.kmer_mask(seq, sa, mask, len(selected['kmer']), selected['idx'], selected['suffix_cnt'], 2, '.')
        print(rseq)
        print('\n' * 2)

        if(rseq.count('.') == 0):
            kmer = selected['kmer']
            if rseq.index(kmer):
                seq, marked_pos = pylibsais.kmer_mask_simple(seq, selected['kmer'], '#')
                marked_positions.extend([(e, selected['kmer']) for e in marked_pos])
                n += 1
                continue
            else:
                raise RuntimeError('No masked positions found')
        #mask sequence with # symbol. The '2' indicates that only stretches of at least 2 consecutive kmers are masked.
        seq, marked_pos = pylibsais.kmer_mask(seq, sa, mask, len(selected['kmer']), selected['idx'], selected['suffix_cnt'], 2, '#')
        marked_positions.extend([(e, selected['kmer']) for e in marked_pos])
        n += 1

    for ref in ref_motifs_list:
        seq, marked_pos = pylibsais.kmer_mask_simple(seq, ref, '#')
        if marked_pos != []:
            print(ref)
            marked_positions.extend([(e, ref) for e in marked_pos])
            if selected_kmers not in selected_kmers:
                selected_kmers += [{'kmer': ref}]


    for selected in selected_kmers:
        #mask sequence with # symbol
        seq, marked_pos = pylibsais.kmer_mask_simple(seq, selected['kmer'], '#')
        marked_positions.extend([(e, selected['kmer']) for e in marked_pos])
           
    return selected_kmers, marked_positions, seq

def mask_all_seq(selected_kmers, marked_positions, seq):
    for s in selected_kmers:
        print(s['kmer'])

    marked_positions.sort(key=lambda x:x[0])
    return marked_positions


def transform_to_df(index, marked_positions, seq_dict):
    seq_df_list = []
    n = 0
    start = 0
    for seq_key, seq_value in seq_dict.items():
        seq_array = pd.DataFrame(list(seq_value), columns=[seq_key])

        pos = [(t[0], t[1]) for t in marked_positions if start <= t[0] < index[n]]

        for p in pos:
            seq_array.iloc[p[0] - start: p[0] - start + len(p[1]), 0] = p[1]

        start = index[n]
        n += 1

        seq_df_list.append(seq_array)

    seq_df = pd.concat(seq_df_list, axis=1)
    return seq_df

def generate_hex_chr_list():
    usable_hex = list(range(1,256))
    usable_hex.remove(60)
    usable_hex.remove(61)
    usable_hex.remove(62)
    usable_hex.remove(45)
    usable_hex.remove(32)
    usable_hex.remove(10)
    usable_hex.remove(13)
    usable_hex = ['{:02x}'.format(x) for x in usable_hex]
    return usable_hex

def assign_chr_to_motif(all_unique_motif, printable_asii):
    if len(all_unique_motif) > len(printable_asii):
        return None, "too many motifs"        
    else:
        printable_asii_to_include = printable_asii[:len(all_unique_motif)]
        motif_dict = dict(zip(all_unique_motif, printable_asii_to_include))
        return motif_dict


def write_seq_in_hex_chr(all_seq_motif_df, motif_dict, cwd_input, random_num):
    all_seq_motif_chr_df = all_seq_motif_df.replace(motif_dict)
    # write into fasta file
    chr_fasta_file = cwd_input + "_motif_in_hex_" + str(random_num) + ".hex"
    seq_chr_dict = {}    
    with open (chr_fasta_file, 'w') as chr_fasta:
        for column in all_seq_motif_df:
            seq_description = column
            seq_motif_in_chr = all_seq_motif_chr_df[column].dropna().tolist()
            seq_motif = all_seq_motif_df[column].dropna().tolist()
            seq_chr = []
            i = 0
            while i < len(seq_motif):
                seq_chr += [seq_motif_in_chr[i]]
                i += len(seq_motif[i])
            seq_chr_string = " ".join(seq_chr)
            chr_fasta.write(">" + seq_description + "\n")
            chr_fasta.write(seq_chr_string + "\n")
            seq_chr_dict[column] = seq_chr_string
            
    motif_fasta_file = cwd_input + "_motif_to_hex_" + str(random_num) + ".fa"
    with open (motif_fasta_file, 'w') as dict_fasta:
        for item in motif_dict:
            dict_fasta.write(">" + motif_dict[item] + "\n")
            dict_fasta.write(item + "\n")        
    return seq_chr_dict
            


def compress_string(string, motif_dict):
    parts = string.split(" ")
    transformed_string = ""
    for key, group in groupby(parts):
        count = len(list(group))
        transformed_string += f"{motif_dict[key]}{count} "
    return transformed_string.strip()


def write_compressed_seq(seq_chr_dict, motif_dict, cwd_input, random_num):
    compressed_fasta_file = cwd_input + "_compressed_representation_" + str(random_num) + ".fa"
    motif_dict_r = {value: key for key, value in motif_dict.items()}
    with open (compressed_fasta_file, 'w') as compressed_fasta_file:
        for seq in seq_chr_dict:
            compressed_fasta_file.write(">" + seq + "\n")
            compressed_fasta_file.write(compress_string(seq_chr_dict[seq], motif_dict_r) + "\n")

def hex_to_bytes(hex_string):
    return binascii.unhexlify(hex_string)


def normalize_edit_distance(distance, seq_1_length, seq_2_length):
    if max(seq_1_length, seq_2_length) != 0:
        return distance / max(seq_1_length, seq_2_length)
    else:
        return 0


def edit_distance_between_seq_byte(input_seq_chr_dict):
    seq_chr_dict = {}
    for seq in input_seq_chr_dict:
        if input_seq_chr_dict[seq] == "":
            seq_chr_dict[seq] = str.encode(" ")
        else:
            seq_chr_dict[seq] = hex_to_bytes(input_seq_chr_dict[seq].replace(" ", ""))
    
    seq_distance_df = pd.DataFrame(list(combinations_with_replacement(list(seq_chr_dict.keys()), 2)), columns=['nseq_1', 'nseq_2'])
    seq_distance_df["seq_1"] = seq_distance_df.apply(lambda x: seq_chr_dict[x["nseq_1"]], axis = 1)
    seq_distance_df["seq_2"] = seq_distance_df.apply(lambda x: seq_chr_dict[x["nseq_2"]], axis = 1)
    ls = pool.starmap(Levenshtein.distance, zip(seq_distance_df.seq_1, seq_distance_df.seq_2))
    seq_distance_df['distance'] = ls
    seq_distance_df['distance'] = seq_distance_df.apply(lambda x: normalize_edit_distance(x["distance"], len(x["seq_1"]), len(x["seq_2"])), axis = 1)

    seq_distance_df = seq_distance_df.drop(columns = ["seq_1", "seq_2"])

    seq_distance_df_2 = seq_distance_df.copy()
    seq_distance_df_2["nseq_2"] = seq_distance_df["nseq_1"]
    seq_distance_df_2["nseq_1"] = seq_distance_df["nseq_2"]
    seq_distance_df = pd.concat([seq_distance_df, seq_distance_df_2], ignore_index = True)
    seq_distance_df = seq_distance_df.drop_duplicates().pivot(index = 'nseq_1', columns = 'nseq_2', values = 'distance').astype(float)

    return seq_distance_df

def get_motif_pairwise_distance(motif_list):
    score_df = pd.DataFrame(list(combinations_with_replacement(motif_list, 2)), columns=['motif_1', 'motif_2'])
    
    score_df['tuple'] = list(zip(score_df.motif_1, score_df.motif_2))
    ls = pool.starmap(pairwise_alignment, score_df['tuple'].tolist())
    score_df['score'] = ls
    score_df = score_df.drop(columns = ['tuple'])

    score_df_2 = score_df.copy()
    score_df_2["motif_2"] = score_df["motif_1"]
    score_df_2["motif_1"] = score_df["motif_2"]
    score_df = pd.concat([score_df, score_df_2], ignore_index = True)
    score_df = score_df.drop_duplicates().pivot(index = 'motif_1', columns = 'motif_2', values = 'score').astype(float)

    return score_df

def edit_distance_between_seq_byte(input_seq_chr_dict):
    seq_chr_dict = {}
    for seq in input_seq_chr_dict:
        if input_seq_chr_dict[seq] == "":
            seq_chr_dict[seq] = str.encode(" ")
        else:
            seq_chr_dict[seq] = hex_to_bytes(input_seq_chr_dict[seq].replace(" ", ""))
    
    seq_distance_df = pd.DataFrame(list(combinations_with_replacement(list(seq_chr_dict.keys()), 2)), columns=['nseq_1', 'nseq_2'])
    seq_distance_df["seq_1"] = seq_distance_df.apply(lambda x: seq_chr_dict[x["nseq_1"]], axis = 1)
    seq_distance_df["seq_2"] = seq_distance_df.apply(lambda x: seq_chr_dict[x["nseq_2"]], axis = 1)
    
    ls = pool.starmap(Levenshtein.distance, zip(seq_distance_df.seq_1, seq_distance_df.seq_2))
    seq_distance_df['distance'] = ls
    seq_distance_df['distance'] = seq_distance_df.apply(lambda x: normalize_edit_distance(x["distance"], len(x["seq_1"]), len(x["seq_2"])), axis = 1)

    seq_distance_df = seq_distance_df.drop(columns = ["seq_1", "seq_2"])

    seq_distance_df_2 = seq_distance_df.copy()
    seq_distance_df_2["nseq_2"] = seq_distance_df["nseq_1"]
    seq_distance_df_2["nseq_1"] = seq_distance_df["nseq_2"]
    seq_distance_df = pd.concat([seq_distance_df, seq_distance_df_2], ignore_index = True)
    seq_distance_df = seq_distance_df.drop_duplicates().pivot(index = 'nseq_1', columns = 'nseq_2', values = 'distance').astype(float)

    return seq_distance_df

def get_motif_pairwise_distance_new(motif_list):
    score_df = pd.DataFrame(list(combinations_with_replacement(motif_list, 2)), columns=['motif_1', 'motif_2'])
    
    score_df['tuple'] = list(zip(score_df.motif_1, score_df.motif_2))
    ls = pool.starmap(pairwise_alignment, score_df['tuple'].tolist())
    score_df['score'] = ls
    score_df = score_df.drop(columns = ['tuple'])
    
    score_df_2 = score_df.copy()
    score_df_2["motif_2"] = score_df["motif_1"]
    score_df_2["motif_1"] = score_df["motif_2"]
    score_df = pd.concat([score_df, score_df_2], ignore_index = True)
    score_df = score_df.drop_duplicates().pivot(index = 'motif_1', columns = 'motif_2', values = 'score').astype(float)

    return score_df

def run_umap(dm, method='UMAP', rank=0.5, norm=True):
    """Maps motif distance matrix 'dm' to a 1-dimensional space using UMAP or MDS.
       :method: string, either 'UMAP' or 'MDS'.
       :rank: float, between 0 and 1. Transforms between original embedding (0.0) and a fully ranked embedding (1.0, i.e. only the ordering of the motifs is preserved).
       :norm: bool, whether to normalize the distance matrix by the geometric mean of the sequence lengths.
    """
    if alignment_score_matrix.shape[0] == 1:
        X_transform_L2 = pd.DataFrame(columns = ["dimension_reduction", "motifs"])
        X_transform_L2["motif"] = dm.columns
        X_transform_L2["dimension_reduction"] = 0
    else:
        n_neighbors = min(10,max(5,dm.shape[0]/2)) # 10 or half of the number of motifs, whichever is smaller, with a minimum of 5
        n_neighbors = int(min(n_neighbors,dm.shape[0])) # cannot be larger than the number of motifs

        data = (-dm).to_numpy(copy=True,dtype=np.float32)

        if norm:
            dnorm = [len(e) for e in dm.index]
            snorm = np.sqrt(dnorm)
            data = data / (snorm[:,np.newaxis] * snorm[np.newaxis,:])
        
        if method == 'UMAP':
            #UMAP has a lot of overhead for small N.
            #This is much faster than builtin method for finding nearest neighbors for small N
            if n_neighbors * 5 < dm.shape[0]: #top-k sort for very large motif sets
                idx = np.argpartition(data,np.arange(n_neighbors),axis=1)[:,:n_neighbors]
                dist = np.take_along_axis(data,idx,axis=1)
            else:
                idx = np.argsort(data,axis=1)[:,:n_neighbors]
                dist = np.take_along_axis(data,idx,axis=1)
            
            #NN-descent only needed for transform of new data (https://github.com/lmcinnes/umap/issues/848)
            import pynndescent
            class DummyNNDescent(pynndescent.NNDescent):
                def __init__(self):
                    return
            precomputed_knn = (idx,dist, DummyNNDescent())
        
            manifold_f = umap.UMAP(n_components = 1, metric = "precomputed", n_neighbors = n_neighbors, min_dist = 0.5, random_state = 0, precomputed_knn=precomputed_knn, force_approximation_algorithm=True)


        elif method=="MDS":
            from sklearn.manifold import MDS
            manifold_f = MDS(n_components=1, n_init=50, metric=False, dissimilarity='precomputed')
        else:
            raise RuntimeError(f'Unknown manifold method {method}. Choose UMAP or MDS.')
    
        result = manifold_f.fit_transform(data) 
        
        if rank:
            result = result.ravel()
            idx = np.argsort(result)
            df =np.diff(result[idx])
            r = np.max(result) - np.min(result)
            w = np.cumsum(df * (1 - rank) + rank * (r / (dm.shape[0] - 1)))
            result[idx[1:]] = w
            result[idx[0]] = 0.0

        X_transform_L2 = pd.DataFrame(result)
        X_transform_L2.columns = ["dimension_reduction"]
        X_transform_L2["motif"] = dm.index
        X_transform_L2 = X_transform_L2.astype({"dimension_reduction": float})

    return X_transform_L2

def map_score_to_alignment(all_seq_motifs, X_transform_L2):
    motif_dict = dict(zip(X_transform_L2.motif, map(float, X_transform_L2.dimension_reduction)))
    motif_dict['nan'] = np.nan
    df = all_seq_motifs.map(lambda x: motif_dict[str(x)])
    return df.transpose()

def summarize_motif(all_seq_motifs, reads_name, input_fasta, random_num, seq_distance_df, figtitle):
    summary = all_seq_motifs.apply(lambda x: x.value_counts(), axis = 0).transpose().reset_index()
    summary["index"] = summary.apply(lambda x: x["index"].split("_")[-1], axis = 1)
    summary = summary.merge(reads_name, left_on = "index", right_on = "nseq", how = "left").drop(columns = ["index", "nseq"]).set_index("read_name").fillna(0)
    motif_file = input_fasta + "_umap_" + str(random_num) + "_only_consecutive_same_string_motif_count" + ".txt"
    summary.to_csv(motif_file, sep = "\t")
    seq_distance_df = seq_distance_df.reset_index().merge(reads_name, left_on = "nseq_1", right_on = "nseq", how = "left").drop(columns = ["nseq_1", "nseq"]).set_index("read_name")
    clusters = shc.linkage(squareform(seq_distance_df), method='average', metric="euclidean")
    
    plt.figure(figsize=(25, 50), dpi = 200)
    d = shc.dendrogram(Z = clusters, labels = seq_distance_df.index, orientation = "right")
    plt.title(figtitle)
    figfile = input_fasta + "_umap_" + str(random_num) + "_only_consecutive_same_string_cluster" + f".{args.format}"
    plt.savefig(figfile, bbox_inches = "tight")
    plt.clf()
    return list(reversed(d["ivl"]))


def plot_df(df, dm, all_seq_motifs, seq_distance_df, figname, figtitle, population_df):
    fig = plt.figure(figsize=(min(max(50, 0.015 * df.shape[1]), 120), min(120,0.5 * df.shape[0])), dpi = 300)
    spec = fig.add_gridspec(ncols=4, nrows=1, width_ratios=[4,1,40,10], height_ratios=[1], wspace=0.02)

    for col in range(3):
        axs = fig.add_subplot(spec[0, col])
    all_axes = fig.get_axes()    

    summary = all_seq_motifs.apply(lambda x: x.value_counts(), axis = 0).transpose().reset_index()
    summary = summary.set_index("index")
    clusters = shc.linkage(squareform(seq_distance_df), method='average', metric="euclidean")

    d = shc.dendrogram(Z = clusters, ax = all_axes[0], labels = seq_distance_df.index, orientation = "left")
    all_axes[0].tick_params(right=False, left = False, top=False, labelright=False, labelleft=False,labeltop=False)
    seq_order_list = list(reversed(d["ivl"]))

    df = df.reindex(seq_order_list)
    all_seq_motifs = all_seq_motifs.transpose()
    all_seq_motifs = all_seq_motifs.reindex(seq_order_list)

    cmap2 = copy.copy(plt.get_cmap('YlGnBu_r'))
    #cmap2 = copy.copy(plt.get_cmap("Spectral"))
    cmap2.set_over('none')

    cbar_ax = fig.add_axes([0.975, .2, .02, .6], title="motif")
    cbaxes = fig.add_axes([0.925, 0.7, .02, 0.1], title="population")
    
    sns.heatmap(df, cmap=cmap2, ax=all_axes[2], cbar_ax = cbar_ax, cbar_kws={"ticks": list(map(float, dm.dimension_reduction))}, yticklabels=True)
   
    all_axes[2].set(ylabel="")
    xticks_positions = np.arange(0, df.shape[1], 50)  # Adjust 100 to your desired interval
    xticks_labels = [str(x) for x in xticks_positions]
    all_axes[2].tick_params(right=True, left = False, top=False, labelright=True, labelleft=False,labeltop=False, rotation=0)
    all_axes[2].set_xticks(xticks_positions)
    all_axes[2].set_xticklabels(xticks_labels, fontsize = 20, rotation=90)
    #all_axes[2].set_xticklabels(xticks_labels, fontsize = 8)
    #all_axes[2].set_yticklabels(all_axes[2].get_yticklabels(), size = 5)
    all_axes[2].tick_params(axis ='x', which ='major')
    cbar_ax.set_yticklabels(dm.motif.tolist()) 
    
    all_seq_motifs = all_seq_motifs.transpose()
    
    for i in range(all_seq_motifs.shape[1]):
        j = 0
        while j < all_seq_motifs.shape[0]:
            if str(all_seq_motifs.iloc[j, i]) != "nan":
                all_axes[2].add_patch(Rectangle((j, i), len(all_seq_motifs.iloc[j, i]), 1, fill=False, edgecolor='black', lw=0.05, clip_on=False))
                j += len(all_seq_motifs.iloc[j, i])
            else:
                j += 1


    seq_population = pd.DataFrame({'seq': seq_order_list})
    seq_population["Sample"] = seq_population.apply(lambda x: x["seq"].split("#")[0], axis = 1)
        
    seq_population = seq_population.merge(population_df[["Sample", "Group"]], on = "Sample", how = "left")
    g = len(seq_population.Group.unique())
    color = ["skyblue", "silver", "lightgreen", "orange", "violet"]
    color_df = pd.DataFrame({'Group': list(seq_population.Group.unique()), 'color': color[:g], 'n': list(range(0, 2*g, 2))})
    seq_population = seq_population.merge(color_df, on = "Group", how = "left")       
    cmap4 = (matplotlib.colors.ListedColormap(color[:g]))
    sns.heatmap(seq_population.set_index("seq")[["n"]], ax=all_axes[1], cmap=cmap4, cbar=False, xticklabels = False, yticklabels = False, linewidths=0.1, linecolor='black')
    all_axes[1].set(ylabel="")

    bounds = list(range(-1, 2*g, 2))
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap4.N)
    cb = plt.colorbar(matplotlib.cm.ScalarMappable(cmap=cmap4, norm=norm), cax=cbaxes, ticks=list(range(0, 2*g, 2)), spacing='uniform', orientation='vertical')
    cbaxes.set_yticklabels(list(seq_population.Group.unique()))
    
    
    cb.outline.set_linewidth(0.1)
    cb.outline.set_edgecolor('black')

    all_axes[2].title.set_text(figtitle)
    plt.yticks(rotation=0)
    plt.savefig(figname, bbox_inches = "tight")


def plot_df_reads(df, dm, all_seq_motifs, seq_distance_df, figname, figtitle):
    fig = plt.figure(figsize=(min(max(50, 0.015 * df.shape[1]), 120), min(120,0.5 * df.shape[0])), dpi = 300)
    spec = fig.add_gridspec(ncols=3, nrows=1, width_ratios=[4,40,5], height_ratios=[1], wspace=0.02)

    for col in range(2):
        axs = fig.add_subplot(spec[0, col])
    all_axes = fig.get_axes()    

    summary = all_seq_motifs.apply(lambda x: x.value_counts(), axis = 0).transpose().reset_index()
    summary = summary.set_index("index")
    clusters = shc.linkage(squareform(seq_distance_df), method='average', metric="euclidean")

    d = shc.dendrogram(Z = clusters, ax = all_axes[0], labels = seq_distance_df.index, orientation = "left")
    all_axes[0].tick_params(right=False, left = False, top=False, labelright=False, labelleft=False,labeltop=False)
    seq_order_list = list(reversed(d["ivl"]))

    df = df.reindex(seq_order_list)
    all_seq_motifs = all_seq_motifs.transpose()
    all_seq_motifs = all_seq_motifs.reindex(seq_order_list)

    cmap2 = copy.copy(plt.get_cmap('YlGnBu_r'))
    #cmap2 = copy.copy(plt.get_cmap("Spectral"))
    cmap2.set_over('none')

    cbar_ax = fig.add_axes([0.975, .2, .02, .6], title="motif")    

    sns.heatmap(df, cmap=cmap2, ax=all_axes[1], cbar_ax = cbar_ax, cbar_kws={"ticks": list(map(float, dm.dimension_reduction))})
    all_axes[1].set(ylabel="")
    xticks_positions = np.arange(0, df.shape[1], 50)  # Adjust to your desired interval
    xticks_labels = [str(x) for x in xticks_positions]
    all_axes[1].tick_params(right=True, left = False, top=False, labelright=True, labelleft=False,labeltop=False,rotation=0)
    all_axes[1].set_xticks(xticks_positions)
    all_axes[1].set_xticklabels(xticks_labels, fontsize = 20,rotation=90)
    all_axes[1].tick_params(axis ='x', which ='major')
    cbar_ax.set_yticklabels(dm.motif.tolist()) 
    
    all_seq_motifs = all_seq_motifs.transpose()
    
    for i in range(all_seq_motifs.shape[1]):
        j = 0
        while j < all_seq_motifs.shape[0]:
            if str(all_seq_motifs.iloc[j, i]) != "nan":
                all_axes[1].add_patch(Rectangle((j, i), len(all_seq_motifs.iloc[j, i]), 1, fill=False, edgecolor='black', lw=0.05, clip_on=False))
                j += len(all_seq_motifs.iloc[j, i])
                
            else:
                j += 1

    all_axes[1].title.set_text(figtitle)
    plt.yticks(rotation=0)
    plt.savefig(figname, bbox_inches = "tight")



def msa_with_characters(cwd_input, random_num):
    chr_fasta_file = cwd_input + "_motif_in_hex_" + str(random_num) + ".hex"
    os.system(f"{mafft_path}/hex2maffttext %s > %s" %(chr_fasta_file, chr_fasta_file.replace(".hex", ".ASCII")))
    os.system(f"mafft --text --op 2.0 --ep 0.1 %s > %s" %(chr_fasta_file.replace(".hex", ".ASCII"), chr_fasta_file.replace(".hex", "_mafft_output.ASCII")))
    os.system(f"{mafft_path}/maffttext2hex %s > %s" % (chr_fasta_file.replace(".hex", "_mafft_output.ASCII"), chr_fasta_file.replace(".hex", "_mafft_output.hex")))
    msa_result = list(SeqIO.parse(chr_fasta_file.replace(".hex", "_mafft_output.hex"), "fasta"))
    os.system("rm %s" %(chr_fasta_file.replace(".hex", ".ASCII")))
    os.system("rm %s" %(chr_fasta_file.replace(".hex", "_mafft_output.ASCII")))
    return msa_result

def prepare_for_plotting(msa_result, motif_dict, dm):
    motif_dict_r = {y: x for x, y in motif_dict.items()}
    motif_dict_r["--"] = ''
    
    msa_df = pd.DataFrame()
    for i in range(len(msa_result)):
        single_seq = str(msa_result[i].seq)
        msa_df = pd.concat([msa_df, pd.DataFrame([single_seq[i:i+2] for i in range(0, len(single_seq), 2)], columns = [msa_result[i].description])], axis = 1)
    
    msa_df["motifs"] = msa_df.apply(lambda x: x.unique().tolist(), axis = 1)
    msa_df["max_length"] = msa_df.apply(lambda x: max(len(motif_dict_r[i]) for i in x["motifs"]), axis = 1)
    msa_df = msa_df.drop(columns = ["motifs"])
    max_motif_length = msa_df["max_length"].tolist()
    msa_adjusted_df = msa_df.loc[msa_df.index.repeat(msa_df.max_length)]
    msa_adjusted_df = msa_adjusted_df.drop(columns = ["max_length"])
    msa_adjusted_df = msa_adjusted_df.transpose()
    msa_adjusted_df = msa_adjusted_df.map(lambda x: motif_dict_r[x])
    msa_adjusted_df.columns = range(msa_adjusted_df.columns.size)

    score_dict = dict(zip(dm.motif, dm.dimension_reduction))
    score_dict[''] = np.nan
    df = msa_adjusted_df.map(lambda x: score_dict[x])
    
    return df, max_motif_length
    

def plot_msa_df(df, dm, all_seq_motifs, seq_distance_df, figname, figtitle, population_df, max_motif_length):
    fig = plt.figure(figsize=(min(max(50, 0.015 * df.shape[1]), 120), min(120,0.5 * df.shape[0])), dpi = 300)
    spec = fig.add_gridspec(ncols=4, nrows=1, width_ratios=[4,1,40,10], height_ratios=[1], wspace=0.01)

    for col in range(3):
        axs = fig.add_subplot(spec[0, col])
    all_axes = fig.get_axes()    

    clusters = shc.linkage(squareform(seq_distance_df), method='average', metric="euclidean")

    d = shc.dendrogram(Z = clusters, ax = all_axes[0], labels = seq_distance_df.index, orientation = "left")
    all_axes[0].tick_params(right=False, left = False, top=False, labelright=False, labelleft=False,labeltop=False)
    seq_order_list = list(reversed(d["ivl"]))


    df = df.reindex(seq_order_list)
    all_seq_motifs = all_seq_motifs.reindex(seq_order_list)

    cmap2 = copy.copy(plt.get_cmap('YlGnBu_r'))
    cmap2.set_over('none')

    
    cbar_ax = fig.add_axes([0.975, .2, .02, .6], title="motif")
    cbaxes = fig.add_axes([0.925, 0.7, .02, 0.1], title="population")

    sns.heatmap(df, cmap=cmap2, ax=all_axes[2], cbar_ax = cbar_ax, cbar_kws={"ticks": list(map(float, dm.dimension_reduction))})
    all_axes[2].set(ylabel="")
    xticks_positions = np.arange(0, df.shape[1], 50)  # Adjust 100 to your desired interval
    xticks_labels = [str(x) for x in xticks_positions]
    all_axes[2].tick_params(right=True, left = False, top=False, labelright=True, labelleft=False,labeltop=False,rotation=0)
    all_axes[2].set_xticks(xticks_positions)
    all_axes[2].set_xticklabels(xticks_labels, fontsize = 20,rotation=90)
    all_axes[2].tick_params(axis ='x', which ='major', rotation=90)
    cbar_ax.set_yticklabels(dm.motif.tolist()) 
    

    for i in range(df.shape[0]):
        j = 0
        m = 0
        while j < df.shape[1]:
            if str(df.iloc[i, j]) != "nan":
                all_axes[2].add_patch(Rectangle((j, i), max_motif_length[m], 1, fill=False, edgecolor='black', lw=0.05, clip_on=False))
            j += max_motif_length[m]
            m += 1

    seq_population = pd.DataFrame({'seq': seq_order_list})
    seq_population["Sample"] = seq_population.apply(lambda x: x["seq"].split("#")[0], axis = 1)
    
    seq_population = seq_population.merge(population_df[["Sample", "Group"]], on = "Sample", how = "left")
    g = len(seq_population.Group.unique())
    color = ["skyblue", "silver", "lightgreen", "orange", "violet"]
    color_df = pd.DataFrame({'Group': list(seq_population.Group.unique()), 'color': color[:g], 'n': list(range(0, 2*g, 2))})
    seq_population = seq_population.merge(color_df, on = "Group", how = "left")       
    cmap4 = (matplotlib.colors.ListedColormap(color[:g]))
    sns.heatmap(seq_population.set_index("seq")[["n"]], ax=all_axes[1], cmap=cmap4, cbar=False, xticklabels = False, yticklabels = False, linewidths=0.1, linecolor='black')
    all_axes[1].set(ylabel="")

    bounds = list(range(-1, 2*g, 2))
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap4.N)
    cb = plt.colorbar(matplotlib.cm.ScalarMappable(cmap=cmap4, norm=norm), cax=cbaxes, ticks=list(range(0, 2*g, 2)), spacing='uniform', orientation='vertical')
    cbaxes.set_yticklabels(list(seq_population.Group.unique()))
    
    cb.outline.set_linewidth(0.1)
    cb.outline.set_edgecolor('black')

    all_axes[2].title.set_text(figtitle)
    plt.yticks(rotation=0)
    plt.savefig(figname, bbox_inches = "tight")


def plot_msa_df_reads(df, dm, all_seq_motifs, seq_distance_df, figname, figtitle, max_motif_length):

    fig = plt.figure(figsize=(min(max(50, 0.015 * df.shape[1]), 120), min(120,0.5 * df.shape[0])), dpi = 300)
    spec = fig.add_gridspec(ncols=3, nrows=1, width_ratios=[4,40,5], height_ratios=[1], wspace=0.01)

    for col in range(2):
        axs = fig.add_subplot(spec[0, col])
    all_axes = fig.get_axes()    

    clusters = shc.linkage(squareform(seq_distance_df), method='average', metric="euclidean")

    d = shc.dendrogram(Z = clusters, ax = all_axes[0], labels = seq_distance_df.index, orientation = "left")
    all_axes[0].tick_params(right=False, left = False, top=False, labelright=False, labelleft=False,labeltop=False)
    seq_order_list = list(reversed(d["ivl"]))

    df = df.reindex(seq_order_list)
    all_seq_motifs = all_seq_motifs.reindex(seq_order_list)

    cmap2 = copy.copy(plt.get_cmap('YlGnBu_r'))
    cmap2.set_over('none')

    cbar_ax = fig.add_axes([0.975, .2, .02, .6], title="motif")

    sns.heatmap(df, cmap=cmap2, ax=all_axes[1], cbar_ax = cbar_ax, cbar_kws={"ticks": list(map(float, dm.dimension_reduction))})
    all_axes[1].set(ylabel="")
    xticks_positions = np.arange(0, df.shape[1], 10) 
    xticks_labels = [str(x) for x in xticks_positions]
    all_axes[1].tick_params(right=True, left = False, top=False, labelright=True, labelleft=False,labeltop=False,rotation=0)
    all_axes[1].set_xticks(xticks_positions)
    all_axes[1].set_xticklabels(xticks_labels, fontsize = 20,rotation=90)
    all_axes[1].tick_params(axis ='x', which ='major', rotation=90)
    cbar_ax.set_yticklabels(dm.motif.tolist()) 
    


    for i in range(df.shape[0]):
        j = 0
        m = 0
        while j < df.shape[1]:
            if str(df.iloc[i, j]) != "nan":
                all_axes[1].add_patch(Rectangle((j, i), max_motif_length[m], 1, fill=False, edgecolor='black', lw=0.05, clip_on=False))
            j += max_motif_length[m]
            m += 1

    all_axes[1].title.set_text(figtitle)
    plt.yticks(rotation=0)
    plt.savefig(figname, bbox_inches = "tight")

def plot_df_single_read(df, dm, all_seq_motifs, figname, figtitle):
    fig = plt.figure(figsize=(min(max(50, 0.015 * df.shape[1]), 120), 1), dpi = 300)
    spec = fig.add_gridspec(ncols=2, nrows=1, width_ratios=[40,5], height_ratios=[1], wspace=0.02)

    for col in range(1):
        axs = fig.add_subplot(spec[0, col])
    all_axes = fig.get_axes()    

    cmap2 = copy.copy(plt.get_cmap('YlGnBu_r'))
    #cmap2 = copy.copy(plt.get_cmap("Spectral"))
    cmap2.set_over('none')

    cbar_ax = fig.add_axes([0.95, .2, .02, .6], title="motif")    

    sns.heatmap(df, cmap=cmap2, ax=all_axes[0], cbar_ax = cbar_ax, cbar_kws={"ticks": list(map(float, dm.dimension_reduction))})
    all_axes[0].set(ylabel="")
    xticks_positions = np.arange(0, df.shape[1], 50)  # Adjust to your desired interval
    xticks_labels = [str(x) for x in xticks_positions]
    all_axes[0].tick_params(right=True, left = False, top=False, labelright=True, labelleft=False,labeltop=False,rotation=0)
    all_axes[0].set_xticks(xticks_positions)
    all_axes[0].set_xticklabels(xticks_labels, fontsize = 20,rotation=90)
    all_axes[0].tick_params(axis ='x', which ='major')
    cbar_ax.set_yticklabels(dm.motif.tolist()) 
    
    i = 0
    j = 0
    while j < all_seq_motifs.shape[0]:
        if str(all_seq_motifs.iloc[j, i]) != "nan":
            all_axes[0].add_patch(Rectangle((j, i), len(all_seq_motifs.iloc[j, i]), 1, fill=False, edgecolor='black', lw=0.05, clip_on=False))
            j += len(all_seq_motifs.iloc[j, i])
            
        else:
            j += 1

    all_axes[0].title.set_text(figtitle)
    plt.yticks(rotation=0)
    plt.savefig(figname, bbox_inches = "tight")


#######################


parser = argparse.ArgumentParser(description='MotifScope')
parser.add_argument('--sequence-type', dest = 'sequence_type', metavar = "[assembly / reads / single]", 
                    help ='Type of input sequences [assembly / reads / single].', type = str,
                    required = True)

parser.add_argument('-i', '--input', default = None, dest='input_fasta_to_count', required = True, 
                    metavar = "input.fa", type = str,
                    help ='Input fasta file to analyze')

parser.add_argument('-mink', '--min_kmer', default = 2, dest = 'min_kmer_size',
                    metavar = 2, type = int,
                    help ='Minimum k to count')

parser.add_argument('-maxk', '--max_kmer', default = 10, dest='max_kmer_size',
                    metavar=10, type=int,
                    help ='Maximum k to count')

parser.add_argument('-t', '--title', default = None, dest='title',
                    metavar = "title", type=str,
                    help ='Title of the plot, default: None')

parser.add_argument('-msa', '--msa', dest = 'run_msa', type = str, metavar = "False",
                    help = 'Boolean (True/False).', default = 'False')

parser.add_argument('-m', '--m', dest = 'motif_guided', type = str, metavar = "False",
                    help = 'Boolean (True/False).', default = 'False')

parser.add_argument('-motif', '--motif', dest = 'ref_motifs', type = str, required = False,
                    help = 'To use established motifs for analysis by taking a file with ref motifs separated by tabs.', metavar = 'motifs.txt')

parser.add_argument('-p', '--population', dest = 'population', required = False, 
                    metavar = "metadata.txt", type=str,
                    help = 'Path to population metadata file')

parser.add_argument('-mp', '--mafft_path', default = None, dest = 'mafft_path', required = False, 
                    metavar = "mafft", type=str,
                    help = 'Path to mafft, default: $PREFIX_CONDA/libexec/mafft')

parser.add_argument('-prof', '--profile', action = 'store_true', dest = 'profile', default = False,
                     help = 'Enable profiling (stored in stats.txt)')

parser.add_argument('-e', '--embed_motif_method', default = 'UMAP', dest= 'embed_motif_method', 
                     help = 'Embedding method for motif color scale (option: MDS or UMAP), default: UMAP')

parser.add_argument('-r', '--motif_rank_embed', default = 0.5, dest = 'motif_rank_embed', type = float,
                     help = 'Hold to original embedding (value=0.0) or only preserve order and place motifs equidistant on color map (value=1.0). Default: 0.5')

parser.add_argument('-f', '--format', default = 'png', dest = 'format',
                    help = 'Image output format (png, pdf, ...). Default: png')

sys.getrecursionlimit()
args = parser.parse_args()

if args.run_msa == "True":
    mafft_path = os.path.join(os.getenv('CONDA_PREFIX'), 'libexec', 'mafft')
    if not os.path.exists(os.path.join(mafft_path, 'hex2maffttext')):
        print(f"Mafft binary hex2maffttext not found within the specified folder '{mafft_path}'. Please provide the correct path to mafft binaries with -mp.")
        sys.exit(1)

if args.mafft_path != None:
    mafft_path = args.mafft_path
title = args.title
sequence_type = args.sequence_type
input_fasta_to_count = args.input_fasta_to_count
max_kmer_size = args.max_kmer_size
min_kmer_size = args.min_kmer_size
run_msa = args.run_msa
motif_guided = args.motif_guided
ref_motifs = args.ref_motifs

if args.profile:
    #imports for profiling
    import cProfile, pstats
    from pstats import SortKey
    print('Profiling enabled. MultiProcessing disabled.')
    pr = cProfile.Profile()
    pr.enable()
    pool_class = DummyPool
else:
    pool_class = Pool

random_num = str(random()).replace('.', '')
cwd_input = os.path.join(os.getcwd(), os.path.basename(input_fasta_to_count))
all_seq_dict = parse_fasta(input_fasta_to_count)
seq_concat, seq_index = prepare_suffix_string(all_seq_dict)

with pool_class() as pool:
    if motif_guided == "False":
        candidate_kmer, masked_postion, masked_seq = select_all_kmer(seq_concat, seq_index, min_kmer_size, max_kmer_size)
        
    elif motif_guided == "True":
        with open(ref_motifs, 'r') as ref:
            ref_motifs = ref.readlines()
        ref_motifs_list = [i.strip().split("\t") for i in ref_motifs]
        ref_motifs_list = [i for l in ref_motifs_list for i in l]
        candidate_kmer, masked_postion, masked_seq = select_all_kmer_motif_guided(seq_concat, seq_index, min_kmer_size, max_kmer_size, ref_motifs_list)

        

    masked_postion = mask_all_seq(candidate_kmer, masked_postion, masked_seq)
    all_seq_df = transform_to_df(seq_index, masked_postion, all_seq_dict)

    unique_motifs = pd.concat([all_seq_df[col] for col in all_seq_df.columns]).unique().tolist()
    if np.nan in unique_motifs:
        unique_motifs.remove(np.nan)

    printable_asii_list = generate_hex_chr_list()
    motif_chr_dict = assign_chr_to_motif(unique_motifs, printable_asii_list)
    all_seq_chr_dict = write_seq_in_hex_chr(all_seq_df, motif_chr_dict, cwd_input, random_num)
    write_compressed_seq(all_seq_chr_dict, motif_chr_dict, cwd_input, random_num)

    all_seq_distance_df = edit_distance_between_seq_byte(all_seq_chr_dict)
    all_seq_chr_dict.clear()
    alignment_score_matrix = get_motif_pairwise_distance(unique_motifs)
    dimension_reduction_result = run_umap(alignment_score_matrix, method=args.embed_motif_method, rank=args.motif_rank_embed)


    if sequence_type == "assembly":
        population = args.population
        population_metadata = pd.read_csv(population, sep = "\t", header = None)
        population_metadata.columns = ["Sample", "Group"]
        if run_msa == "True":
            msa = msa_with_characters(cwd_input, random_num)
            msa_df_to_plot, motif_length = prepare_for_plotting(msa, motif_chr_dict, dimension_reduction_result)
            figfile = cwd_input + "_" + str(random_num) + "_msa" + f".{args.format}"
            plot_msa_df(msa_df_to_plot, dimension_reduction_result, all_seq_df, all_seq_distance_df, figfile, title, population_metadata, motif_length)
        elif run_msa == "False":
            df_to_plot = map_score_to_alignment(all_seq_df, dimension_reduction_result)
            figfile = cwd_input + "_" + str(random_num) + f".{args.format}"
            plot_df(df_to_plot, dimension_reduction_result, all_seq_df, all_seq_distance_df, figfile, title, population_metadata)

    elif sequence_type == "reads":
        if run_msa == "True":
            msa = msa_with_characters(cwd_input, random_num)
            msa_df_to_plot, motif_length = prepare_for_plotting(msa, motif_chr_dict, dimension_reduction_result)
            figfile = cwd_input + "_reads_" + str(random_num) + "_msa" + f".{args.format}"
            plot_msa_df_reads(msa_df_to_plot, dimension_reduction_result, all_seq_df, all_seq_distance_df, figfile, title, motif_length)
        elif run_msa == "False":
            df_to_plot = map_score_to_alignment(all_seq_df, dimension_reduction_result)
            figfile = cwd_input + "_reads_" + str(random_num) + f".{args.format}"
            plot_df_reads(df_to_plot, dimension_reduction_result, all_seq_df, all_seq_distance_df, figfile, title)
    
    elif sequence_type == "single":
        df_to_plot = map_score_to_alignment(all_seq_df, dimension_reduction_result)
        figfile = cwd_input + "_" + str(random_num) + f".{args.format}"
        plot_df_single_read(df_to_plot, dimension_reduction_result, all_seq_df, figfile, title)


if args.profile:
    pr.disable()
    print('Storing profiling results in stats.txt...')
    with open('stats.txt','w') as f:
        sortby = SortKey.CUMULATIVE
        ps = pstats.Stats(pr, stream=f).sort_stats(sortby)
        ps.print_stats()
        ps.print_callers()
        ps.print_callees()

print('DONE')
        




