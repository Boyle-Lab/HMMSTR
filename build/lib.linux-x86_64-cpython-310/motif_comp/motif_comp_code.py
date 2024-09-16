import csv
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
import subprocess
import pandas as pd
import os
import argparse
import gzip
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd
import matplotlib.ticker as ticker
from collections import Counter
import mappy

def load_tandem_repeats_info(tandem_repeats_info_file, target):
    """
    Extract tandem repeat information from the TSV file output by HMMSTR.

    Parameters:
    tandem_repeats_info_file (str): Path to the HMMSTR output TSV file containing tandem repeat information.
    target (str): The target name to filter the rows in the TSV file.

    Returns:
    tandem_repeats_info: A list of dictionaries, each containing information about a tandem repeat.
    set of haplotypes: A set of unique haplotypes found in the tandem repeats.
    """
   
    tandem_repeats_info = []
    haplotypes = set() # Set to store unique haplotypes
    with open(tandem_repeats_info_file, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if row['name'] == target and row['outlier'] != "True":
                haplotype = row['cluster_assignments']
                try:
                    haplotype = haplotype.split(".")[0]
                except ValueError:
                    continue  
                haplotypes.add(haplotype)
                tandem_repeats_info.append({
                    "seq_id": row['read_id'],
                    "align_start": int(row['align_start']),
                    "align_end": int(row['align_end']),
                    "strand": row['strand'],
                    "haplotype": str(haplotype),  
                    "repeat_start": int(row['repeat_start']),
                    "repeat_end": int(row['repeat_end']),
                    # "outlier": row['outlier'] == "True"
                })
    return tandem_repeats_info, list(haplotypes)

def extract_tandem_repeats(read_assignments, sample_file, output_folder, target):
    """
    Extract the tandem repeats sequence for a specified target using information 
    from the read assignments file from HMMSTR.

    Parameters:
    read_assignments (str): Path to the read assignments file.
    sample_file (str): Path to the sample FASTA file.
    output_folder (str): Directory where output files will be saved.
    target (str): The target for which tandem repeats are to be extracted.

    Returns:
    None
    """
    sequence_dict = {name: seq for name, seq, _ in mappy.fastx_read(sample_file)}
        
    print(f"Loaded {len(sequence_dict)} sequences from {sample_file}")

    tandem_repeats_info, haplotypes = load_tandem_repeats_info(read_assignments, target)
    print(f"Processing target {target} with {len(tandem_repeats_info)} tandem repeats")
    
    haplotype_files = {}
    for haplotype in haplotypes:
        haplotype_file = os.path.join(output_folder, f"{target}_H{haplotype}.fa")
        haplotype_files[haplotype] = open(haplotype_file, 'w')

    # Calculate coverage for each haplotype using Counter
    haplotype_counts = Counter(info['haplotype'] for info in tandem_repeats_info if info['haplotype'] in haplotypes)
    # print(haplotype_counts)
    
    for repeat_info in tandem_repeats_info:
        seq_id = repeat_info["seq_id"]
        if seq_id in sequence_dict:
            # if repeat_info["outlier"]:
            #     continue
            align_sequence = sequence_dict[seq_id][repeat_info["align_start"]:repeat_info["align_end"]]

            repeat_start = repeat_info['repeat_start'] - repeat_info["align_start"] 
            repeat_end = repeat_start + (repeat_info["repeat_end"] - repeat_info['repeat_start'])

            if repeat_info["strand"] == "reverse":
                align_sequence = mappy.revcomp(align_sequence)
                seq_length = len(align_sequence)
                repeat_start, repeat_end = seq_length - repeat_end, seq_length - repeat_start

            description = f"repeat_start={repeat_start} repeat_end={repeat_end} coverage={haplotype_counts.get(repeat_info['haplotype'], 0)}"

            new_record = SeqIO.SeqRecord(Seq(align_sequence), id=seq_id, description=description)

            if repeat_info["haplotype"] in haplotypes:
                SeqIO.write(new_record, haplotype_files[repeat_info["haplotype"]], "fasta")
    
    for file in haplotype_files.values():
        file.close()


def update_repeat_positions(msa, repeat_positions):
    """
    Update repeat positions to account for gaps introduced during alignment.

    Parameters:
    msa (MultipleSeqAlignment): The multiple sequence alignment object.
    repeat_positions (dict): A dictionary mapping sequence IDs to their repeat start and end positions.

    Returns:
    dict: A dictionary with updated repeat positions accounting for gaps introduced during alignment.
    """
    updated_positions = {}

    # Create a dictionary to map sequence IDs to their indices in the MSA
    id_to_index = {record.id[:30]: i for i, record in enumerate(msa)}
    
    for seq_id, (repeat_start, repeat_end) in repeat_positions.items():
        # Get the sequence record by ID
        record_index = id_to_index.get(seq_id[:30])     
        if record_index is None:
            continue  
        sequence = str(msa[record_index].seq)
        # Initialize variables for new positions
        new_start = new_end = None
        gap_count = 0
        # Iterate over the sequence to adjust the start position of the repeat
        for i, base in enumerate(sequence):
            if i >= repeat_start + gap_count:
                new_start = i
                break
            if base == '-':
                gap_count += 1
        # Iterate over the sequence to adjust the end position of the repeat
        gap_count = 0
        for i, base in enumerate(sequence):
            if i >= repeat_end + gap_count:
                new_end = i
                break
            if base == '-':
                gap_count += 1
        if new_start is not None and new_end is not None:
            updated_positions[seq_id] = (new_start, new_end)
   
    return updated_positions

def find_matching_key(dictionary, search_string):
    """
    Find a key in the dictionary that matches the search string from characters 1 to 30.
    """
    search_substring = search_string[:30]
    for key in dictionary.keys():
        if key[:30] == search_substring:
            return key
    return None

def concensus_gen(msa, repeat_positions, threshold=0.7):
    """
    Generate a consensus sequence from a multiple sequence alignment (MSA).

    Parameters:
    msa (MultipleSeqAlignment): The multiple sequence alignment object.
    repeat_positions (dict): A dictionary mapping sequence IDs to their repeat start and end positions.
    threshold (float): The threshold for consensus generation.

    Returns:
    list: A list representing the consensus sequence.
    """
    consensus = ""
    # Extract repeat regions from each sequence in the alignment
    repeat_sequences = []
    for i, record in enumerate(msa):
        id = find_matching_key(repeat_positions, record.id)
        repeat_start, repeat_end = repeat_positions[id]
        repeat_sequences.append(record.seq[repeat_start:repeat_end])

    equal_repeat_sequences = []
    max_length = max(len(seq) for seq in repeat_sequences)
    for sequence in repeat_sequences:
        if len(sequence) < max_length:
            sequence += ("-" * (max_length -len(sequence)))
            equal_repeat_sequences.append(sequence)
        else:
            equal_repeat_sequences.append(sequence)
    
    # Build concensus from extracted repeats 
    alignment_length = len(equal_repeat_sequences[0])
    for i in range(alignment_length):
        counts = {}
        for sequence in equal_repeat_sequences:
            base = sequence[i]
            if base in counts:
                counts[base] += 1
            else:
                counts[base] = 1

        consensus_base = max(counts, key=counts.get)
        consensus_percentage = counts[consensus_base] / len(equal_repeat_sequences)
        
        if consensus_base != '-' and consensus_percentage >= threshold:
            consensus += consensus_base
    return consensus

def generate_consensus(read_assignments, sample_file, output_folder, targets="all", clustalw_path='clustalw2', gap_open=8, gap_extension=0.5):
    """
    Align reads using ClustalW and build consensus sequences for each haplotype.

    Parameters:
    read_assignments (str): Path to the TSV file containing read assignments.
    sample_file (str): Path to the sample file (FASTA/FASTQ).
    output_folder (str): Path to the output folder where results will be saved.
    targets (str or list): Target names to process. Default is "all".
    clustalw_path (str): Path to the ClustalW executable. Default is 'clustalw2'.
    gap_open (float): Gap opening penalty for ClustalW. Default is 8.
    gap_extension (float): Gap extension penalty for ClustalW. Default is 0.5.

    Returns:
    None but generates .fa file with all consensus sequences. 
    """
    read_assignments_df = pd.read_csv(f"{read_assignments}", sep="\t")
    # If targets is "all", get all unique target names from the DataFrame
    if targets == "all":
        targets = read_assignments_df['name'].unique()

    for target in targets:        
        # Extract tandem repeats for the current target using the read assignments
        extract_tandem_repeats(read_assignments, sample_file, output_folder, target)

        # Load the unique haplotypes for each target
        _, haplotypes = load_tandem_repeats_info(read_assignments, target)

        consensus_sequences = []

        # Iterate over each haplotype for the current target
        for haplotype in haplotypes:
            target_haplotype_file = os.path.join(output_folder, f"{target}_H{haplotype}.fa")

            if os.path.exists(target_haplotype_file):
                hap_sequences = list(SeqIO.parse(target_haplotype_file, "fasta"))

                # If there is only one sequence, extract the repeat region directly
                if len(hap_sequences) == 1:
                    for rec in hap_sequences:
                        id = rec.description.split()[0]
                        desc_dict = {item.split('=')[0]: int(item.split('=')[1]) for item in rec.description.split() if '=' in item}
                    consensus_hap = hap_sequences[0].seq
                    consensus_hap = consensus_hap[desc_dict['repeat_start']:desc_dict['repeat_end']]
                else:
                    # If there are multiple sequences, align them and extract the consensus
                    repeat_positions = {}
                    for rec in hap_sequences:
                        id = rec.description.split()[0]
                        desc_dict = {item.split('=')[0]: int(item.split('=')[1]) for item in rec.description.split() if '=' in item}
                        repeat_positions[id] = (desc_dict['repeat_start'], desc_dict['repeat_end'])

                    # Run ClustalW to align the sequences         
                    clustalw_command = f"{clustalw_path} -INFILE={target_haplotype_file} -gapopen={gap_open} -gapext={gap_extension}"
                    subprocess.run(clustalw_command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    
                    # Read the alignment file
                    alignment_file = target_haplotype_file.replace(".fa", ".aln")
                    alignment = AlignIO.read(alignment_file, "clustal")
                    msa = MultipleSeqAlignment(alignment)

                    # Update repeat positions based on the alignment
                    updated_positions = update_repeat_positions(msa, repeat_positions)

                    # Generate the consensus sequence
                    consensus_hap = concensus_gen(msa, updated_positions, 0.7)
                    coverage = desc_dict['coverage']

                # Append the consensus sequence to the list
                consensus_sequences.append((f">{target}.H{haplotype}.{coverage}", str(consensus_hap)))
                
        # Define the output file for consensus sequences
        consensus_output = os.path.join(output_folder, "consensus_sequence_file.fa")

        # Write the consensus sequence to the final output file
        with open(consensus_output, 'a') as file:
            for header, sequence in consensus_sequences:
                file.write(f"{header}\n{sequence}\n")

        # Clean up intermediate files
        files_to_delete = [os.path.join(output_folder, f"{target}_H{haplotype}.fa") for haplotype in haplotypes]
        files_to_delete += [file.replace(".fa", ".aln") for file in files_to_delete]
        files_to_delete += [file.replace(".fa", ".dnd") for file in files_to_delete]

        for file_path in files_to_delete:
            if os.path.exists(file_path):
                os.remove(file_path)

def process_and_run_motifscope(motif_script, input_fasta, output_dir, motif_data, targets):
    """
    Process sequences and run the motifscope tool for each target.

    Parameters:
    motif_script (str): Path to the motifscope script.
    input_fasta (str): Path to the input FASTA file containing consensus sequences for each target.
    output_dir (str): Directory where output files will be saved.
    motif_data (str): Path to the TSV file containing motif information for each target.
    targets (str or list): List of target names to process or "all" to process all targets.

    Returns:
    None but saves the motifscope output files.
    """
    # Read the input FASTA file
    sequences = {}
    with open(input_fasta, 'r') as file:
        header = None
        seq = []
        for line in file:
            if line.startswith(">"):
                if header:
                    sequences[header] = ''.join(seq)
                header = line.strip()
                seq = []
            else:
                seq.append(line.strip())
        if header:
            sequences[header] = ''.join(seq)
    
    # Read the motif information from the TSV file
    motif_df = pd.read_csv(motif_data, sep="\t")

    if targets != "all":
        motif_df = motif_df[motif_df['Target'].isin(targets)]

    # Process each target in the list
    for index, row in motif_df.iterrows():
        target_name = row['Target']
        motifs = row['Motifs']
        max_kmer = row['Max_length']
        
        # Create a motif file for the current target
        motif_file = os.path.join(output_dir, f"{target_name}_motifs.txt")
        with open(motif_file, 'w') as file:
            file.write('\n'.join(motifs.split(',')))

        # Filter sequences for the current target
        target_sequences = {}
        for header, seq in sequences.items():
            base_name = header.split('.H')[0][1:]
            if target_name == base_name:
                target_sequences[header] = seq

        # If there are sequences for the current target, create a FASTA file
        if target_sequences:
            target_fasta_file = os.path.join(output_dir, f"{target_name}.fa")
            with open(target_fasta_file, 'w') as fasta_file:
                for header, seq in target_sequences.items():
                    fasta_file.write(f"{header}\n{seq}\n")

            # Dynamically construct the motifscope command based on the number of haplotypes
            num_haplotypes = len(set(header.split('.H')[1] for header in target_sequences.keys()))
            if num_haplotypes > 0:
                motifscope_command = f"""
                python {motif_script} --sequence-type reads -i {target_fasta_file} -mink 2 -maxk {max_kmer} -m True -motif {motif_file}
                """
                # Run motifscope command
                process = subprocess.run(motifscope_command, shell=True, executable='/bin/bash', stdout=subprocess.DEVNULL) 

                # Delete the motif file and temporary fasta file after running motifscope
                os.remove(motif_file)
                os.remove(target_fasta_file)

                print(f"Ran motifscope correctly for target {target_name} with {num_haplotypes} haplotypes.")
            else:
                print(f"No haplotypes found for target {target_name}. Skipping.")
        else:
            print(f"No sequences found for target {target_name}. Skipping.")


def graphing(file_path, targets):
    """
    Generate a graphical representation of the motif composition for the specified targets.

    Parameters:
    file_path (str): Path to the motifscope output file.
    targets (str or list): List of target names to process or "all" to process all targets.

    Returns:
    None but saves the plot as a PNG file.
    """
    # Parsing the motifscope output file
    haplotypes = {}
    with open(file_path, 'r') as file:
        current_haplotype = None
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                current_haplotype = line[1:]
                haplotypes[current_haplotype] = []
            else:
                motifs = line.split()
                for motif in motifs:
                    seq = ''.join(filter(str.isalpha, motif))
                    count = int(''.join(filter(str.isdigit, motif)))
                    haplotypes[current_haplotype].extend([seq] * count)
    
    # Filtering for specified targets
    if targets == "all":
        filtered_motifs = haplotypes
    else:
        filtered_motifs = {}
        for target in targets:  
            target_motifs = {k: v for k, v in haplotypes.items() if k.startswith(f"{target}.H")}
            filtered_motifs.update(target_motifs)
    
    # Sorting haplotypes to ensure H1 is on top of H2 (reverse order)
    sorted_haplotypes = dict(sorted(filtered_motifs.items(), key=lambda x: x[0], reverse=True))

    # Generating colors for motifs
    all_motifs = [motif for motifs in sorted_haplotypes.values() for motif in motifs]
    unique_motifs = list(set(all_motifs))
    palette = sns.color_palette("colorblind", len(unique_motifs))
    palette = [color for color in palette if color != 'grey']  # Removing grey as it is used to indicate non-repeat sequence
    colors = {motif: color for motif, color in zip(unique_motifs, palette)}

    # Creating the plot
    bar_height = 0.4  
    if len(sorted_haplotypes) == 1:
        height = (3+len(sorted_haplotypes))/len(sorted_haplotypes)
        fig, ax = plt.subplots(figsize=(16, 2.5))  
    else:
        fig, ax = plt.subplots(figsize=(16, 3 + len(sorted_haplotypes)))

    # Adjust genome length display for sequences longer than 1000
    max_length = 0
    for haplotype, motifs in sorted_haplotypes.items():
        max_length += sum(len(motif) for motif in motifs)

    if max_length > 1000:
        genome_length = 40
    else:
        genome_length = 3
    
    max_length = 0  # Track the maximum length of the haplotype bars for setting x-axis limits
    motif_start_offset = int(genome_length + bar_height / 2)

    for i, (haplotype, motifs) in enumerate(sorted_haplotypes.items()):
        start = 0

        # Draw left semicircle
        left_semi = mpatches.Wedge(
            (start, i), bar_height / 2, 90, 270,
            linewidth=1, edgecolor='grey', facecolor='grey'
        )
        ax.add_patch(left_semi)
        
        # Draw left rectangle 
        left_rect = mpatches.Rectangle(
            (start, i - bar_height / 2), genome_length, bar_height,
            linewidth=1, edgecolor='grey', facecolor='grey'
        )
        ax.add_patch(left_rect)
        start += genome_length

        motif_length = sum(len(motif) for motif in motifs)

        for motif in motifs:
            color = colors[motif]
            ax.barh(i, len(motif), height=bar_height, left=start, color=color, edgecolor='grey', linewidth=0.5)
            start += len(motif)
        
        # Draw right rectangle 
        right_rect = mpatches.Rectangle(
            (start, i - bar_height / 2), genome_length, bar_height,
            linewidth=1, edgecolor='grey', facecolor='grey'
        )
        ax.add_patch(right_rect)
        start += genome_length
        
        # Draw right semicircle
        right_semi = mpatches.Wedge(
            (start, i), bar_height / 2, 270, 90,
            linewidth=1, edgecolor='grey', facecolor='grey'
        )
        ax.add_patch(right_semi)

        # Add shine effect across the entire bar including genome sections
        shine_length = motif_length + genome_length * 2
        shine_rect = mpatches.Rectangle(
            (0, i + bar_height / 40), shine_length, bar_height / 3.5,
            linewidth=0, edgecolor='none', facecolor='white', alpha=0.2
        )
        ax.add_patch(shine_rect)

        # Update max_length for x-axis scaling
        if start > max_length:
            max_length = start

    # Setting plot labels and grid
    ax.set_yticks(range(len(sorted_haplotypes)))
    if len(sorted_haplotypes) == 1:
        labels = [f"coverage: {hap.split('.')[2]}  H1" for hap in sorted_haplotypes.keys()]
    else:
        labels = [f"coverage: {hap.split('.')[2]}  {hap.split('.')[1]}" for hap in sorted_haplotypes.keys()]

    ax.set_yticklabels(labels)
    ax.set_ylim(-0.5, len(sorted_haplotypes) - 0.5)

    # Adjust x-axis to start from -50 for sequences longer than 1000
    if max_length + motif_start_offset > 1000:
        ax.set_xlim(left=-50, right=max_length + motif_start_offset)
    else:
        ax.set_xlim(left=-motif_start_offset, right=max_length + motif_start_offset)
    
    # Adjust x-axis ticks and labels
    if max_length + motif_start_offset <= 200:
        # Major ticks every 10
        x_ticks = range(motif_start_offset, max_length + motif_start_offset + 10, 10)
        ax.set_xticks(x_ticks)
        ax.set_xticklabels([tick - motif_start_offset for tick in x_ticks])
        
        # Minor ticks every 1
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
        ax.tick_params(axis='x', which='minor', length=4)
    elif max_length + motif_start_offset < 1000:
        # Major ticks every 50
        x_major_ticks = range(motif_start_offset, max_length + motif_start_offset + 50, 50)
        ax.set_xticks(x_major_ticks)
        ax.set_xticklabels([tick - motif_start_offset for tick in x_major_ticks])

        # Minor ticks every 25
        x_minor_ticks = range(motif_start_offset, max_length + motif_start_offset + 25, 25)
        ax.set_xticks(x_minor_ticks, minor=True)
    else:
        # Major ticks every 200
        x_major_ticks = range(motif_start_offset, max_length + motif_start_offset + 200, 200)
        ax.set_xticks(x_major_ticks)
        ax.set_xticklabels([tick - motif_start_offset for tick in x_major_ticks])

        # Minor ticks every 50
        x_minor_ticks = range(motif_start_offset, max_length + motif_start_offset + 50, 50)
        ax.set_xticks(x_minor_ticks, minor=True)

    # Set title 
    if len(targets) == 1:
        title = f"{targets[0]} Motif Composition"
    else:
        title = "Motif Composition"
    ax.set_title(title, fontsize=11)

    # Creating the legend
    patches_legend = [mpatches.Patch(color=color, label=label) for label, color in colors.items()]
    ax.legend(handles=patches_legend, title='Motifs', bbox_to_anchor=(0.5, -0.5), loc='upper center', ncol=5)
    plt.subplots_adjust(bottom=0.6)

    # Save the figure
    if len(targets) == 1:
        plt.savefig(f"{targets[0]}_motif_comp.png")
    else:
        plt.savefig("motif_comp.png")
    
    # plt.show()