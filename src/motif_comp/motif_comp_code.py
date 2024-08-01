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

def load_tandem_repeats_info(tandem_repeats_info_file, target):
    """
    Extract tandem repeat information from the TSV file output by HMMSTR.
    """
    tandem_repeats_info = []
    with open(tandem_repeats_info_file, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if row['name'] == target:
                haplotype = row['cluster_assignments']
                if haplotype in ['1', '1.0']:
                    haplotype = '1'
                elif haplotype in ['2', '2.0']:
                    haplotype = '2'
                tandem_repeats_info.append({
                    "seq_id": row['read_id'],
                    "align_start": int(row['align_start']),
                    "align_end": int(row['align_end']),
                    "strand": row['strand'],
                    "haplotype": haplotype,
                    "repeat_start": int(row['repeat_start']),
                    "repeat_end": int(row['repeat_end']),
                    "outlier": row['outlier'] == "True"
                })
    return tandem_repeats_info

def extract_tandem_repeats(read_assignments, sample_file, output_folder, target):
    """
    Extract the tandem repeats sequence for a specified target using information 
    from the read assignments file from HMMSTR.
    """
    with gzip.open(sample_file, "rt") as handle:
        sequence_dict = SeqIO.to_dict(SeqIO.parse(handle, "fastq"))
    print(f"Loaded {len(sequence_dict)} sequences from {sample_file}")

    output_fasta_file_1 = os.path.join(output_folder, f"{target}_H1.fa")
    output_fasta_file_2 = os.path.join(output_folder, f"{target}_H2.fa")

    with open(output_fasta_file_1, 'w') as output_1:
        tandem_repeats_info = load_tandem_repeats_info(read_assignments, target)
        
        print(f"Processing target {target} with {len(tandem_repeats_info)} tandem repeats")

        # Calculate coverage for each haplotype using Counter
        haplotype_counts = Counter(info['haplotype'] for info in tandem_repeats_info if info['haplotype'] in ['1', '2'])
        coverage_H1 = haplotype_counts.get('1', 0)
        coverage_H2 = haplotype_counts.get('2', 0)
        
        print(f"Coverage for Haplotype 1: {coverage_H1}")
        print(f"Coverage for Haplotype 2: {coverage_H2}")

        hap2_exists = coverage_H2 > 0
        print(f"Haplotype 2 exists: {hap2_exists}")

        if hap2_exists:
            output_2 = open(output_fasta_file_2, 'w')
        
        for repeat_info in tandem_repeats_info:
            seq_id = repeat_info["seq_id"]
            if seq_id in sequence_dict:
                if repeat_info["outlier"]:
                    continue
                align_sequence = sequence_dict[seq_id].seq[repeat_info["align_start"]:repeat_info["align_end"]]

                repeat_start = repeat_info['repeat_start']-repeat_info["align_start"] 
                repeat_end = repeat_start + (repeat_info["repeat_end"]-repeat_info['repeat_start'])

                if repeat_info["strand"] == "reverse":
                    align_sequence = align_sequence.reverse_complement()
                    seq_length = len(align_sequence)
                    repeat_start, repeat_end = seq_length - repeat_end, seq_length - repeat_start

                description = f"repeat_start={repeat_start} repeat_end={repeat_end} coverage={coverage_H1 if repeat_info['haplotype'] == '1' else coverage_H2}"

                new_record = SeqIO.SeqRecord(align_sequence, id=seq_id, description=description)

                if repeat_info["haplotype"] == "":
                    continue
                if repeat_info["haplotype"] in ["1", "1.0"]:
                    SeqIO.write(new_record, output_1, "fasta")
                elif hap2_exists and repeat_info["haplotype"] in ["2", "2.0"]:
                    SeqIO.write(new_record, output_2, "fasta")
        
        if hap2_exists:
            output_2.close()

def update_repeat_positions(msa, repeat_positions):
    """
    Update repeat positions to account for gaps introduced during alignment.
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
    Calculate the consensus sequence considering only the repeat region.
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
    """
    read_assignments_df = pd.read_csv(f"{read_assignments}", sep="\t")
    if targets == "all":
        targets = read_assignments_df['name'].unique()

    for target in targets:
        extract_tandem_repeats(read_assignments, sample_file, output_folder, target)
        target_H1 = os.path.join(output_folder, f"{target}_H1.fa")
        target_H2 = os.path.join(output_folder, f"{target}_H2.fa")
        
        consensus_sequences = []

        if os.path.exists(target_H1):
            hap1_sequences = list(SeqIO.parse(target_H1, "fasta"))
            if len(hap1_sequences) == 1:
                for rec in hap1_sequences:
                    id = rec.description.split()[0]
                    desc_dict = {item.split('=')[0]: int(item.split('=')[1]) for item in rec.description.split() if '=' in item}
                consensus_H1 = hap1_sequences[0].seq
                consensus_H1 = consensus_H1[desc_dict['repeat_start']:desc_dict['repeat_end']]
            else:
                repeat_positions = {}
                for rec in hap1_sequences:
                    id = rec.description.split()[0]
                    desc_dict = {item.split('=')[0]: int(item.split('=')[1]) for item in rec.description.split() if '=' in item}
                    repeat_positions[id] = (desc_dict['repeat_start'], desc_dict['repeat_end'])
                                             
                clustalw_command_1 = f"{clustalw_path} -INFILE={target_H1} -gapopen={gap_open} -gapext={gap_extension}"
                subprocess.run(clustalw_command_1, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                alignment_file_1 = target_H1.replace(".fa", ".aln")
                alignment_1 = AlignIO.read(alignment_file_1, "clustal")
                msa_1 = MultipleSeqAlignment(alignment_1)
                updated_positions_1 = update_repeat_positions(msa_1, repeat_positions)
                consensus_H1 = concensus_gen(msa_1, updated_positions_1, 0.7)
                coverage = desc_dict['coverage']
            consensus_sequences.append((f">{target}.H1.{coverage}", str(consensus_H1)))

        if os.path.exists(target_H2):
            hap2_sequences = list(SeqIO.parse(target_H2, "fasta"))
            if len(hap2_sequences) == 1:
                for rec in hap1_sequences:
                    id = rec.description.split()[0]
                    desc_dict = {item.split('=')[0]: int(item.split('=')[1]) for item in rec.description.split() if '=' in item}
                consensus_H2 = hap1_sequences[0].seq
                consensus_H2 = consensus_H2[desc_dict['repeat_start']:desc_dict['repeat_end']]
            else:
                repeat_positions = {}
                for rec in hap2_sequences:
                    id = rec.description.split()[0]
                    desc_dict = {item.split('=')[0]: int(item.split('=')[1]) for item in rec.description.split() if '=' in item}
                    repeat_positions[id] = (desc_dict['repeat_start'], desc_dict['repeat_end'])

                clustalw_command_2 = f"{clustalw_path} -INFILE={target_H2} -gapopen={gap_open} -gapext={gap_extension}"
                subprocess.run(clustalw_command_2, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                alignment_file_2 = target_H2.replace(".fa", ".aln")
                alignment_2 = AlignIO.read(alignment_file_2, "clustal")
                msa_2 = MultipleSeqAlignment(alignment_2)
                updated_positions_2 = update_repeat_positions(msa_2, repeat_positions)
                consensus_H2 = concensus_gen(msa_2, updated_positions_2, 0.7)
                coverage = desc_dict['coverage']
                consensus_sequences.append((f">{target}.H2.{coverage}", str(consensus_H2)))

        consensus_output = os.path.join(output_folder, "consensus_sequence_file.fa")
        with open(consensus_output, 'a') as file:
            for header, sequence in consensus_sequences:
                file.write(f"{header}\n{sequence}\n")

        # Clean up intermediate files
        files_to_delete = [
            target_H1,
            target_H2,
            target_H1.replace(".fa", ".aln"),
            target_H2.replace(".fa", ".aln"),
            target_H1.replace(".fa", ".dnd"),
            target_H2.replace(".fa", ".dnd")]

        for file_path in files_to_delete:
            if os.path.exists(file_path):
                os.remove(file_path)

def process_and_run_motifscope(input_fasta, output_dir, motif_data, targets):
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
    
    # Process each target in the DataFrame
    motif_df = pd.read_csv(motif_data, sep="\t")

    if targets != "all":
        motif_df = motif_df[motif_df['name'].isin(targets)]
    
    for index, row in motif_df.iterrows():
        target_name = row['name']
        motif = row['repeat']
        
        # Create a motif file for the current target
        motif_file = os.path.join(output_dir, f"{target_name}_motifs.txt")
        with open(motif_file, 'w') as file:
            file.write(f"{motif}")

        # Filter sequences for the current target
        target_sequences = {}
        for header, seq in sequences.items():
            # Extract the base target name from the header
            base_name = header.split('.H')[0][1:]
            if target_name == base_name:
                target_sequences[header] = seq

        # If there are sequences for the current target, create a FASTA file
        if target_sequences:
            target_fasta_file = os.path.join(output_dir, f"{target_name}.fa")
            with open(target_fasta_file, 'w') as fasta_file:
                for header, seq in target_sequences.items():
                    fasta_file.write(f"{header}\n{seq}\n")

            # Determine the command based on the number of haplotypes
            if len(target_sequences) == 1:
                env_name = "motifscope"
                # Create the motifscope command
                motifscope_command = f"""
                source /home/imf/miniconda3/etc/profile.d/conda.sh && \
                conda activate {env_name} && \
                motifscope --sequence-type reads -i {target_fasta_file} -mink 2 -maxk 20 -m True -motif {motif_file}
                """

            elif len(target_sequences) == 2:
                env_name = "motifscope"
                # Create the motifscope command
                motifscope_command = f"""
                source /home/imf/miniconda3/etc/profile.d/conda.sh && \
                conda activate {env_name} && \
                motifscope --sequence-type reads -i {target_fasta_file} -mink 2 -maxk 20 -m True -motif {motif_file}
                """

            # Run motifscope command
            process = subprocess.run(motifscope_command, shell=True, executable='/bin/bash') 
               
            # Delete the motif file and temporary fasta file after running motifscope
            os.remove(motif_file)
            os.remove(target_fasta_file)

            deactivate_command = """
            source /home/imf/miniconda3/etc/profile.d/conda.sh && \
            conda deactivate
            """
            subprocess.run(deactivate_command, shell=True, executable='/bin/bash')

def graphing(file_path, targets):
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
    colors = {motif: color for motif, color in zip(unique_motifs, palette)}

    # Creating the plot
    bar_height = 0.4  # Set the height of the bars
    fig, ax = plt.subplots(figsize=(16, 3)) #0.4 * len(sorted_haplotypes)))  

    max_length = 0  # Track the maximum length of the haplotype bars for setting x-axis limits
    motif_start_offset = int(3 + bar_height / 2)  # Convert to integer

    for i, (haplotype, motifs) in enumerate(sorted_haplotypes.items()):
        start = 0

        # Draw left semicircle
        left_semi = mpatches.Wedge(
            (start, i), bar_height / 2, 90, 270, 
            linewidth=1, edgecolor='grey', facecolor='grey'
        )
        ax.add_patch(left_semi)
        
        # Draw left rectangle (3 units)
        left_rect = mpatches.Rectangle(
            (start, i - bar_height / 2), 3, bar_height,
            linewidth=1, edgecolor='grey', facecolor='grey'
        )
        ax.add_patch(left_rect)
        start += 3

        motif_length = sum(len(motif) for motif in motifs)

        for motif in motifs:
            color = colors[motif]
            ax.barh(i, len(motif), height=bar_height, left=start, color=color, edgecolor='grey', linewidth=0.5) 
            start += len(motif)
        
        # Draw right rectangle (3 units)
        right_rect = mpatches.Rectangle(
            (start, i - bar_height / 2), 3, bar_height,
            linewidth=1, edgecolor='grey', facecolor='grey'
        )
        ax.add_patch(right_rect)
        start += 3
        
        # Draw right semicircle
        right_semi = mpatches.Wedge(
            (start, i), bar_height / 2, 270, 90, 
            linewidth=1, edgecolor='grey', facecolor='grey'
        )
        ax.add_patch(right_semi)

        # Add shine effect across the entire bar including genome sections
        shine_length = motif_length + 6  
        shine_rect = mpatches.Rectangle(
            (0, i + bar_height / 40), shine_length, bar_height / 3.5,
            linewidth=0, edgecolor='none', facecolor='white', alpha=0.2
        )
        ax.add_patch(shine_rect)

        # Update max_length for x-axis scaling
        if start > max_length:
            max_length = start

    # Setting plot labels and grid
    if len(sorted_haplotypes) == 1:
        ax.set_yticks([0])
        labels = [f"coverage: {hap.split('.')[2]}  {hap.split('.')[1]}" for hap in sorted_haplotypes.keys()]
        ax.set_yticklabels(labels)
        # ax.set_yticklabels(f"coverage: {hap.split('.')[2]}" for hap in sorted_haplotypes.keys() [list(sorted_haplotypes.keys())[0].split(".")[1]])
        ax.set_ylim(-1, 1)
    else:
        ax.set_yticks(range(len(sorted_haplotypes)))
        labels = [f"coverage: {hap.split('.')[2]}  {hap.split('.')[1]}" for hap in sorted_haplotypes.keys()]
        ax.set_yticklabels(labels)
        ax.set_ylim(-0.5, 1.5)
    
    
    # Adjust x-axis to start from 0 at the start of the repeat region
    ax.set_xlim(left=-motif_start_offset, right=max_length + motif_start_offset)
    
    # Adjust x-axis ticks and labels
    if max_length + motif_start_offset <= 200:
        x_ticks = range(motif_start_offset, max_length + motif_start_offset + 10, 10)
        ax.set_xticks(x_ticks)
        ax.set_xticklabels([tick - motif_start_offset for tick in x_ticks])
        # Add minor ticks without labels
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
        ax.tick_params(axis='x', which='minor', length=4)
    elif max_length + motif_start_offset < 1000:
        # Major ticks every 50
        x_major_ticks = range(motif_start_offset, max_length + motif_start_offset + 50, 50)
        ax.set_xticks(x_major_ticks)
        ax.set_xticklabels([tick - motif_start_offset for tick in x_major_ticks])

        # Minor ticks every 25
        x_minor_ticks = range(motif_start_offset, max_length + motif_start_offset + 25, 10)
        ax.set_xticks(x_minor_ticks, minor=True)
    else:
        x_major_ticks = range(motif_start_offset, max_length + motif_start_offset + 200, 200)
        ax.set_xticks(x_major_ticks)
        ax.set_xticklabels([tick - motif_start_offset for tick in x_major_ticks])

        # Minor ticks every 50
        x_minor_ticks = range(motif_start_offset, max_length + motif_start_offset + 25, 25)
        ax.set_xticks(x_minor_ticks, minor=True)

    # Add minor ticks without labels
    # ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    # ax.tick_params(axis='x', which='minor', length=4)

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
    #plt.show()