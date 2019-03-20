import generic as gen
import seq_ops as seqo
import sequence_ops as sequo
import sim_ops_containers as soc
import main_tests_ops as mto
import ops
import numpy as np
import collections
import re
import os
import pandas as pd
# import scipy.stats as stats
# import matplotlib.pyplot as plt
import conservation
import os
import zipfile
from useful_motif_sets import stops, codon_map
from regex_patterns import codon_pattern
import containers as cont
import itertools as it
from Bio import SeqIO
import scipy.stats




# motifs = "source_data/motif_sets/int3.txt"
# motifs = "source_data/motif_sets/RESCUE.txt"
# seqs_file = "clean_run/lincrna/lincRNA.multi_exon.exons.fasta"
# seqs_file = "clean_run/genome_sequences/lincrna/cabili/multi_exons.fasta"
#

# alignments_file = "source_data/cabili_full_alignments.fasta"
# bed_file = "source_data/cabili_clean_transcripts.bed"
#
# output_exon_file = "clean_run/genome_sequences/lincrna/cabili/exon_alignments.fasta"
# output_intron_file = "clean_run/genome_sequences/lincrna/cabili/intron_alignments.fasta"



bed_entries = gen.read_many_fields(bed_file, "\t")
entries = collections.defaultdict(lambda: collections.defaultdict())
for entry in bed_entries:
    entries[int(entry[1])][int(entry[2])] = entry


with open(output_exon_file, "w") as exon_outfile:
    with open(output_intron_file, "w") as intron_outfile:
        with open(alignments_file, "r") as input_file:
            lines = input_file.readlines()
            lines = [i.strip("\n") for i in lines]
            for entry_id in range(0, len(lines), 5):
                if entry_id:
                    human_name = lines[entry_id]
                    human_seq = lines[entry_id+1]
                    mac_name = lines[entry_id+2]
                    mac_seq = lines[entry_id+3]

                    # ensure that there isnt too many indels
                    if np.divide(len([i for i in human_seq if i == "-"]), len(human_seq)) < 0.15 and np.divide(len([i for i in mac_seq if i == "-"]), len(mac_seq)) < 0.15:
                        coordinates = [int(i) for i in human_name.split(":")[-1].split("-")]
                        bed_entry = entries[coordinates[0]][coordinates[1]]
                        id = bed_entry[3]
                        exon_sizes = [int(i) for i in bed_entry[-2].split(",") if len(i)]
                        exon_starts = [int(i) for i in bed_entry[-1].split(",") if len(i)]
                        exon_coordinates = [[start, start + exon_sizes[i]] for i, start in enumerate(exon_starts)]
                        intron_coordinates = [[start + exon_sizes[i], exon_starts[i+1]] for i, start in enumerate(exon_starts) if i+1 < len(exon_starts)]

                        human_exon_sequences = [human_seq[i[0]:i[1]].upper() for i in exon_coordinates]
                        mac_exon_sequences = [mac_seq[i[0]:i[1]].upper() for i in exon_coordinates]
                        [exon_outfile.write(">{0}.{1}\n{2},{3}\n".format(id, i+1, human_exon_sequences[i], mac_exon_sequences[i])) for i in range(len(human_exon_sequences))]

                        human_intron_sequences = [human_seq[i[0]:i[1]].upper() for i in intron_coordinates]
                        mac_intron_sequences = [mac_seq[i[0]:i[1]].upper() for i in intron_coordinates]
                        [intron_outfile.write(">{0}.{1}\n{2},{3}\n".format(id, i+1, human_intron_sequences[i], mac_intron_sequences[i])) for i in range(len(human_intron_sequences))]





# alignment_entries = gen.read_many_fields(alignments_file, "\t")
# for i in range(0, len(alignment_entries), 5):
#     if i < 5:
#         print(alignment_entries[i])
