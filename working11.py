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




motifs = "source_data/motif_sets/int3.txt"
# motifs = "source_data/motif_sets/RESCUE.txt"
# seqs_file = "clean_run/lincrna/lincRNA.multi_exon.exons.fasta"
# seqs_file = "clean_run/genome_sequences/lincrna/cabili/multi_exons.fasta"

families_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt"
alignments_file = "clean_run/genome_sequences/lincrna/cabili/exon_alignments.fasta"
# intron_file = "clean_run/genome_sequences/lincrna/cabili/intron_alignments.fasta"

families_file = "clean_run/genome_sequences/human/human.cds.families.bed"
alignments_file = "clean_run/genome_sequences/human.macaque.alignments.fasta"

exon_list = collections.defaultdict(lambda: collections.defaultdict())
exon_names, exon_seqs = gen.read_fasta(alignments_file)
for i, name in enumerate(exon_names):
    alignments = exon_seqs[i].split(",")
    if len(alignments[0]) > 207:
        exon_list[name.split(".")[0]] = alignments

# intron_list = collections.defaultdict(lambda: collections.defaultdict())
# intron_names, intron_seqs = gen.read_fasta(intron_file)
# for i, name in enumerate(intron_names):
#     alignments = intron_seqs[i].split(",")
#     id = name.split(".")[0]
#     if len(alignments[0]) and id in exon_list:
#         intron_list[id][int(name.split(".")[-1])] = alignments

# now retain only the exons with retained introns
# exon_list = {id: exon_list[id] for id in exon_list if id in intron_list}
exon_list = sequo.pick_random_family_member(families_file, exon_list)
print(len(exon_list))

with open("temp_files/core_rates.csv", "w") as outfile:
    outfile.write("id,exon_flank,exon_core,stop_f_rate,stop_c_rate\n")

    human_flanks = []
    mac_flanks = []

    human_cores = []
    mac_cores = []

    for i, id in enumerate(exon_list):
        human_exon_flanks = []
        human_exon_cores = []
        mac_exon_flanks = []
        mac_exon_cores = []

        human_exon = exon_list[id][0]
        mac_exon = exon_list[id][1]

        human_flanks.extend([human_exon[2:69], human_exon[-69:-2]])
        mac_flanks.extend([mac_exon[2:69], mac_exon[-69:-2]])

        midpoint = int(len(human_exon) / 2)
        human_cores.append(human_exon[midpoint-33:midpoint+34])
        mac_cores.append(mac_exon[midpoint-33:midpoint+34])

    joined_human_exon_flanks = "".join(human_flanks)
    joined_human_exon_cores = "".join(human_cores)
    joined_mac_exon_flanks = "".join(mac_flanks)
    joined_mac_exon_cores = "".join(mac_cores)

    # print(joined_human_exon_flanks)

    # if len(joined_human_exon_flanks) > 100 and len(joined_human_exon_cores) > 100:

    exon_flank_rate = sequo.get_sub_rate(joined_human_exon_flanks, joined_mac_exon_flanks)
    exon_core_rate = sequo.get_sub_rate(joined_human_exon_cores, joined_mac_exon_cores)

    human_flank_stop_overlaps = sequo.sequence_overlap_indicies(joined_human_exon_flanks, stops)
    human_flank_non_stop_overlaps = [i for i in range(len(joined_human_exon_flanks)) if i not in human_flank_stop_overlaps]
    human_flank_stops = "".join([joined_human_exon_flanks[i] for i in human_flank_stop_overlaps])
    human_flank_non_stops = "".join([joined_human_exon_flanks[i] for i in human_flank_non_stop_overlaps])
    mac_flank_stops = "".join([joined_mac_exon_flanks[i] for i in human_flank_stop_overlaps])
    mac_flank_non_stops = "".join([joined_mac_exon_flanks[i] for i in human_flank_non_stop_overlaps])

    human_core_stop_overlaps = sequo.sequence_overlap_indicies(joined_human_exon_cores, stops)
    human_core_non_stop_overlaps = [i for i in range(len(joined_human_exon_cores)) if i not in human_core_stop_overlaps]
    human_core_stops = "".join([joined_human_exon_cores[i] for i in human_core_stop_overlaps])
    human_core_non_stops = "".join([joined_human_exon_cores[i] for i in human_core_non_stop_overlaps])
    mac_core_stops = "".join([joined_mac_exon_cores[i] for i in human_core_stop_overlaps])
    mac_core_non_stops = "".join([joined_mac_exon_cores[i] for i in human_core_non_stop_overlaps])

    print(exon_flank_rate, exon_core_rate)

    flank_stop_rate = sequo.get_sub_rate(human_flank_stops, mac_flank_stops)
    flank_non_stop_rate = sequo.get_sub_rate(human_flank_non_stops, mac_flank_non_stops)
    core_stop_rate = sequo.get_sub_rate(human_core_stops, mac_core_stops)
    core_non_stop_rate = sequo.get_sub_rate(human_core_non_stops, mac_core_non_stops)

    print(flank_stop_rate, flank_non_stop_rate, core_stop_rate, core_non_stop_rate, np.divide(flank_stop_rate,flank_non_stop_rate), np.divide(core_stop_rate, core_non_stop_rate))
