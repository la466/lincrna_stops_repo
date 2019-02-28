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

eses = sequo.read_motifs("source_data/motif_sets/int3.txt")
# eses = sequo.read_motifs("source_data/motif_sets/RESCUE.txt")
# eses = sequo.read_motifs("source_data/motif_sets/ke400.txt")

fileset = {
    "fileset_prefix": "pc",
    "exons_file": "clean_run/genome_sequences/human/human.cds.clean_coding_exons.fasta",
    "introns_file": "clean_run/genome_sequences/human/human.clean_introns.fasta",
    "families_file": "clean_run/genome_sequences/human/human.cds.families.bed" ,
}

lincrna_files = {
    "fileset_prefix": "lincrna",
    "exons_file": "clean_run/genome_sequences/lincrna/cabili/exons.fasta",
    "introns_file": "clean_run/genome_sequences/lincrna/cabili/introns.fasta",
    "families_file": "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt" ,
}



exon_names, exon_seqs = gen.read_fasta(lincrna_files["exons_file"])
exon_list = collections.defaultdict(lambda: [])
for i, name in enumerate(exon_names):
    exon_list[name.split(".")[0]].append(exon_seqs[i])

exon_list = sequo.pick_random_family_member(lincrna_files["families_file"], exon_list)
temp_exon_list = []
for id in exon_list:
    for seq in exon_list[id]:
        # temp_exon_list.append(seq[2:69])
        # temp_exon_list.append(seq[-69:-2])
        temp_exon_list.append(seq)
# [temp_exon_list.extend(exon_list[i]) for i in exon_list]
exon_list = temp_exon_list


stop_eses = [i for i in eses if len(re.findall("(?=(TAA|TAG|TGA))", i))]
non_stop_eses = [i for i in eses if i not in stop_eses]

stop_motif_hits = gen.flatten([sequo.return_overlap_motifs(seq, stop_eses) for seq in exon_list])
non_stop_motif_hits = gen.flatten([sequo.return_overlap_motifs(seq, non_stop_eses) for seq in exon_list])

stop_overlap_hits = [i for i in stop_motif_hits if len(i) > 6]
stop_overlap_hits_length = len([i for i in stop_motif_hits if len(i) > 6])
non_stop_overlap_hits = [i for i in non_stop_motif_hits if len(i) > 6]
non_stop_overlap_hits_length = len([i for i in non_stop_motif_hits if len(i) > 6])
total_overlap_hits = stop_overlap_hits_length + non_stop_overlap_hits_length
total_hits = len(stop_motif_hits) + len(non_stop_motif_hits)

prop_overlaps = np.divide(total_overlap_hits, total_hits)

expected_stops = prop_overlaps*len(stop_motif_hits)
expected_non_stops = prop_overlaps*len(non_stop_motif_hits)


p = []
groups = [[expected_stops, stop_overlap_hits_length],[expected_non_stops, non_stop_overlap_hits_length]]
for group in groups:
    oe = np.divide((group[1] - group[0])**2, group[0])
    p.append(oe)

print(sum(p))
print(stop_overlap_hits_length, expected_stops, np.divide(stop_overlap_hits_length, expected_stops))
print(non_stop_overlap_hits_length, expected_non_stops, np.divide(non_stop_overlap_hits_length, expected_non_stops))

print(np.mean([len(i) for i in stop_overlap_hits]))
print(np.mean([len(i) for i in non_stop_overlap_hits]))
# print(np.divide(len(stop_overlap_hits), len(stop_motif_hits)))
# print(np.divide(len(non_stop_overlap_hits), len(non_stop_motif_hits)))
