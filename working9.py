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
lseqs_file = "clean_run/genome_sequences/lincrna/cabili/multi_exons.fasta"
lseqs_families_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt"

names, seqs = gen.read_fasta(lseqs_file)
lincrna_seq_list = collections.defaultdict(lambda: [])
[lincrna_seq_list[name.split(".")[0]].append(seqs[i]) for i, name in enumerate(names) if len(seqs[i])]

lincrna_seq_list = sequo.pick_random_family_member(lseqs_families_file, lincrna_seq_list)

seqs_file = "clean_run/genome_sequences/human/human.cds.clean_coding_exons.fasta"
families_file = "clean_run/genome_sequences/human/human.cds.families.bed"

motifs = sequo.read_motifs(motifs)

names, seqs = gen.read_fasta(seqs_file)
seq_list = collections.defaultdict(lambda: [])
[seq_list[name.split(".")[0]].append(seqs[i]) for i, name in enumerate(names) if len(seqs[i])]

seq_list = sequo.pick_random_family_member(families_file, seq_list)


stop_motifs = [i for i in motifs if len(re.findall("(?=(TAA|TAG|TGA))", i))]
non_stop_motifs = [i for i in motifs if i not in stop_motifs]

stop_densities = []
non_stop_densities = []
relative_densities = []
for id in seq_list:
    stop_density = seqo.calc_motif_density(seq_list[id], stop_motifs)
    non_stop_density = seqo.calc_motif_density(seq_list[id], non_stop_motifs)
    stop_densities.append(stop_density)
    non_stop_densities.append(non_stop_density)
    relative_densities.append(np.divide(stop_density, non_stop_density))

lstop_densities = []
lnon_stop_densities = []
lrelative_densities = []
for id in lincrna_seq_list:
    stop_density = seqo.calc_motif_density(lincrna_seq_list[id], stop_motifs)
    non_stop_density = seqo.calc_motif_density(lincrna_seq_list[id], non_stop_motifs)
    lstop_densities.append(stop_density)
    lnon_stop_densities.append(non_stop_density)
    lrelative_densities.append(np.divide(stop_density, non_stop_density))

# print(scipy.stats.wilcoxon(stop_densities, non_stop_densities))
print(relative_densities)
print(lrelative_densities)
print(scipy.stats.ranksums(relative_densities, lrelative_densities))
print(np.nanmean(relative_densities))
print(np.nanmean(lrelative_densities))
