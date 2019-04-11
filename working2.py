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
import scipy.stats
# import matplotlib.pyplot as plt
import conservation
import os
import zipfile
from useful_motif_sets import stops
from regex_patterns import codon_pattern
import containers as cont
import itertools as it
import random
import copy



# file = "clean_run/genome_sequences/human/human.cds.clean_coding_exons.fasta"
# cds_file = "clean_run/genome_sequences/human/human.cds.multi_exons.fasta"
# families_file = "clean_run/genome_sequences/human/human.cds.families.bed"
file = "clean_run/genome_sequences/lincrna/cabili/multi_exons.fasta"
families_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt"


motif_set = "int3"
motif_file = "source_data/motif_sets/{0}.txt".format(motif_set)
motifs = sequo.read_motifs(motif_file)


sim_path = "clean_run/dinucleotide_controls/{0}_dinucleotide_controls_matched_stops".format(motif_set)
sims = gen.get_filepaths(sim_path)[:100]


names, seqs = gen.read_fasta(file)
families = gen.read_many_fields(families_file, "\t")
seqs = {name.split(".")[0]: seqs[i] for i, name in enumerate(names)}

seqs = sequo.group_family_results(seqs, families)
seqs = [seqs[i][0] for i in seqs]
seqs = seqs


def calc_overlaps(file):
    motifs = sequo.read_motifs(file)
    new_seqs = []
    for seq in seqs:
        overlaps = sequo.sequence_overlap_indicies(seq, motifs)
        new_seq = "".join([nt if i not in overlaps else "C" for i, nt in enumerate(list(seq))])
        new_seqs.append(new_seq)
    longest_orfs = seqo.get_longest_orfs(new_seqs)
    return longest_orfs

real_longest_orfs = seqo.get_longest_orfs(seqs)

# sim_orfs = []
# for file in sims:
#     sim_orfs.append(calc_overlaps(file))
#
# zs = []
# for i, length in enumerate(real_longest_orfs):
#     sims = [sim[i] for sim in sim_orfs]
#     z = np.divide(length - np.mean(sims), np.std(sims))
#     zs.append(z)

print(real_longest_orfs)

# print(len([z for z in zs if z < 0]), len(zs))
# print(len([z for z in zs if z < 0 and z < -1.96]), len(zs))
