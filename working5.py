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


ese_set = "combined_eses"
motif_file = "source_data/motif_sets/{0}.txt".format(ese_set)
controls_dir = "clean_run/dinucleotide_controls/{0}_dinucleotide_controls_matched_stops".format(ese_set)
seqs_file = "clean_run/genome_sequences/lincrna/cabili/multi_exons.fasta"
families_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt"

names, seqs = gen.read_fasta(seqs_file)
exons = collections.defaultdict(lambda: [])
[exons[name.split(".")[0]].append(seqs[i]) for i, name in enumerate(names)]

exons = sequo.pick_random_family_member(families_file, exons)

exon_list = []
[[exon_list.append(i) for i in exons[id] if len(i) >= 211] for id in exons]

cores = []
flanks = []

for exon in exon_list:
    flanks.append(exon[2:69])
    flanks.append(exon[-69:-2])

    midpoint = int(len(exon) / 2)
    core = exon[midpoint-33:midpoint+34]
    cores.append(core)

core_count = seqo.calc_motif_counts(cores, stops)
core_nts = len("".join(cores))
print(core_count)
print(core_nts)
core_rate = np.divide(core_count*3, core_nts)
print(core_rate)

flank_count = seqo.calc_motif_counts(flanks, stops)
flank_nts = len("".join(flanks))
expected_flanks = core_rate*flank_nts
print(flank_count, expected_flanks)
print(np.divide(flank_count - expected_flanks, expected_flanks))
