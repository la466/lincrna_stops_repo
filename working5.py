import sim_ops as simo
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


ese_set = "int3"
motif_file = "source_data/motif_sets/{0}.txt".format(ese_set)
controls_dir = "clean_run/dinucleotide_controls/{0}_dinucleotide_controls_matched_stops".format(ese_set)
seqs_file = "clean_run/genome_sequences/human/human.cds.clean_non_coding_exons.fasta"
families_file = "clean_run/genome_sequences/human/human.cds.families.bed"
# seqs_file = "clean_run/genome_sequences/lincrna/cabili/single_exons.fasta"
# families_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt"

motifs = sequo.read_motifs(motif_file)

names, seqs = gen.read_fasta(seqs_file)
exons = collections.defaultdict(lambda: [])
[exons[name.split(".")[0]].append(seqs[i]) for i, name in enumerate(names)]

# exons = sequo.pick_random_family_member(families_file, exons)
exons = {i: exons[i][0] for i in exons}

# print()
#
# exon_list = []
# [[exon_list.append(i) for i in exons[id] if len(i)] for id in exons]
#
# density = seqo.calc_motif_density(exon_list, stops)


simulations = ["real"] + list(range(10))
output_directory = "temp_files_working_single"
gen.remove_directory(output_directory)
gen.create_output_directories(output_directory)
output_filelist = []

outputs = simo.simulate_lincrna_stop_codon_density(simulations, exons, output_directory, output_filelist, inverse = None)

sims = []
for id in outputs:
    # print(id)
    density = float(gen.read_many_fields(outputs[id], ",")[0][-1])
    # print(density)
    if id == "real":
        # print("here")
        # print(density)
        real = density
        # print(density)
    else:
        sims.append(density)

print(real)
print(sims)
print(np.divide(real - np.mean(sims), np.mean(sims)))
gen.remove_directory(output_directory)

#
# cores = []
# flanks = []
#
# for exon in exon_list:
#     # flanks.append(exon[2:69])
#     flanks.append(exon[-69:-2])
#
#     midpoint = int(len(exon) / 2)
#     core = exon[midpoint-33:midpoint+34]
#     cores.append(core)
#
# core_density = seqo.calc_motif_density(cores, motifs)
# flank_density = seqo.calc_motif_density(flanks, motifs)
#
# core_stop_density = seqo.calc_motif_density(cores, stops)
# flank_stop_density = seqo.calc_motif_density(flanks, stops)
#
# print(core_density, flank_density)
# print(core_stop_density, flank_stop_density)


# core_count = seqo.calc_motif_counts(cores, stops)
# core_nts = len("".join(cores))
# print(core_count)
# print(core_nts)
# core_rate = np.divide(core_count*3, core_nts)
# print(core_rate)
#
# flank_count = seqo.calc_motif_counts(flanks, stops)
# flank_nts = len("".join(flanks))
# expected_flanks = core_rate*flank_nts
# print(flank_count, expected_flanks)
# print(np.divide(flank_count - expected_flanks, expected_flanks))
