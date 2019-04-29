import generic as gen
import seq_ops as seqo
import sequence_ops as sequo
import sim_ops_containers as soc
import sim_ops as simo
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
# seqs_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_transcript_sequences.fasta"
seqs_file = "clean_run/genome_sequences/lincrna/cabili/multi_exons.fasta"
introns_file = "clean_run/genome_sequences/lincrna/cabili/introns.fasta"
families_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt"

names, seqs = gen.read_fasta(seqs_file)
exon_list = collections.defaultdict(lambda: [])
[exon_list[name.split(".")[0]].append(seqs[i]) for i, name in enumerate(names)]
# exon_list = {name.split(".")[0]: seqs[i] for i, name in enumerate(names)}
exon_list = sequo.pick_random_family_member(families_file, exon_list)

# seq_list = []
# [seq_list.extend(exon_list[i]) for i in exon_list]
# seq_list = seq_list[]
# print(seq_list)

intron_names, intron_seqs = gen.read_fasta(introns_file)
intron_list = collections.defaultdict(lambda: [])
[intron_list[name.split(".")[0]].append(intron_seqs[i]) for i, name in enumerate(intron_names) if name.split(".")[0] in exon_list]


intron_list = {i: "".join(intron_list[i]) for i in intron_list}


# real_density = seqo.calc_motif_density(introns, stops)
#
# greater_than = []
#
# for id in exon_list:
#     # exons = "".join(exon_list[id])
#     # introns = "".join(intron_list[id])]
#     exon_density = seqo.calc_motif_density(exon_list[id], stops)
#     intron_density = seqo.calc_motif_density(intron_list[id], stops)
#     # print(exon_density, intron_density)
#     if intron_density > real_density:
#         greater_than.append(id)
#
# print(len(greater_than))

simulations = ["real", 1, 2, 3, 4, 5]
output_directory = "temp_files_working"
gen.remove_directory(output_directory)
gen.create_output_directories(output_directory)
output_filelist = []

outputs = soc.run_simulation_function(simulations, [intron_list, output_directory, output_filelist], simo.simulate_lincrna_stop_codon_density, sim_run = False)
print(outputs)
