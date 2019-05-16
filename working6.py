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

### masking eses for orf lengths



# file = "clean_run/genome_sequences/human/human.cds.clean_coding_exons.fasta"
# cds_file = "clean_run/genome_sequences/human/human.cds.multi_exons.fasta"
# families_file = "clean_run/genome_sequences/human/human.cds.families.bed"
file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_transcript_sequences.fasta"
families_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt"

# test_file = "clean_run/tests/lincrna/orf_length_sim/cabili_sim_orf_lengths_zs_grouped.csv"

# entries = gen.read_many_fields(test_file, ",")[1:]
# entries = {i[0]: float(i[7]) for i in entries if float(i[7]) > 0}



# motif_set = "int3"
# motif_file = "source_data/motif_sets/{0}.txt".format(motif_set)
# motifs = sequo.read_motifs(motif_file)
#

# sim_path = "clean_run/dinucleotide_controls/{0}_dinucleotide_controls_matched_stops".format(motif_set)
# sims = gen.get_filepaths(sim_path)[:100]


names, seqs = gen.read_fasta(file)
seqs = {name.split(".")[0]: seqs[i] for i, name in enumerate(names)}
# families = gen.read_many_fields(families_file, "\t")
# seqs = {name.split(".")[0]: seqs[i] for i, name in enumerate(names)}

# lengths = {i: seqo.get_longest_orfs([seqs[i]])[0] for i in seqs}
# lengths = sequo.group_family_results(lengths, families)
# lengths = {i: np.nanmedian(lengths[i]) for i in lengths}
# print(seqs)


# seqs = {i: seqs[i] for i in seqs if i in entries}
# seqs = [seqs[i][0] for i in seqs]
# print(len(seqs))
# seqs = seqs[:500]
seqs_string = "".join([seqs[i] for i in seqs])
seqs = sequo.pick_random_family_member(families_file, seqs)
seqs = [seqs[i] for i in seqs]
seqs_lengths = [len(i) for i in seqs]

longest_orfs = seqo.get_longest_orfs(seqs)


def calc_orfs(sims, seqs_string, seqs_lengths):
    outputs = {}
    if len(sims):
        np.random.seed()
        for i, sim in enumerate(sims):
            gen.print_parallel_status(i, sims)
            sim_nts = list(seqs_string)
            np.random.shuffle(sim_nts)
            sim_string = "".join(sim_nts)
            new_seqs = []
            index = 0
            for length in seqs_lengths:
                new_seqs.append(sim_string[index:index+length])
            sim_orfs = seqo.get_longest_orfs(new_seqs)
            outputs[sim] = sim_orfs

    return outputs

outputs = soc.run_simulation_function(list(range(1000)), [seqs_string, seqs_lengths], calc_orfs, sim_run = False)
#
output_dir = "clean_run/tests/lincrna/total_lengths1"
gen.create_output_directories(output_dir)

for threshold in list(range(200,610,10)):

    real_greater = len([i for i in longest_orfs if i >= threshold])
    output_file = "{0}/{1}.csv".format(output_dir, threshold)
    with open(output_file, "w") as outfile:
        outfile.write("id,greater,max_orf,total\n")
        outfile.write("real,{0},{1},{2}\n".format(real_greater, max(longest_orfs), len(longest_orfs)))
        for i in outputs:
            sim_output = outputs[i]
            outfile.write("sim_{0},{1},{2},{3}\n".format(i, len([j for j in sim_output if j >= threshold]), max(sim_output), len(sim_output)))
