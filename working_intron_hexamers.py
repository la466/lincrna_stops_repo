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
from useful_motif_sets import stops
from regex_patterns import codon_pattern
import containers as cont
import itertools as it
import random

purine = ["A", "G"]
pyrimidine = ["C", "T"]

exons_file = "clean_run/genome_sequences/human/human.cds.clean_coding_exons.fasta"
introns_file = "clean_run/genome_sequences/human/human.clean_introns.fasta"
families_file = "clean_run/genome_sequences/human/human.cds.families.bed"
ese_file = "source_data/motif_sets/int3.txt"

final_output_dir = "clean_run/tests/introns"
gen.create_output_directories(final_output_dir)
output_file = "{0}/intron_hexamers.csv".format(final_output_dir)

eses = sequo.read_motifs(ese_file)

exon_names, exon_seqs = gen.read_fasta(exons_file)
intron_names, intron_seqs = gen.read_fasta(introns_file)

introns = collections.defaultdict(lambda: [])
for i, name in enumerate(intron_names):
    introns[name.split(".")[0]].append(intron_seqs[i])


introns = sequo.pick_random_family_member(families_file, introns)

intron_seqs = []
[intron_seqs.extend(introns[i]) for i in introns]

intron_seq = "X".join(intron_seqs)
choices = list(range(len(intron_seq)))

def locate_random_hexamers(iterations, seq, output_directory):
    output_files = []
    for i, iteration in enumerate(iterations):
        np.random.seed()
        gen.print_parallel_status(i, iterations)

        test_set = []
        generated = False
        while not generated:
            required = len(eses) - len(test_set)
            chosen = np.random.choice(choices, len(eses) - len(test_set))
            test_set.extend(list(set([intron_seq[i:i+6] for i in chosen if "X" not in intron_seq[i:i+6] and intron_seq[i:i+6] not in test_set])))
            if len(test_set) == len(eses):
                generated = True

        output_file = "{0}/{1}.txt".format(output_directory, random.random())
        with open(output_file, "w") as outfile:
            [outfile.write("{0}\n".format(i)) for i in test_set]

    return output_files


output_dir = "clean_run/intron_hexamers"
gen.create_output_directories(output_dir)

simulations = 1000
if len(os.listdir(output_dir)) < simulations:
    required = list(range(simulations - len(os.listdir(output_dir))))
    soc.run_simulation_function(required, [intron_seq, output_dir], locate_random_hexamers, sim_run = False)

real_purine_content = sequo.calc_purine_content(eses)
real_nt_content = sequo.calc_nucleotide_content(eses)

test_purine_content = []
test_nt_content = []
for file in os.listdir(output_dir):
    filepath = "{0}/{1}".format(output_dir, file)
    motifs = sequo.read_motifs(filepath)
    test_purine_content.append(sequo.calc_purine_content(motifs))
    test_nt_content.append(sequo.calc_nucleotide_content(motifs))

with open(output_file, "w") as outfile:
    outfile.write("id,purine_content,a_content,c_content,g_content,t_content\n")
    outfile.write("real,{0},{1}\n".format(real_purine_content, ",".join(gen.stringify([real_nt_content[i] for i in sorted(real_nt_content)]))))
    for i in range(len(test_purine_content)):
        outfile.write("{0},{1},{2}\n".format(i+1, test_purine_content[i], ",".join(gen.stringify([test_nt_content[i][j] for j in sorted(test_nt_content[i])]))))
