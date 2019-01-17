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


def locate_random_hexamers(iterations, seq, motifs, output_directory):
    choices = list(range(len(seq)))
    output_files = []
    for i, iteration in enumerate(iterations):
        np.random.seed()
        gen.print_parallel_status(i, iterations)

        test_set = []
        generated = False
        while not generated:
            required = len(motifs) - len(test_set)
            chosen = np.random.choice(choices, len(motifs) - len(test_set))
            test_set.extend(list(set([seq[i:i+6] for i in chosen if "X" not in seq[i:i+6] and seq[i:i+6] not in test_set])))
            if len(test_set) == len(motifs):
                generated = True

        output_file = "{0}/{1}.txt".format(output_directory, random.random())
        with open(output_file, "w") as outfile:
            [outfile.write("{0}\n".format(i)) for i in test_set]

    return output_files




def run_hexamer_test(simulatioms, motif_file, sequences_file, controls_dir, output_file, families_file = None):

    gen.create_output_directories(controls_dir)

    motifs = sequo.read_motifs(motif_file)

    if len(os.listdir(controls_dir)) < simulations:
        gen.create_output_directories(controls_dir)
        required = list(range(simulations - len(os.listdir(controls_dir))))

        names, seqs = gen.read_fasta(sequences_file)

        seqs_list = collections.defaultdict(lambda: [])
        for i, name in enumerate(names):
            seqs_list[name.split(".")[0]].append(seqs[i])

        if families_file:
            seqs_list = sequo.pick_random_family_member(families_file, seqs_list)

        all_seqs = []
        [all_seqs.extend(seqs_list[i]) for i in seqs_list]
        full_seq = "X".join(all_seqs)
        soc.run_simulation_function(required, [full_seq, motifs, controls_dir], locate_random_hexamers, sim_run = False)


    real_purine_content = sequo.calc_purine_content(motifs)
    real_nt_content = sequo.calc_nucleotide_content(motifs)

    test_purine_content = []
    test_nt_content = []
    for file in os.listdir(controls_dir):
        filepath = "{0}/{1}".format(controls_dir, file)
        test_motifs = sequo.read_motifs(filepath)
        test_purine_content.append(sequo.calc_purine_content(test_motifs))
        test_nt_content.append(sequo.calc_nucleotide_content(test_motifs))

    with open(output_file, "w") as outfile:
        outfile.write("id,purine_content,a_content,c_content,g_content,t_content\n")
        outfile.write("real,{0},{1}\n".format(real_purine_content, ",".join(gen.stringify([real_nt_content[i] for i in sorted(real_nt_content)]))))
        for i in range(len(test_purine_content)):
            outfile.write("{0},{1},{2}\n".format(i+1, test_purine_content[i], ",".join(gen.stringify([test_nt_content[i][j] for j in sorted(test_nt_content[i])]))))

simulations = 1000

# intron_hexamer_dir = "clean_run/intron_hexamers"
# intron_output_dir = "clean_run/tests/introns"
# gen.create_output_directories(intron_output_dir)
# intron_output_file = "{0}/intron_hexamers.csv".format(intron_output_dir)
# run_hexamer_text(ese_file, introns_file, families_file = families_file)

exon_hexamer_dir = "clean_run/exon_hexamers"
exon_output_dir = "clean_run/tests/exons"
gen.create_output_directories(exon_output_dir)
exon_output_file = "{0}/exon_hexamers.csv".format(exon_output_dir)
run_hexamer_test(simulations, ese_file, exons_file, exon_hexamer_dir, exon_output_file, families_file = families_file)



# run_hexamer_text(simulations, ese_file, exons_file, intron_hexamer_dir, output_file, families_file = families_file)
