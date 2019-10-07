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
from useful_motif_sets import stops, codon_map
from regex_patterns import codon_pattern
import containers as cont
import itertools as it
import random
import copy
from scipy.stats import chisquare
import collections


seqs_file = "clean_run/genome_sequences/lincrna/cabili/transcript_sequences.fasta"
families_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt"
gc_sets = "clean_run/gc_matched_sets.bed"

codon_sets = gen.read_many_fields(gc_sets, "\t")
# codon_sets = codon_sets[:50]

names, seqs = gen.read_fasta(seqs_file)
seq_list = {name: seqs[i] for i, name in enumerate(names) if "N" not in seqs[i]}
seq_list = sequo.pick_random_family_member(families_file, seq_list)


def randomise_seq(seq):
    nts = list(seq)
    np.random.shuffle(nts)
    return "".join(nts)


output_directory = "temp_shuffle_linc"
gen.create_output_directories(output_directory)

def run_simulations(iterations, seq_list, codon_sets, output_directory):
    outputs = []
    if len(iterations) > 0:
        np.random.seed()

        for i, iteration in enumerate(iterations):
            print("{0}/{1}".format(i + 1, len(iterations)))

            if iteration != "real":
                new_seqs = []
                for id in seq_list:
                    new_seqs.append(randomise_seq(seq_list[id]))
            else:
                new_seqs = [seq_list[i] for i in seq_list]

            iteration_densities = {}
            for codon_set in codon_sets:
                iteration_densities["_".join(codon_set)] = seqo.calc_motif_density(new_seqs, codon_set)

            output_file = "{0}/{1}.csv".format(output_directory, iteration)
            with open(output_file, "w") as outfile:
                for codon_set in codon_sets:
                    outfile.write("{0},{1}\n".format("_".join(codon_set), iteration_densities["_".join(codon_set)]))
    return outputs

iterations = ["real"] + list(range(1000))
soc.run_simulation_function(iterations, [seq_list, codon_sets, output_directory], run_simulations, sim_run = False)
