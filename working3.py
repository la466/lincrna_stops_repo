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


seq_file = "clean_run/genome_sequences/lincrna/cabili/multi_exons.fasta"
seqs_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_transcript_sequences.fasta"
families_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt"

def calc_sim_densities(sims, sim_seqs):
    outputs = []
    if len(sims):
        for sim in sims:
            sim_seqs = randomise_seqs(kept)
            outputs.append(seqo.calc_motif_density(sim_seqs, stops))
    return outputs

def randomise_seqs(seqs):
    new_seqs = []
    for seq in seqs:
        nts = list(seq)
        np.random.shuffle(nts)
        new_seqs.append("".join(nts))

    return new_seqs

names, seqs = gen.read_fasta(seqs_file)
# all = {name.split(".")[0]: seqs[i] for i, name in enumerate(names)}

all = dict(zip([name.split(".")[0] for name in names], seqs))

all = sequo.pick_random_family_member(families_file, all)
# print(len(tnames))

with open("clean_run/tests/lincrna/stop_density/upstream_atg_stop_density.csv", "w") as outfile:
    outfile.write("threshold,real_density,median_sim_density,p_value\n")
    for length in list(range(10, 110, 10)):
        print(length)
        kept = []
        for id in all:
            seq = all[id]
            try:
                atg = seq.index("ATG")
                if atg > length:
                    kept.append(seq[:atg])
            except:
                pass



        real_density = seqo.calc_motif_density(kept, stops)

        sim_densities  = []
        sims = list(range(1000))

        sim_densities = soc.run_simulation_function(sims, [all], calc_sim_densities, sim_run = False)
        less_than = [i for i in sim_densities if i <= real_density]

        output = [length, real_density, np.median(sim_densities), np.divide(len(less_than) + 1, len(sims) + 1)]
        outfile.write("{0}\n".format(",".join(gen.stringify(output))))
