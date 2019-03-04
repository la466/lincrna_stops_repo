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
from scipy.stats import chisquare


motifs = "source_data/motif_sets/int3.txt"
seqs_file = "clean_run/lincrna/lincRNA.multi_exon.fasta"

motifs = sequo.read_motifs(motifs)

names, seqs = gen.read_fasta(seqs_file)
seq_list = {name: seqs[i] for i, name in enumerate(names[:50])}

def calc_non_hit_densities(sequences):

    non_hit_motifs = []
    for sequence in sequences:
        non_hits = sequo.return_overlap_motifs(sequence, motifs, inverse = True)
        non_hit_motifs.extend(non_hits)
    density = seqo.calc_motif_density(non_hit_motifs, stops)
    print(density)


sim_seqs = {}

for sim in list(range(3)):
    simulated_seqs = []
    for id in seq_list:
        seq = seq_list[id]
        nts = list(seq)
        np.random.shuffle(nts)
        simulated_seqs.append("".join(nts))
    sim_seqs[sim] = simulated_seqs

calc_non_hit_densities([seq_list[i] for i in seq_list])
[calc_non_hit_densities(sim_seqs[sim]) for sim in sim_seqs]
