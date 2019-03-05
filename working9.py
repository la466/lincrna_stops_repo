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
# seqs_file = "clean_run/lincrna/lincRNA.multi_exon.exons.fasta"
seqs_file = "clean_run/genome_sequences/human/human.cds.clean_coding_exons.fasta"

motifs = sequo.read_motifs(motifs)

names, seqs = gen.read_fasta(seqs_file)
seq_list = {name: seqs[i] for i, name in enumerate(names)}


nts = 0
all_overlaps1 = []
all_overlaps2 = []

for id in seq_list:
    sequence = seq_list[id]
    three_prime_flank = sequence[:200]
    nts += len(three_prime_flank)
    overlaps = sequo.sequence_overlap_indicies(three_prime_flank, ["CCC", "TTT", "CGT"])
    # overlaps = [len(three_prime_flank) - i for i in overlaps]
    all_overlaps2.extend(overlaps)

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

densities = []
ilist = []
for i in range(200):
    ilist.append(i)
    density = np.divide(all_overlaps2.count(i), nts)
    densities.append(density)

fig, ax = plt.subplots()
ax.plot(ilist, densities)
ax.set(xlabel='Position', ylabel='Density')
ax.grid()
plt.show()
