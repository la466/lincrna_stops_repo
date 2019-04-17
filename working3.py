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
t_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_transcript_sequences1.fasta"
families_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt"

all = collections.defaultdict(lambda: [])
names, seqs = gen.read_fasta(seq_file)
[all[name.split(".")[0]].append(seqs[i]) for i, name in enumerate(names)]

exons = collections.defaultdict(lambda: collections.defaultdict())
for i, name in enumerate(names):
    exons[name.split(".")[0]][int(name.split(".")[1].split("(")[0])] = seqs[i]

tall = collections.defaultdict(lambda: [])
tnames, tseqs = gen.read_fasta(t_file)
[tall[name.split(".")[0]].append(tseqs[i]) for i, name in enumerate(tnames)]


tall = sequo.pick_random_family_member(families_file, tall)
print(len(tall))

# retained = []
# for name in all:
#     test = [i for i in all[name] if "N" in i]
#     if len(test) == 0:
#         retained.append(name)

# output_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_transcript_sequences1.fasta"
# with open(output_file, "w") as outfile:
#     for name in tall:
#         if name in retained:
#             # for exon in sorted(tall[name]):
#             outfile.write(">{0}\n{1}\n".format(name, tall[name]))
