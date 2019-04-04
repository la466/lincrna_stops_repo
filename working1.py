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



file = "clean_run/genome_sequences/human/human.cds.clean_coding_exons.fasta"
cds_file = "clean_run/genome_sequences/human/human.cds.multi_exons.fasta"
families_file = "clean_run/genome_sequences/human/human.cds.families.bed"
# file = "clean_run/genome_sequences/lincrna/cabili/multi_exons.fasta"
# families_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt"


motif_file = "source_data/motif_sets/int3.txt"
motifs = sequo.read_motifs(motif_file)
stop_motifs = [i for i in motifs if len(re.findall("(?=TAA|TAG|TGA)", i)) > 0]
non_stop_motifs = [i for i in motifs if i not in stop_motifs]


print(len(stop_motifs), len(non_stop_motifs), np.divide(len(non_stop_motifs), len(stop_motifs)))
