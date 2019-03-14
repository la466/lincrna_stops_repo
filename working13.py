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
import scipy.stats


dir = "clean_run/dinucleotide_controls/int3_dinucleotide_controls_matched_stops"

for file in os.listdir(dir):
    path = dir + "/" + file
    motifs = sequo.read_motifs(path)
    print(len([i for i in motifs if len(re.findall("(?=(TAA|TAG|TGA))", i)) == 0]))
