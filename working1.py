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




motif_file = "source_data/motif_sets/int3.txt"
motifs = sorted(sequo.read_motifs(motif_file))

cdir = "clean_run/dinucleotide_controls/int3_dinucleotide_controls_matched_stops"
for file in os.listdir(cdir)[:10]:
    path = "{0}/{1}".format(cdir, file)
    motifs = sequo.read_motifs(path)
    test = calc_subs(motifs)
    print(sum(test))
