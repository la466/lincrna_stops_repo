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


file = "source_data/motif_sets/iss.txt"
motifs = sequo.read_motifs(file)

stop_motifs = [i for i in motifs if len(re.findall("(?=(TAA|TAG|TGA))", i)) > 0]

print(np.divide(len(stop_motifs), len(motifs)))
