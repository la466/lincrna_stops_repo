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




file = "clean_run/tests/lincrna/substitution_rates/lincrna_int3_substitution_rates_motif.csv"

data = pd.read_csv(file)

data.rename(columns={'relative_ese_rate':'relative_ese_diff'}, inplace=True)
data.rename(columns={'relative_non_ese_rate':'relative_non_ese_diff'}, inplace=True)


data.to_csv(path_or_buf = "clean_run/tests/lincrna/substitution_rates/lincrna_int3_substitution_rates_motif1.csv", index = False)
