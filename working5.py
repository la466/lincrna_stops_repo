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


ese_file = "source_data/motif_sets/int3.txt"
eses = sequo.read_motifs(ese_file)

output_directory = "clean_run/tests/ese_densities"
gen.create_output_directories(output_directory)

fileset = {
    "fileset_prefix": "pc",
    "exons_file": "clean_run/genome_sequences/human/human.cds.clean_coding_exons.fasta",
    "introns_file": "clean_run/genome_sequences/human/human.clean_introns.fasta",
    "families_file": "clean_run/genome_sequences/human/human.cds.families.bed" ,
}



exon_names, exon_seqs = gen.read_fasta(fileset["exons_file"])
exon_list = collections.defaultdict(lambda: [])
[exon_list[name.split(".")[0]].append(exon_seqs[i]) for i, name in enumerate(exon_names[:100])]

for id in exon_list:

    for seq in exon_list[id]:
        overlap_motifs = sequo.return_overlap_motifs(seq, eses)
        print(overlap_motifs)
