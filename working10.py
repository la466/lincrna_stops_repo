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

nts = ["A", "C", "G", "T"]

motifs = "source_data/motif_sets/int3.txt"
# motifs = "source_data/motif_sets/RESCUE.txt"
# seqs_file = "clean_run/lincrna/lincRNA.multi_exon.exons.fasta"
# families_file = "clean_run/lincrna/lincRNA.multi_exon_families.bed"
# seqs_file = "clean_run/genome_sequences/lincrna/cabili/multi_exons.fasta"
# families_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt"

# seqs_file = "clean_run/genome_sequences/human/human.cds.clean_coding_exons.fasta"
families_file = "clean_run/genome_sequences/human/human.cds.families.bed"

seqs_file = "clean_run/genome_sequences/human.macaque.alignments.fasta"
# seqs_file = "source_data/cabili_clean_alignments.fasta"

motifs = sequo.read_motifs(motifs)

names, seqs = gen.read_fasta(seqs_file)
seq_list = collections.defaultdict(lambda: [])
[seq_list[name.split(".")[0]].append(seqs[i].split(",")) for i, name in enumerate(names)]

seq_list = sequo.pick_random_family_member(families_file, seq_list)

sub_ese = 0
ese = 0
sub_non_ese = 0
non_ese = 0

ese_stop_hits = 0
ese_stop_subs = 0
non_ese_stop_hits = 0
non_ese_stop_subs = 0

for id in seq_list:
    sequences = seq_list[id][0]
    human = sequences[0]
    mac = sequences[1]

    # print(seqs[0])
    hits = sequo.sequence_overlap_indicies(human, motifs)
    stop_hits = sequo.sequence_overlap_indicies(human, stops)
    # print(hits)
    for i in range(len(human)):
        if human[i] in nts and mac[i] in nts:
            if i in hits:
                ese += 1
                if human[i] != mac[i]:
                    sub_ese += 1
                if i in stop_hits:
                    ese_stop_hits += 1
                    if human[i] != mac[i]:
                        ese_stop_subs += 1
            else:
                non_ese += 1
                if human[i] != mac[i]:
                    sub_non_ese += 1
                if i in stop_hits:
                    non_ese_stop_hits += 1
                    if human[i] != mac[i]:
                        non_ese_stop_subs += 1



print(np.divide(sub_ese, ese))
print(np.divide(sub_non_ese, non_ese))


print(np.divide(ese_stop_subs, ese_stop_hits))
print(np.divide(non_ese_stop_subs, non_ese_stop_hits))
