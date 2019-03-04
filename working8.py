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


def clean_sequences(alignment_file, sequence_file, output_file, indel_threshold = 0.15):

    names, seqs = gen.read_fasta(sequence_file)
    seq_list = {name: seqs[i] for i, name in enumerate(names)}

    entries = gen.read_many_fields(alignment_file, "\t")

    kept_list = {}

    for i in range(0, len(entries), 5):
        alignment_entry = entries[i:i+4]

        id = alignment_entry[0][0].split(".")[-1]
        seq1 = alignment_entry[1][0]
        seq2 = alignment_entry[3][0]

        if id in seq_list:
            if len(seq1) == len(seq_list[id]):
                indel_prop1 = np.divide(seq1.count("-"), len(seq1))
                indel_prop2 = np.divide(seq2.count("-"), len(seq2))
                if indel_prop1 < indel_threshold and indel_prop2 < indel_threshold:
                    kept_list[id] = [seq1, seq2]


    with open(output_file, "w") as outfile:
        for id in sorted(kept_list):
            outfile.write(">{0}\n{1},{2}\n".format(id, kept_list[id][0].upper(), kept_list[id][1].upper()))


file = "source_data/cabili_alignments.fasta"
sequences = "clean_run/genome_sequences/lincrna/cabili/transcript_sequences.fasta"
clean_file = "source_data/cabili_clean_alignments.fasta"
ese_file = "source_data/motif_sets/RESCUE.txt"

eses = sequo.read_motifs(ese_file)

# clean_sequences(file, sequences, clean_file)

names, alignments = gen.read_fasta(clean_file)
alignments = {name: alignments[i].split(",") for i, name in enumerate(names)}

sub_rate = []
sub_rate_non = []
for id in names[:1]:
    alignment_seqs = alignments[id]
    seq1 = alignment_seqs[0]
    seq2 = alignment_seqs[1]

    overlaps = sequo.sequence_overlap_indicies(seq1, eses)
    non_overlaps = [i for i in range(len(seq1)) if i not in overlaps]
    seq1_motifs = [seq1[i] for i in overlaps]
    seq2_motifs = [seq2[i] for i in overlaps]
    seq1_non = [seq1[i] for i in non_overlaps]
    seq2_non = [seq2[i] for i in non_overlaps]

    print(overlaps)
    print(non_overlaps)

    subs = 0
    for i in range(len(seq1_motifs)):
        if seq1_motifs[i] != seq2_motifs[i]:
            subs += 1

    subs_non = 0
    for i in range(len(seq1_non)):
        if seq1_non[i] != seq2_non[i]:
            subs_non += 1

    nts = len(seq1_motifs)
    sub_rate.append(np.divide(subs, nts))
    sub_rate_non.append(np.divide(subs_non, nts))

print(np.nanmedian(sub_rate))
print(np.nanmedian(sub_rate_non))
