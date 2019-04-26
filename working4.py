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


ese_set = "combined_eses"
motif_file = "source_data/motif_sets/{0}.txt".format(ese_set)
controls_dir = "clean_run/dinucleotide_controls/{0}_dinucleotide_controls_matched_stops".format(ese_set)
seqs_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_transcript_sequences.fasta"
families_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt"

names, seqs = gen.read_fasta(seqs_file)
all = dict(zip([name.split(".")[0] for name in names], seqs))
all = sequo.pick_random_family_member(families_file, all)
seqs = [all[i] for i in all]

print(len(seqs))

# seqs = seqs[:500]


# print(seqs)

sims = gen.get_filepaths(controls_dir)[:15]

def return_overlap_motifs(seq, overlaps, inverse = None):
    if inverse:
        overlaps = [i for i in list(range(len(seq))) if i not in overlaps]
    else:
        overlaps = [i for i in list(range(len(seq))) if i in overlaps]
    chunked_overlaps = sequo.chunk_overlaps(overlaps)
    overlap_motifs = ["".join([seq[i] for i in chunk]) for chunk in chunked_overlaps]
    return overlap_motifs

def calc_density(motif_set):
    motifs = sequo.read_motifs(motif_set)
    motif_gc = seqo.calc_seq_gc("".join(motifs))
    non_hits = []
    hits = []
    for seq in seqs:
        overlaps = sequo.sequence_overlap_indicies(seq, motifs)
        hits.extend(return_overlap_motifs(seq, overlaps))
        non_hits.extend(return_overlap_motifs(seq, overlaps, inverse = True))
    hit_density = seqo.calc_motif_density(hits, stops)
    non_hit_density = seqo.calc_motif_density(non_hits, stops)
    hit_gc = seqo.calc_seq_gc("".join(non_hits))
    print(hit_density, non_hit_density, np.divide(hit_density, non_hit_density), motif_gc, hit_gc)

calc_density(motif_file)
[calc_density(i) for i in sims]

#
# if simulation_id != "real":
#     hit_lengths = [len(i) for i in hits]
#     all_nts = list("".join(hits))
#     np.random.shuffle(all_nts)
#     nt_index = 0
#     new_hits = []
#     for length in hit_lengths:
#         new_hits.append("".join(all_nts[nt_index:nt_index+length]))
#         nt_index += length
#     hits = new_hits
#
# # [test_kept.extend(sequo.return_overlap_motifs(i, motifs, inverse = True)) for i in test_seqs]
# density = seqo.calc_motif_density(hits, stops)
