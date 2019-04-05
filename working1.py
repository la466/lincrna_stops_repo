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


motif_set = "int3"
motif_file = "source_data/motif_sets/{0}.txt".format(motif_set)
sim_path = "clean_run/dinucleotide_controls/{0}_dinucleotide_controls_matched_stops".format(motif_set)
sims = gen.get_filepaths(sim_path)[:1000]


def calc_overlaps(motifs):
    overlaps = collections.defaultdict(lambda: [])
    for motif in sorted(motifs):
        for overlap_motif in motifs:

            if overlap_motif[-2:] == motif[:2]:
                overlaps[motif].append(overlap_motif)
            if overlap_motif[-3:] == motif[:3]:
                overlaps[motif].append(overlap_motif)
            if overlap_motif[-4:] == motif[:4]:
                overlaps[motif].append(overlap_motif)
            if overlap_motif[-5:] == motif[:5]:
                overlaps[motif].append(overlap_motif)
            if motif[1:] == overlap_motif[:5]:
                overlaps[motif].append(overlap_motif)
            if motif[2:] == overlap_motif[:4]:
                overlaps[motif].append(overlap_motif)
            if motif[3:] == overlap_motif[:3]:
                overlaps[motif].append(overlap_motif)
            if motif[4:] == overlap_motif[:2]:
                overlaps[motif].append(overlap_motif)

    overlaps = {i: len(list(set(overlaps[i]))) for i in overlaps}
    return overlaps

def randomise_overlaps(overlaps):
    ids = list(overlaps)
    counts = [overlaps[i] for i in overlaps]
    np.random.shuffle(counts)
    overlaps = dict(zip(ids, counts))
    return overlaps

def calc_stop_overlaps(overlaps, stop_motifs, non_stop_motifs):
    stop_overlaps = {i: overlaps[i] for i in overlaps if i in stop_motifs}
    non_stop_overlaps = {i: overlaps[i] for i in overlaps if i in non_stop_motifs}

    stop_mean = np.mean([stop_overlaps[i] for i in stop_overlaps])
    non_stop_mean = np.mean([non_stop_overlaps[i] for i in non_stop_overlaps])
    diff = non_stop_mean - stop_mean

    return stop_mean, non_stop_mean, diff

def get_motifs(motif_file):
    motifs = sequo.read_motifs(motif_file)
    stop_motifs = [i for i in motifs if len(re.findall("(?=(TAA|TAG|TGA))", i)) > 0]
    non_stop_motifs = [i for i in motifs if i not in stop_motifs]
    return motifs, stop_motifs, non_stop_motifs

motifs, stop_motifs, non_stop_motifs = get_motifs(motif_file)
overlaps = calc_overlaps(motifs)
real_stop, real_non_stop, real_diff = calc_stop_overlaps(overlaps, stop_motifs, non_stop_motifs)


def get_overlap_stats(filepaths):
    outputs = []
    if len(filepaths):
        for i, filepath in enumerate(filepaths):
            gen.print_parallel_status(i, filepaths)
            motifs, stop_motifs, non_stop_motifs = get_motifs(filepath)
            overlaps = calc_overlaps(motifs)
            sim_stop, sim_non_stop, sim_diff = calc_stop_overlaps(overlaps, stop_motifs, non_stop_motifs)
            outputs.append([sim_stop, sim_non_stop, sim_diff])
    return outputs



def calc_ps(real, sim_outputs):
    real_stop = real[0]
    real_non_stop = real[1]
    real_diff = real[2]
    sim_stop_outputs = [output[0] for output in sim_outputs]
    sim_non_stop_outputs = [output[1] for output in sim_outputs]
    sim_diff_outputs = [output[2] for output in sim_outputs]

    stop_p = np.divide(len([i for i in sim_stop_outputs if i <= real_stop]) + 1, len(sim_stop_outputs) + 1)
    non_stop_p = np.divide(len([i for i in sim_non_stop_outputs if i >= real_non_stop]) + 1, len(sim_non_stop_outputs) + 1)
    diff_p = np.divide(len([i for i in sim_diff_outputs if i >= real_diff]) + 1, len(sim_diff_outputs) + 1)

    nd = np.divide(real_stop - np.mean(sim_stop_outputs), np.mean(sim_stop_outputs))
    non_nd = np.divide(real_non_stop - np.mean(sim_non_stop_outputs), np.mean(sim_non_stop_outputs))

    print(np.mean(sim_diff_outputs))
    print(stop_p, non_stop_p, diff_p)
    print(nd, non_nd)


dint_outputs = soc.run_simulation_function(sims, [], get_overlap_stats, sim_run = False)

print(real_stop, real_non_stop, real_diff)
calc_ps([real_stop, real_non_stop, real_diff], dint_outputs)
