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


motif_set = "RESCUE"
motif_file = "source_data/motif_sets/{0}.txt".format(motif_set)
sim_path = "clean_run/dinucleotide_controls/{0}_dinucleotide_controls_matched_stops".format(motif_set)
sims = gen.get_filepaths(sim_path)[:10]


names, exons = gen.read_fasta(file)
exon_list = collections.defaultdict(lambda: [])
[exon_list[name.split(".")[0]].append(exons[i]) for i, name in enumerate(names)]

exon_list = sequo.pick_random_family_member(families_file, exon_list)
exons = []
[exons.extend(exon_list[i]) for i in exon_list]


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
real_all = np.mean([overlaps[i] for i in overlaps])
real_stop, real_non_stop, real_diff = calc_stop_overlaps(overlaps, stop_motifs, non_stop_motifs)


def get_overlap_stats(filepaths):
    outputs = []
    if len(filepaths):
        for i, filepath in enumerate(filepaths):
            gen.print_parallel_status(i, filepaths)
            motifs, stop_motifs, non_stop_motifs = get_motifs(filepath)
            overlaps = calc_overlaps(motifs)
            sim_all = np.mean([overlaps[i] for i in overlaps])
            sim_stop, sim_non_stop, sim_diff = calc_stop_overlaps(overlaps, stop_motifs, non_stop_motifs)
            outputs.append([sim_stop, sim_non_stop, sim_diff, sim_all])
    return outputs


def calc_diff(motif_file):

    motifs, stop_motifs, non_stop_motifs = get_motifs(motif_file)
    overlaps = calc_overlaps(motifs)
    stop, non_stop, diff = calc_stop_overlaps(overlaps, stop_motifs, non_stop_motifs)

    # print(stop, non_stop)
    percentage_diff = np.divide(stop - non_stop, non_stop)


    # return diff


    all_motifs = []
    for exon in exons:
        sequence_motifs = sequo.return_overlap_motifs(exon, motifs)
        all_motifs.extend(sequence_motifs)

    overlap_motifs = [i for i in all_motifs if len(i) > 6]
    non_overlap_motifs = [i for i in all_motifs if len(i) == 6]

    stop_overlaps = seqo.calc_motif_density(overlap_motifs, stops)
    stop_non_overlaps = seqo.calc_motif_density(non_overlap_motifs, stops)


    overlap_diff = np.divide(stop_non_overlaps - stop_overlaps, stop_overlaps)
    pd = percentage_diff*100
    od = overlap_diff*100
    print(stop_overlaps, stop_non_overlaps, pd, od, od - pd)


real = calc_diff(motif_file)
diffs = [calc_diff(i) for i in sims]

print(len([i for i in diffs if i >= real]))
