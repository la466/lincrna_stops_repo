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


motif_file = "source_data/motif_sets/int3.txt"
simdir = "clean_run/dinucleotide_controls/int3_dinucleotide_controls_matched_stops"
sims = gen.get_filepaths(simdir)[:30]



names, seqs = gen.read_fasta(file)
seq_list = collections.defaultdict(lambda: [])
[seq_list[name.split(".")[0]].append(seqs[i]) for i, name in enumerate(names)]
seq_list = sequo.pick_random_family_member(families_file, seq_list)
exons = []
[exons.extend(seq_list[i]) for i in seq_list]


cds_names, cds_seqs = gen.read_fasta(cds_file)
cds_list = {name.split(".")[0]: cds_seqs[i] for i, name in enumerate(cds_names)}
# cds_list = sequo.pick_random_family_member(families_file, cds_list)
#
# cds = [cds_list[i] for i in cds_list]


def calc_density(motifs, exons):
    outputs = {}
    if len(motifs):
        for motif in motifs:
            outputs[motif] = seqo.calc_motif_density(exons, [motif])
    return outputs
def calc_densities(motif_file):
    motifs = sequo.read_motifs(motif_file)
    outputs = soc.run_simulation_function(motifs, [exons], calc_density, sim_run = False, workers = 18)

    stop_motifs = [i for i in motifs if len(re.findall("(?=TAA|TAG|TGA)", i)) > 0]
    non_stop_motifs = [i for i in motifs if i not in stop_motifs]
    stop_densities = [outputs[i] for i in stop_motifs]
    non_stop_densities = [outputs[i] for i in non_stop_motifs]
    non_stop_densities = [np.divide(i, 3)*2 for i in non_stop_densities]

    print(np.mean(stop_densities))
    print(np.mean(non_stop_densities))

    file = "temp_files/_output1.csv"

    d = {'stop': np.array(stop_densities), 'non_stop': np.array(non_stop_densities)}
    data = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in d.items() ]))

    # data = pd.DataFrame({'stop': stop_densities, 'non_stop': non_stop_densities})
    # data = data.transpose()
    data.to_csv(path_or_buf = file)

    return [np.mean(stop_densities), np.mean(non_stop_densities)]
def calc_all_density(motif_file):
    motifs = sequo.read_motifs(motif_file)
    stop_motifs = [i for i in motifs if len(re.findall("(?=TAA|TAG|TGA)", i)) > 0]
    non_stop_motifs = [i for i in motifs if i not in stop_motifs]

    stop_density = seqo.calc_motif_density(exons, stop_motifs)
    non_stop_density = seqo.calc_motif_density(exons, non_stop_motifs)

    # print(stop_density, non_stop_density, len(stop_motifs), seqo.calc_motif_density(motifs, stops))
    print(stop_density, non_stop_density*np.divide(2,3)*np.divide(len(stop_motifs), len(non_stop_motifs)))

    stop_hits = seqo.calc_motif_counts(exons, stop_motifs)
    non_stop_hits = seqo.calc_motif_counts(exons, non_stop_motifs)
    total_hits = stop_hits + non_stop_hits

    # normalised_non_stop_hits = non_stop_hits*np.divide(2,3)*np.divide(len(stop_motifs), len(non_stop_motifs))
    expected_stop_hits = total_hits*np.divide(len(stop_motifs), len(motifs))
    expected_non_stop_hits = total_hits*np.divide(len(non_stop_motifs), len(motifs))

    print(total_hits, expected_stop_hits, stop_hits, expected_non_stop_hits, non_stop_hits)

    chi = scipy.stats.chisquare([stop_hits, non_stop_hits], f_exp = [expected_stop_hits, expected_non_stop_hits])
    print(chi)
# real = calc_all_density(motif_file)

motifs = sequo.read_motifs(motif_file)
stop_motifs = [i for i in motifs if len(re.findall("(?=TAA|TAG|TGA)", i)) > 0]
non_stop_motifs = [i for i in motifs if i not in stop_motifs]

print(stop_motifs)
#
#

def calc_hits(motif_files):
    outputs = []

    if len(motif_files):
        for file in motif_files:
            motifs = sequo.read_motifs(file)
            stop_motifs = [i for i in motifs if len(re.findall("(?=TAA|TAG|TGA)", i)) > 0]
            non_stop_motifs = [i for i in motifs if i not in stop_motifs]

            s = []
            ns = []

            for frame in [0, 1, 2]:
                frame_hits = seqo.calc_motif_count_frame(exons, motifs, frame = frame)
                stop_hits = [i for i in frame_hits if len(re.findall("(?=TAA|TAG|TGA)", i)) > 0]
                non_stop_hits = [i for i in frame_hits if len(re.findall("(?=TAA|TAG|TGA)", i)) == 0]
                s.append(len(stop_hits))
                ns.append(len(non_stop_hits))

            s_adj = copy.deepcopy(s)
            s_adj[0] = np.divide(s_adj[0], len(stop_motifs))*np.divide(len(stop_motifs), 5)
            s_adj[1] = np.divide(s_adj[1], len(stop_motifs))*np.divide(len(stop_motifs), 6)
            s_adj[2] = np.divide(s_adj[2], len(stop_motifs))*np.divide(len(stop_motifs), 7)
            n_adj = [np.divide(i, len(non_stop_motifs)) for i in ns]

            total_s = sum(s_adj)
            total_n = sum(n_adj)

            outputs.append([total_s, total_n])
    return outputs
    # print(total_n)

    # print(s, sum(s))
    # print(ns, sum(ns))
    # print(s_adj)
    # print(n_adj)


real = calc_hits([motif_file])[0]
print(real)
outputs = soc.run_simulation_function(sims, [], calc_hits, sim_run = False)


nd_s = np.divide(real[0] - np.mean([i[0] for i in outputs]), np.mean([i[0] for i in outputs]))
nd_n = np.divide(real[1] - np.mean([i[1] for i in outputs]), np.mean([i[1] for i in outputs]))

print(nd_s)
print(nd_n)
