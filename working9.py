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


motifs = "source_data/motif_sets/int3.txt"
# motifs = "source_data/motif_sets/RESCUE.txt"
# seqs_file = "clean_run/lincrna/lincRNA.multi_exon.exons.fasta"
# seqs_file = "clean_run/genome_sequences/lincrna/cabili/multi_exons.fasta"

seqs_file = "clean_run/genome_sequences/human/human.cds.clean_coding_exons.fasta"
families_file = "clean_run/genome_sequences/human/human.cds.families.bed"

motifs = sequo.read_motifs(motifs)

names, seqs = gen.read_fasta(seqs_file)
seq_list = collections.defaultdict(lambda: [])
[seq_list[name.split(".")[0]].append(seqs[i]) for i, name in enumerate(names[:3000]) if len(seqs[i]) > 211]

# seq_list = sequo.pick_random_family_member(families_file, seq_list)


def get_density(sequences):
    hits = sum([len(sequo.sequence_overlap_indicies(i, motifs)) for i in sequences])
    lengths = sum([len(i) for i in sequences])
    density = np.divide(hits, lengths)
    return density

flank_densities = []
core_densities = []

flanks = []
cores = []
for id in seq_list:
    exons = seq_list[id]
    for exon in exons:
        flanks.append(exon[2:69])
        flanks.append(exon[-69:-2])
        midpoint = int(len(exon) // 2)
        cores.append(exon[midpoint-33:midpoint+34])

def get_real(flanks, cores, flank_excess):

    flank_motifs = []
    [flank_motifs.extend(sequo.return_overlap_motifs(flank, motifs)) for flank in flanks]
    flank_stops = []
    [flank_stops.extend(sequo.return_overlap_motifs(flank, stops)) for flank in flank_motifs]
    flank_stop_density = np.divide(len(flank_stops), len(flank_motifs))

    core_motifs = []
    [core_motifs.extend(sequo.return_overlap_motifs(core, motifs)) for core in cores]
    core_stops = []
    [core_stops.extend(sequo.return_overlap_motifs(core, stops)) for core in core_motifs]
    core_stop_density = np.divide(len(core_stops), len(core_motifs))

    flank_density = np.divide(sum([len(i) for i in flank_motifs]), sum([len(i) for i in flanks]))
    core_density = np.divide(sum([len(i) for i in core_motifs]), sum([len(i) for i in cores]))
    print(flank_density)
    print(core_density)
    print(flank_stop_density)
    print(core_stop_density)

    # if flank_excess > 0:
    #     factor = 1+flank_excess
    # else:
    #     factor = 1- flank_excess
    #
    # expected_flank = factor*core_stop_density
    #
    # print(flank_stop_density)
    # print(core_stop_density)
    #
    # # print(core_stop_density)
    # # print(flank_stop_density, expected_flank)
    # diff = flank_stop_density - expected_flank
    # return diff*100

def get_sims(iterations, flanks, cores, flank_excess):
    outputs = []
    if len(iterations):
        np.random.seed()
        for i, it in enumerate(iterations):
            gen.print_parallel_status(i, iterations)

            flank_motifs = []
            for seq in flanks:
                overlaps = sequo.sequence_overlap_indicies(seq, motifs)
                chunked_overlaps = sequo.chunk_overlaps(overlaps)
                new_overlaps = []
                for chunk in chunked_overlaps:
                    choice = np.random.choice(list(range(len(seq)-len(chunk))))
                    new_overlaps.append(list(range(choice, choice+len(chunk))))
                overlap_motifs = ["".join([seq[i] for i in chunk]) for chunk in new_overlaps]
                flank_motifs.extend(overlap_motifs)
            core_motifs = []
            for seq in cores:
                overlaps = sequo.sequence_overlap_indicies(seq, motifs)
                chunked_overlaps = sequo.chunk_overlaps(overlaps)
                new_overlaps = []
                for chunk in chunked_overlaps:
                    choice = np.random.choice(list(range(len(seq)-len(chunk))))
                    new_overlaps.append(list(range(choice, choice+len(chunk))))
                overlap_motifs = ["".join([seq[i] for i in chunk]) for chunk in new_overlaps]
                core_motifs.extend(overlap_motifs)

            flank_stops = []
            [flank_stops.extend(sequo.return_overlap_motifs(flank, stops)) for flank in flank_motifs]
            flank_stop_density = np.divide(len(flank_stops), len(flank_motifs))
            core_stops = []
            [core_stops.extend(sequo.return_overlap_motifs(core, stops)) for core in core_motifs]
            core_stop_density = np.divide(len(core_stops), len(core_motifs))

            flank_stop_density = np.divide(len(flank_stops), len(flank_motifs))
            core_stop_density = np.divide(len(core_stops), len(core_motifs))
            if flank_excess > 0:
                factor = 1+flank_excess
            else:
                factor = 1- flank_excess
            expected_flank = factor*core_stop_density
            diff = (flank_stop_density - expected_flank)*100
            outputs.append(diff)


    return(outputs)


flank_density = seqo.calc_motif_density(flanks, motifs)
core_density = seqo.calc_motif_density(cores, motifs)
flank_excess = np.divide(flank_density - core_density, core_density)


# outputs = soc.run_simulation_function(list(range(5)), [flanks, cores, flank_excess], get_sims, sim_run = False)
real = get_real(flanks, cores, flank_excess)

# less = [i for i in outputs if i <= real]

print(real)
# print(len(less))
