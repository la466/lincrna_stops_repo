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
introns_fasta = "clean_run/genome_sequences/human/human.clean_introns.fasta"
cds_file = "clean_run/genome_sequences/human/human.cds.multi_exons.fasta"
families_file = "clean_run/genome_sequences/human/human.cds.families.bed"
# file = "clean_run/genome_sequences/lincrna/cabili/multi_exons.fasta"
# families_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt"


motif_set = "int3"
motif_file = "source_data/motif_sets/{0}.txt".format(motif_set)
sim_path = "clean_run/dinucleotide_controls/{0}_dinucleotide_controls_matched_stops".format(motif_set)
sims = gen.get_filepaths(sim_path)[:1000]


# get the exon fasta
exon_names, exon_seqs = gen.read_fasta(file)
# get the introns fasta
intron_names, intron_seqs = gen.read_fasta(introns_fasta)

flanks = True
restrict_size = True

motifs = sequo.read_motifs(motif_file)

# if limited to flanks, extract the flanks
if flanks:
    # get all the exon flanks
    exon_list = collections.defaultdict(lambda: collections.defaultdict())
    for i, name in enumerate(exon_names):
        if len(exon_seqs[i]) > 211 and "N" not in exon_seqs[i]:
            exon_list[name.split(".")[0]][int(name.split(".")[-1].split("(")[0])] = exon_seqs[i]
    exon_list = {i: exon_list[i] for i in exon_list}


    intron_list = collections.defaultdict(lambda: [])

    # get only those that flank an included exon
    for i, name in enumerate(intron_names):
        intron_id = name.split(".")[0]
        if intron_id in exon_list:
            intron_flanking_exons = [int(i) for i in name.split(".")[-1].split("(")[0].split("-")]
            if intron_flanking_exons[0] in exon_list[intron_id] or intron_flanking_exons[1] in exon_list[intron_id]:
                intron_list[name.split(".")[0]].append(intron_seqs[i])
    intron_list = {i: intron_list[i] for i in intron_list}

    temp_exon_list = collections.defaultdict(lambda: [])
    for id in exon_list:
        flanks = []
        for exon_id in exon_list[id]:
            flanks.append(exon_list[id][exon_id][2:69])
            flanks.append(exon_list[id][exon_id][-69:-2])
        temp_exon_list[id] = flanks
    exon_list = temp_exon_list
    exon_list = {i: exon_list[i] for i in exon_list}

else:
    # if not flanks or is non coding exons
    exon_list = collections.defaultdict(lambda: collections.defaultdict())
    for i, name in enumerate(exon_names):
        if "N" not in exon_seqs[i]:
            if restrict_size and len(exon_seqs[i]) > 211:
                exon_list[name.split(".")[0]][int(name.split(".")[-1].split("(")[0])] = exon_seqs[i]
        else:
            exon_list[name.split(".")[0]][int(name.split(".")[-1].split("(")[0])] = exon_seqs[i]
    exon_list = {i: exon_list[i] for i in exon_list}

    intron_list = collections.defaultdict(lambda: [])
    [intron_list[name.split(".")[0]].append(intron_seqs[i]) for i, name in enumerate(intron_names) if name.split(".")[0] in exon_list]
    intron_list = {i: intron_list[i] for i in intron_list}

    temp_exon_list = collections.defaultdict(lambda: [])
    for id in exon_list:
        exons = []
        for exon_id in exon_list[id]:
            exons.append(exon_list[id][exon_id])
        temp_exon_list[id] = exons
    exon_list = temp_exon_list
    exon_list = {i: exon_list[i] for i in exon_list}


# print(exon_list)

intron_sizes = {i: np.median([len(intron) for intron in intron_list[i]]) for i in intron_list}

exon_list = sequo.pick_random_family_member(families_file, exon_list)
intron_list = {i: intron_list[i] for i in exon_list if i in intron_list}

with open("temp_files/test_intron.csv", "w") as outfile:
    outfile.write("id,overlap_count,overlap_prop,intron_size,density\n")
    for id in exon_list:
        if id in intron_list:
            exons = exon_list[id]
            introns = intron_list[id]

            motif_overlaps = seqo.calc_motif_density(exons, motifs)
            motif_hits = gen.flatten([sequo.return_overlap_motifs(exon, motifs) for exon in exons])
            overlap_motifs = [i for i in motif_hits if len(i) == 6]
            overlap_prop = np.divide(len(overlap_motifs), len(motif_hits))

            outfile.write("{0},{1},{2},{3},{4}\n".format(id, len(overlap_motifs), overlap_prop, intron_sizes[id], motif_overlaps))
