import generic as gen
import seq_ops as seqo
import sequence_ops as sequo
import sim_ops_containers as soc
import conservation as cons
import file_ops as fo
import sim_ops as simo
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

exon_bed = "source_data/hsAllCompmerge.cage+polyASupported.bed"
single_exon_bed = "clean_run/genome_sequences/lincrna/GENCODE_CLS/single_exons.bed"
genome_fasta = "../source_data/genomes/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
multi_exon_fasta = "clean_run/genome_sequences/lincrna/GENCODE_CLS/hsAllCompmerge.cage+polyASupported.bed.filtered.exons.fasta"
multi_exon_families_file = "clean_run/genome_sequences/lincrna/GENCODE_CLS/families.txt"
single_exon_fasta = "clean_run/genome_sequences/lincrna/GENCODE_CLS/single_exons.fasta"
single_exon_families_file = "clean_run/genome_sequences/lincrna/GENCODE_CLS/single_exons_families.txt"
multi_exon_transcripts = "clean_run/genome_sequences/lincrna/GENCODE_CLS/hsAllCompmerge.cage+polyASupported.full_transcripts.fasta"

motif_file = "source_data/motif_sets/int3.txt"
motifs = sequo.read_motifs(motif_file)


np.random.seed()

def get_seq_list(fasta_file, with_chr = False, full_seq_list = None):
    names, seqs = gen.read_fasta(fasta_file)
    seq_list = collections.defaultdict(lambda: [])
    for i, name in enumerate(names):
        id = name.split("(")[0]
        if with_chr:
            id = ".".join(id.split(".")[:-1])
        if full_seq_list:
            if id in full_seq_list:
                seq_list[id].append(seqs[i])
        else:
            seq_list[id].append(seqs[i])
    seq_list = {id: seq_list[id] for i, id in enumerate(seq_list)}
    return seq_list

def randomise(seq):
    nts = list(seq)
    np.random.shuffle(nts)
    return "".join(nts)

def randomise_densities(ids, seq_list, iterations):

    random_densities = collections.defaultdict()
    if len(ids):

        for i, id in enumerate(ids):
            exons = seq_list[id]
            sim_densities = []

            for j, exon in enumerate(exons):
                sim_densities = []
                for iteration in range(iterations):
                    sim_exon = randomise(exon)
                    sim_density = seqo.calc_motif_density([sim_exon], stops)
                    sim_densities.append(sim_density)
                random_densities['{0}.{1}'.format(id, j)] = sim_densities

    return random_densities

def process_densities(sim_outputs):
    random_densities = collections.defaultdict(lambda: [])
    for output in sim_outputs:
        result = output.get()
        for id in result:
            random_densities[id].extend(result[id])

    random_densities = {i: random_densities[i] for i in random_densities}
    return random_densities

def calc_nds(densities, randomise_densities):
    nds = collections.defaultdict(lambda: [])
    for id in randomise_densities:
        for i, exon_density in enumerate(densities[id]):
            nd = np.divide(exon_density - np.mean(randomise_densities[id][i]), np.mean(randomise_densities[id][i]))
            nds[id].append(nd)
    return nds

def calc_values(seq_list):

    densities = collections.defaultdict(lambda: [])
    gcs = collections.defaultdict(lambda: [])
    ese = collections.defaultdict(lambda: [])

    for id in seq_list:
        for i, exon in enumerate(seq_list[id]):
            density = seqo.calc_motif_density([exon], stops)
            densities['{0}.{1}'.format(id, i)].append(density)
            gcs['{0}.{1}'.format(id, i)].append(seqo.calc_gc_seqs_combined([exon]))
            ese['{0}.{1}'.format(id, i)].append(seqo.calc_motif_density([exon], motifs))


    ids = list(seq_list)
    its = 1000
    sim_outputs = gen.run_in_parallel(ids, ["foo", seq_list, its], randomise_densities)

    randomised_densities = process_densities(sim_outputs)
    nds = calc_nds(densities, randomised_densities)

    return densities, nds, gcs, ese


motif_file = "source_data/motif_sets/int3.txt"

names, seqs = gen.read_fasta(multi_exon_transcripts)
full_seq_list = {".".join(name.split(":")): seqs[i] for i, name in enumerate(names)}

single_exon_seq_list = get_seq_list(single_exon_fasta)
single_exon_seq_list = sequo.pick_random_family_member(single_exon_families_file, single_exon_seq_list)


multi_exon_seq_list = get_seq_list(multi_exon_fasta, full_seq_list = full_seq_list, with_chr = True)
multi_exon_seq_list = sequo.pick_random_family_member(multi_exon_families_file, multi_exon_seq_list)

densities_s, nds_s, gcs_s, ese_s = calc_values(single_exon_seq_list)
densities_m, nds_m, gcs_m, ese_m = calc_values(multi_exon_seq_list)

exon_count_s = {i: len(single_exon_seq_list[i]) for i in single_exon_seq_list}
exon_count_m = {i: len(multi_exon_seq_list[i]) for i in multi_exon_seq_list}


with open("clean_run/gc_output_single.csv", "w") as outfile:
    outfile.write("id,gc,density,fe,ese_density\n")
    for id in nds_s:
        outfile.write("{0},{1},{2},{3},{4}\n".format(id, gcs_s[id][0], densities_s[id][0], nds_s[id][0], ese_s[id][0]))

with open("clean_run/gc_output_multi.csv", "w") as outfile:
    outfile.write("id,gc,density,fe,ese_density\n")
    for id in nds_m:
        outfile.write("{0},{1},{2},{3},{4}\n".format(id, gcs_m[id][0], densities_m[id][0], nds_m[id][0], ese_m[id][0]))
