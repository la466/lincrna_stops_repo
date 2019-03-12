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
lseqs_file = "clean_run/genome_sequences/lincrna/cabili/multi_exons.fasta"
lseqs_families = "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt"

seqs_file = "clean_run/genome_sequences/human/human.cds.clean_coding_exons.fasta"
families_file = "clean_run/genome_sequences/human/human.cds.families.bed"

motifs = sequo.read_motifs(motifs)

names, seqs = gen.read_fasta(seqs_file)
seq_list = collections.defaultdict(lambda: [])
[seq_list[name.split(".")[0]].append(seqs[i]) for i, name in enumerate(names) if len(seqs[i]) > 211]

seq_list = sequo.pick_random_family_member(families_file, seq_list)

seqs_list = []
[seqs_list.extend(seq_list[i]) for i in seq_list]

# seqs_list_flanks = []
# [seqs_list_flanks.extend([i[2:69], i[-69:-2]]) for i in seqs_list]

lnames, lseqs = gen.read_fasta(lseqs_file)
lseq_list = collections.defaultdict(lambda: [])
[lseq_list[name.split(".")[0]].append(lseqs[i]) for i, name in enumerate(lnames) if len(lseqs[i]) > 211]
lseq_list = sequo.pick_random_family_member(lseqs_families, lseq_list)

lseqs_list = []
[lseqs_list.extend(lseq_list[i]) for i in lseq_list]
lseqs_list_flanks = []
# [lseqs_list_flanks.extend([i[2:69], i[-69:-2]]) for i in lseqs_list]

def randomise_seqs(seqs):
    new_seqs = []
    for seq in seqs:
        nts = list(seq)
        np.random.shuffle(nts)
        new_seqs.append("".join(nts))
    return new_seqs

def calc_sims(iterations, seqs, motifs):
    densities = collections.defaultdict(lambda: [])
    if len(iterations):
        np.random.seed()
        for i, iteration in enumerate(iterations):
            gen.print_parallel_status(i, iterations)
            rseqs = randomise_seqs(seqs)
            r = calc_densities(rseqs, motifs)
            for motif in r:
                densities[motif].append(r[motif])
    densities = {i: densities[i] for i in densities}

    return densities



def calc_densities(seqs, motifs):
    densities = {}
    for motif in motifs:
        densities[motif] = seqo.calc_motif_density(seqs, [motif])
    return densities

real = calc_densities(seqs_list, motifs)
lreal = calc_densities(lseqs_list, motifs)

sims = soc.run_simulation_function(list(range(5)), [seqs_list, motifs], calc_sims, sim_run = False)
lsims = soc.run_simulation_function(list(range(5)), [lseqs_list, motifs], calc_sims, sim_run = False)
nd = {motif: np.divide(real[motif] - np.mean(sims[motif]), np.mean(sims[motif])) for motif in real}
lnd = {motif: np.divide(lreal[motif] - np.mean(lsims[motif]), np.mean(lsims[motif])) for motif in lreal}


output_file = "temp_files/motifs.csv"
with open(output_file, "w") as outfile:
    outfile.write("motif,stop_motif,density,ldensity,nd,lnd\n")
    for motif in sorted(real):
        outfile.write("{0},{1},{2},{3},{4},{5}\n".format(motif, 1 if len(re.findall("(?=(TAA|TAG|TGA))", motif)) else 0, real[motif], lreal[motif], nd[motif], lnd[motif]))
