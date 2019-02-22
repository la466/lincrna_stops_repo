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
import _rs_conservation as rsc
import os
import zipfile
from useful_motif_sets import stops, codon_map
from regex_patterns import codon_pattern
import containers as cont
import itertools as it


ese_file = "source_data/motif_sets/int3.txt"

output_directory = "clean_run/tests/ese_densities"
gen.create_output_directories(output_directory)

lincrna_files = {
    "exons_file": "clean_run/genome_sequences/lincrna/cabili/exons.fasta",
    "introns_file": "clean_run/genome_sequences/lincrna/cabili/introns.fasta",
    "families_file": "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt" ,
    "output_file": "{0}/lincrna_ese_densities.csv".format(output_directory) ,
}

pc_files = {
    "exons_file": "clean_run/genome_sequences/human/human.cds.clean_coding_exons.fasta",
    "introns_file": "clean_run/genome_sequences/human/human.clean_introns.fasta",
    "families_file": "clean_run/genome_sequences/human/human.cds.families.txt" ,
    "output_file": "{0}/pc_ese_densities.csv".format(output_directory) ,
}

# filesets = [lincrna_files]
filesets = [lincrna_files, pc_files]

eses = sequo.read_motifs(ese_file)

for fileset in filesets:

    exon_names, exon_seqs = gen.read_fasta(fileset["exons_file"])
    exon_list = collections.defaultdict(lambda: [])
    [exon_list[name.split(".")[0]].append(exon_seqs[i]) for i, name in enumerate(exon_names) if len(exon_seqs[i])]
    exon_list = sequo.pick_random_family_member(fileset["families_file"], exon_list)

    intron_names, intron_seqs = gen.read_fasta(fileset["introns_file"])
    intron_list = collections.defaultdict(lambda: [])
    [intron_list[name.split(".")[0]].append(intron_seqs[i]) for i, name in enumerate(intron_names)]


    intron_sizes = {i: np.median([len(intron) for intron in intron_list[i]]) for i in intron_list}

    stop_eses = [i for i in eses if len(re.findall("(?=(TAA|TAG|TGA))", i))]
    non_stop_eses = [i for i in eses if i not in stop_eses]

    ese_densities = {i: seqo.calc_motif_density(exon_list[i], eses) for i in exon_list}
    stop_ese_densities = {i: seqo.calc_motif_density(exon_list[i], stop_eses) for i in exon_list}
    non_stop_ese_densities = {i: seqo.calc_motif_density(exon_list[i], non_stop_eses) for i in exon_list}

    with open(fileset["output_file"], "w") as outfile:
        outfile.write("id,intron_size,ese_density,stop_ese_density,non_stop_ese_density\n")
        [outfile.write("{0},{1},{2},{3},{4}\n".format(id, intron_sizes[id], ese_densities[id], np.divide(stop_ese_densities[id], len(stop_eses)), np.divide(non_stop_ese_densities[id], len(non_stop_eses)))) for id in ese_densities if id in intron_sizes]
