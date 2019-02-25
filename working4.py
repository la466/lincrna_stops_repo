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


ese_file = "source_data/motif_sets/RESCUE.txt"

output_directory = "clean_run/tests/ese_densities"
gen.create_output_directories(output_directory)

lincrna_files = {
    "fileset_prefix": "lincrna",
    "exons_file": "clean_run/genome_sequences/lincrna/cabili/exons.fasta",
    "introns_file": "clean_run/genome_sequences/lincrna/cabili/introns.fasta",
    "families_file": "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt" ,
    "output_file": "{0}/{1}_lincrna_ese_densities.csv".format(output_directory, ese_file.split("/")[-1].split(".")[0]) ,
}

pc_files = {
    "fileset_prefix": "pc",
    "exons_file": "clean_run/genome_sequences/human/human.cds.clean_coding_exons.fasta",
    "introns_file": "clean_run/genome_sequences/human/human.clean_introns.fasta",
    "families_file": "clean_run/genome_sequences/human/human.cds.families.bed" ,
}

# filesets = [lincrna_files]
# filesets = [pc_files]
filesets = [
    lincrna_files,
    pc_files
]

ese_sets = [
    "source_data/motif_sets/int3.txt",
    "source_data/motif_sets/RESCUE.txt",
    "source_data/motif_sets/ke400.txt",
    "source_data/motif_sets/ESR.txt",
    "source_data/motif_sets/PESE.txt",
    "source_data/motif_sets/combined_eses.txt",
]


def run_densities(ese_sets, exon_list, intron_list, exon_sizes, intron_sizes, fileset, all_sequences = None, flanks = None):
    if ese_sets:
        for i, ese_set in enumerate(ese_sets):
            gen.print_parallel_status(i, ese_sets)
            eses = sequo.read_motifs(ese_set)

            stop_eses = [i for i in eses if len(re.findall("(?=(TAA|TAG|TGA))", i))]
            non_stop_eses = [i for i in eses if i not in stop_eses]

            ese_densities = {i: seqo.calc_motif_density(exon_list[i], eses) for i in exon_list}
            stop_ese_densities = {i: seqo.calc_motif_density(exon_list[i], stop_eses) for i in exon_list}
            non_stop_ese_densities = {i: seqo.calc_motif_density(exon_list[i], non_stop_eses) for i in exon_list}

            if all_sequences:
                if flanks:
                    output_file = "{0}/{1}_{2}_ese_densities_all_sequences_flanks.csv".format(output_directory, ese_set.split("/")[-1].split(".")[0], fileset["fileset_prefix"])
                else:
                    output_file = "{0}/{1}_{2}_ese_densities_all_sequences.csv".format(output_directory, ese_set.split("/")[-1].split(".")[0], fileset["fileset_prefix"])
            else:
                if flanks:
                    output_file = "{0}/{1}_{2}_ese_densities_flanks.csv".format(output_directory, ese_set.split("/")[-1].split(".")[0], fileset["fileset_prefix"])
                else:
                    output_file = "{0}/{1}_{2}_ese_densities.csv".format(output_directory, ese_set.split("/")[-1].split(".")[0], fileset["fileset_prefix"])
            with open(output_file, "w") as outfile:
                outfile.write("id,exon_size,intron_size,ese_density,stop_ese_density,non_stop_ese_density\n")
                [outfile.write("{0},{1},{2},{3},{4},{5}\n".format(id, exon_sizes[id], intron_sizes[id], ese_densities[id], stop_ese_densities[id], non_stop_ese_densities[id])) for id in ese_densities if id in intron_sizes]

    # return for function
    return []

def calculate_densities(filesets, all_sequences = False, flanks = None):

    for fileset in filesets:
        exon_names, exon_seqs = gen.read_fasta(fileset["exons_file"])


        if flanks:
            exon_list = collections.defaultdict(lambda: collections.defaultdict())
            for i, name in enumerate(exon_names):
                if len(exon_seqs[i]) > 207:
                    exon_list[name.split(".")[0]][int(name.split(".")[-1].split("(")[0])] = exon_seqs[i]
            exon_list = {i: exon_list[i] for i in exon_list}
            exon_sizes = {i: np.median([len(exon_list[i][exon_id]) for exon_id in exon_list[i]]) for i in exon_list}

            intron_names, intron_seqs = gen.read_fasta(fileset["introns_file"])
            intron_list = collections.defaultdict(lambda: [])
            # get only those that flank an included exon
            for i, name in enumerate(intron_names):
                intron_id = name.split(".")[0]
                if intron_id in exon_list:
                    if int(name.split(".")[-1].split("(")[0].split("-")[0]) in exon_list[intron_id] or int(name.split(".")[-1].split("(")[0].split("-")[1]) in exon_list[intron_id]:
                        intron_list[name.split(".")[0]].append(intron_seqs[i])
            intron_list = {i: intron_list[i] for i in intron_list}
            intron_sizes = {i: np.median([len(intron) for intron in intron_list[i]]) for i in intron_list}


            temp_exon_list = collections.defaultdict(lambda: [])
            for id in exon_list:
                flanks = []
                for exon_id in exon_list[id]:
                    flanks.append(exon_list[id][exon_id][2:69])
                    flanks.append(exon_list[id][exon_id][-69:-2])
                temp_exon_list[id] = flanks
            exon_list = temp_exon_list

        else:

            [exon_list[name.split(".")[0]].append(exon_seqs[i]) for i, name in enumerate(exon_names)]
            exon_list = {i: exon_list[i] for i in exon_list}
            exon_sizes = {i: np.median([len(exon) for exon in exon_list[i]]) for i in exon_list}

            intron_names, intron_seqs = gen.read_fasta(fileset["introns_file"])
            intron_list = collections.defaultdict(lambda: [])
            [intron_list[name.split(".")[0]].append(intron_seqs[i]) for i, name in enumerate(intron_names) if name.split(".")[0] in exon_list]
            intron_list = {i: intron_list[i] for i in intron_list}
            intron_sizes = {i: np.median([len(intron) for intron in intron_list[i]]) for i in intron_list}

        if not all_sequences:
            exon_list = sequo.pick_random_family_member(fileset["families_file"], exon_list)

        soc.run_simulation_function(ese_sets, [exon_list, intron_list, exon_sizes, intron_sizes, fileset], run_densities, kwargs_dict = {"all_sequences": all_sequences, "flanks": flanks}, sim_run = False)

# calculate_densities(filesets)
# calculate_densities(filesets, all_sequences = True)
calculate_densities(filesets, flanks = True)