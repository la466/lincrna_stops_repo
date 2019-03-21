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
from scipy.stats import chisquare

# eses = sequo.read_motifs("source_data/motif_sets/int3.txt")

ese_sets = {
    # "int3": sequo.read_motifs("source_data/motif_sets/int3.txt"),
    "RESCUE": sequo.read_motifs("source_data/motif_sets/RESCUE.txt"),
}

# eses = sequo.read_motifs("source_data/motif_sets/ke400.txt")

pc_fileset = {
    "fileset_prefix": "protein_coding",
    "exons_file": "clean_run/genome_sequences/human/human.cds.clean_coding_exons.fasta",
    # "introns_file": "clean_run/genome_sequences/human/human.clean_introns.fasta",
    "families_file": "clean_run/genome_sequences/human/human.cds.families.bed" ,
}

lincrna_files = {
    "fileset_prefix": "lincrna",
    "exons_file": "clean_run/genome_sequences/lincrna/cabili/exons.fasta",
    # "introns_file": "clean_run/genome_sequences/lincrna/cabili/introns.fasta",
    "families_file": "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt" ,
}

filesets = [
    # pc_fileset,
    lincrna_files,
]

def get_sequences(fileset):
    exon_names, exon_seqs = gen.read_fasta(fileset["exons_file"])
    exon_list = collections.defaultdict(lambda: [])
    for i, name in enumerate(exon_names):
        exon_list[name.split(".")[0]].append(exon_seqs[i])
    exon_list = {i: exon_list[i] for i in exon_list}
    return exon_list

def calc_overlaps(runs, motifs, exon_list, fileset, sim = None):
    outputs = []
    if len(runs):

        # stop_eses = [i for i in motifs if len(re.findall("(?=(TAA|TAG|TGA))", i))]
        # non_stop_eses = [i for i in motifs if i not in stop_eses]

        for i, run in enumerate(runs):
            np.random.seed()
            gen.print_parallel_status(i, runs)
            data_exon_list = sequo.pick_random_family_member(fileset["families_file"], exon_list)

            temp_exon_list = []
            for id in data_exon_list:
                for seq in data_exon_list[id]:
                    temp_exon_list.append(seq)
            data_exon_list = temp_exon_list

            motif_hits = gen.flatten([sequo.return_overlap_motifs(seq, motifs) for seq in data_exon_list])

            if sim:
                overlap_hit_length = len([i for i in motif_hits if len(i) > 6])
                overlap_hits = np.random.choice(motif_hits, overlap_hit_length, replace = False)
                non_overlap_hits = motif_hits
            else:
                overlap_hits = [i for i in motif_hits if len(i) > 6]
                non_overlap_hits = [i for i in motif_hits if len(i) == 6]

            stop_non_overlaps = seqo.calc_motif_density(overlap_hits, stops)
            stop_overlaps = seqo.calc_motif_density(non_overlap_hits, stops)

            outputs.append([stop_non_overlaps, stop_overlaps])




    return outputs


def run_test(runs, ese_sets, filesets):
    for set in ese_sets:
        for fileset in filesets:
            np.random.seed()

            exon_list = get_sequences(fileset)

            motifs = ese_sets[set]

            real_output = soc.run_simulation_function(list(range(runs)), [motifs, exon_list, fileset], calc_overlaps, sim_run = False)[0]
            outputs = soc.run_simulation_function(list(range(runs)), [motifs, exon_list, fileset], calc_overlaps, kwargs_dict = {"sim": True}, sim_run = False)

            output_directory = "clean_run/tests/motif_overlaps/"
            gen.create_output_directories(output_directory)

            output_file = "{0}/{1}_{2}_motif_overlap_density.csv".format(output_directory, fileset["fileset_prefix"], set)
            with open(output_file, "w") as outfile:
                outfile.write("sim_id,non_overlap_stop_density,overlap_stop_density\n")
                outfile.write("real,{0},{1}\n".format(real_output[0], real_output[1]))
                for i, output in enumerate(outputs):
                    outfile.write("sim_{0},{1},{2}\n".format(i+1,output[0], output[1]))
            # pvals = [i[-1].pvalue for i in outputs]
            # adj_pvals = [i*runs for i in pvals]
            #
            # output_directory = "clean_run/tests/motif_overlaps/"
            # gen.create_output_directories(output_directory)
            #
            # output_file = "{0}/{1}_{2}_motif_overlap_chisquare.csv".format(output_directory, fileset["fileset_prefix"], set)
            # with open(output_file, "w") as outfile:
            #     outfile.write("id,expected_stop_overlaps,observed_stop_overlaps,expected_non_stop_overlaps,observed_non_stop_overlaps,mean_stop_overlap_length,mean_non_stop_overlap_length,,chi_statistic,p_value,adj_p\n")
            #     for i, output in enumerate(outputs):
            #         outfile.write("run_{0},{1},,{2},{3},{4}\n".format(i+1, ",".join(gen.stringify(output[:-1])), output[-1].statistic, output[-1].pvalue, adj_pvals[i]))


run_test(1000, ese_sets, filesets)
