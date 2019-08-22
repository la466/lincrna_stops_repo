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

# single_exon_sequences = "clean_run/genome_sequences/lincrna/cabili/single_exons.fasta"
# ese_file = "source_data/motif_sets/int3.txt"
# motifs = sequo.read_motifs(ese_file)
#
# stop_motifs = [i for i in motifs if len(re.findall("(TAA|TAG|TGA)", i)) > 0]
#
# names, seqs = gen.read_fasta(single_exon_sequences)
# print(names)

exon_bed = "source_data/hsAllCompmerge.cage+polyASupported.bed"
single_exon_bed = "clean_run/genome_sequences/lincrna/GENCODE_CLS/single_exons.bed"
genome_fasta = "../source_data/genomes/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
multi_exon_fasta = "clean_run/genome_sequences/lincrna/GENCODE_CLS/hsAllCompmerge.cage+polyASupported.bed.filtered.exons.fasta"
multi_exon_families_file = "clean_run/genome_sequences/lincrna/GENCODE_CLS/families.txt"
single_exon_fasta = "clean_run/genome_sequences/lincrna/GENCODE_CLS/single_exons.fasta"
blast_file = "clean_run/genome_sequences/lincrna/GENCODE_CLS/single_exon_blast_output.csv"
single_exon_families_file = "clean_run/genome_sequences/lincrna/GENCODE_CLS/single_exons_families.txt"


# entries = gen.read_many_fields(exon_bed, "\t")
# for entry in entries[:4]:
#     print(entry)
#
# exons = [int(i[-3]) for i in entries]
# single_exons = [i for i in exons if i == 1]
#
# with open(single_exon_bed, "w") as outfile:
#     for entry in entries:
#         if int(entry[-3]) == 1:
#             outfile.write("{0}\t{1}\n".format(entry[0].strip("chr"), "\t".join(entry[1:])))
#
# fo.fasta_from_intervals(single_exon_bed, single_exon_fasta, genome_fasta, names = True)
# cons.filter_families(single_exon_fasta, blast_file, single_exon_families_file, database_path = None, clean_run = None)

np.random.seed()
def randomise(seq):
    nts = list(seq)
    np.random.shuffle(nts)
    return "".join(nts)


def get_seq_list(fasta_file, with_chr = False):
    names, seqs = gen.read_fasta(fasta_file)
    print(len(names))
    seq_list = collections.defaultdict(lambda: [])
    for i, name in enumerate(names):
        id = name.split("(")[0]
        if with_chr:
            id = ".".join(id.split(".")[:-1])
        seq_list[id].append(seqs[i])
    return seq_list


def calc_values(seq_list):

    densities = collections.defaultdict(lambda: [])
    for id in seq_list:
        for seq in seq_list[id]:
            density = seqo.calc_motif_density([seq], stops)
            densities[id].append(density)

    # random_densities = collections.defaultdict(lambda: [])
    # for id in seq_list:
    #     for seq in seq_list[id]:
    #         sim_densities = []
    #         for i in range(1):
    #             random_seq = randomise(seq)
    #             sim_density = seqo.calc_motif_density([random_seq], stops)
    #             sim_densities.append(sim_density)
    #         random_densities[id].append(sim_densities)
    #
    # nds = collections.defaultdict(lambda: [])
    # for id in densities:
    #     for i, exon_density in enumerate(densities[id]):
    #         nd = np.divide(exon_density - np.mean(random_densities[id][i]), np.mean(random_densities[id][i]))
    #         nds[id].append(nd)

    return densities
    # return densities, nds


# single_exon_seq_list = get_seq_list(single_exon_fasta)
# single_densities, single_nds = calc_values(single_exon_seq_list)

multi_exon_seq_list = get_seq_list(multi_exon_fasta, with_chr = True)
multi_densities = calc_values(multi_exon_seq_list)

print(len(multi_densities))
[print(i) for i in sorted(multi_densities)]


# single_exon_families = gen.read_many_fields(single_exon_families_file, "\t")
# single_densities = sequo.group_family_results(single_densities, single_exon_families)
# single_nds = sequo.group_family_results(single_nds, single_exon_families)

multi_exon_families = gen.read_many_fields(multi_exon_families_file, "\t")
multi_exon_families = [[".".join(i.split(":")) for i in family] for family in multi_exon_families]
print(len(multi_exon_families))
# # print(multi_exon_families)
multi_densities = sequo.group_family_results(multi_densities, multi_exon_families)
# multi_nds = sequo.group_family_results(multi_nds, multi_exon_families)
#
# print(len(single_densities), len(single_nds))
print(len(multi_densities))


# output_file = "temp_data/compare_densities.csv"
#
# def get_output(entry):
#     output = []
#     for i in entry:
#         output.append(np.median(entry[i]))
#     return output
#
# with open(output_file, "w") as outfile:
#
#     outfile.write("single_density,{0}\n".format(",".join(gen.stringify(get_output(single_densities)))))
#     outfile.write("single_nd,{0}\n".format(",".join(gen.stringify(get_output(single_nds)))))
#     outfile.write("multi_density,{0}\n".format(",".join(gen.stringify(get_output(multi_densities)))))
#     outfile.write("multi_nd,{0}\n".format(",".join(gen.stringify(get_output(multi_nds)))))
