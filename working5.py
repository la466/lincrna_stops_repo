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
blast_file = "clean_run/genome_sequences/lincrna/GENCODE_CLS/single_exon_blast_output.csv"
single_exon_families_file = "clean_run/genome_sequences/lincrna/GENCODE_CLS/single_exons_families.txt"
multi_exon_transcripts = "clean_run/genome_sequences/lincrna/GENCODE_CLS/hsAllCompmerge.cage+polyASupported.full_transcripts.fasta"



# entries = gen.read_many_fields(exon_bed, "\t")
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

np.random.seed()
def randomise(seq):
    nts = list(seq)
    np.random.shuffle(nts)
    return "".join(nts)


def randomise_densities(ids, seq_list, iterations):

    random_densities = collections.defaultdict(lambda: [])
    if len(ids):

        for i, id in enumerate(ids):
            # print("ID {0}/{1}".format(i+1, len(ids)))

            for seq in seq_list[id]:
                sim_densities = []
                for iteration in range(iterations):
                    random_seq = randomise(seq)
                    sim_density = seqo.calc_motif_density([random_seq], stops)
                    sim_densities.append(sim_density)
                random_densities[id].append(sim_densities)

    random_densities = {i: random_densities[i] for i in random_densities}
    # print(random_densities)
    return random_densities

def calc_values(seq_list):

    densities = collections.defaultdict(lambda: [])
    for id in seq_list:
        for seq in seq_list[id]:
            density = seqo.calc_motif_density([seq], stops)
            densities[id].append(density)


    ids = list(seq_list)
    sim_outputs = gen.run_in_parallel(list(seq_list), ["foo", seq_list, 1000], randomise_densities)

    random_densities = collections.defaultdict(lambda: [])
    for output in sim_outputs:
        result = output.get()
        for id in result:
            random_densities[id].extend(result[id])

    random_densities = {i: random_densities[i] for i in random_densities}

    nds = collections.defaultdict(lambda: [])
    for id in densities:
        for i, exon_density in enumerate(densities[id]):
            nd = np.divide(exon_density - np.mean(random_densities[id][i]), np.mean(random_densities[id][i]))
            nds[id].append(nd)

    return densities, nds



names, seqs = gen.read_fasta(multi_exon_transcripts)
full_seq_list = {".".join(name.split(":")): seqs[i] for i, name in enumerate(names)}

single_exon_seq_list = get_seq_list(single_exon_fasta)
single_exon_seq_list = sequo.pick_random_family_member(single_exon_families_file, single_exon_seq_list)
single_densities, single_nds = calc_values(single_exon_seq_list)

multi_exon_seq_list = get_seq_list(multi_exon_fasta, full_seq_list = full_seq_list, with_chr = True)
multi_exon_seq_list = sequo.pick_random_family_member(multi_exon_families_file, multi_exon_seq_list)
multi_densities, multi_nds = calc_values(multi_exon_seq_list)



output_file = "temp_data/compare_densities1.csv"

def get_output(entry):
    output = []
    for id in entry:
        output.extend(entry[id])
    return gen.stringify(output)

with open(output_file, "w") as outfile:


    outfile.write("single_density,{0}\n".format(",".join(get_output(single_densities))))
    outfile.write("single_nd,{0}\n".format(",".join(get_output(single_nds))))
    outfile.write("multi_density,{0}\n".format(",".join(get_output(multi_densities))))
    outfile.write("multi_nd,{0}\n".format(",".join(get_output(multi_nds))))
