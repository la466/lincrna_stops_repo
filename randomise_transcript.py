import generic as gen
import seq_ops as seqo
import sequence_ops as sequo
import sim_ops as simo
import sim_ops_containers as soc
import containers as cont
import ops
import numpy as np
import collections
import pandas as pd
import conservation as cons
import file_ops as fo
import os
import copy
from useful_motif_sets import codon_map, stops
import re
import itertools as it
from scipy import stats
import random

exons_file = "clean_run/genome_sequences/human/human.cds.clean_coding_exons.fasta"
introns_file = "clean_run/genome_sequences/human/human.clean_introns.fasta"
families_file = "clean_run/genome_sequences/human/human.cds.families.bed"

output_file = "clean_run/tests/purine_content/exon_intron_purine_content_with_ese.csv"

exon_names, exon_seqs = gen.read_fasta(exons_file)
intron_names, intron_seqs = gen.read_fasta(introns_file)

introns = collections.defaultdict(lambda: collections.defaultdict())
for i, name in enumerate(intron_names):
    introns[name.split(".")[0]][int(name.split(".")[1].split("(")[0].split("-")[0])] = intron_seqs[i]

exons = collections.defaultdict(lambda: collections.defaultdict())
for i, name in enumerate(exon_names):
    if name.split(".")[0] in introns:
        exons[name.split(".")[0]][int(name.split(".")[1].split("(")[0])] = exon_seqs[i]


exons = {i: exons[i] for i in exons}
introns = {i: introns[i] for i in introns}


def wilcox_test(filepath, col1, col2):
    data = pd.read_csv(filepath)
    return stats.wilcoxon(data[col1], data[col2])

def get_median(filepath, col):
    data = pd.read_csv(filepath)
    return data[col].median()

def get_sequences(simulations, exons, introns, temp_directory):

    outputs = []

    if len(simulations):
        for i, simulation in enumerate(simulations):
            gen.print_parallel_status(i, simulations)

            temp_file1 = "{0}/run_content_{1}.txt".format(temp_directory, simulation)
            temp_file2 = "{0}/run_{1}.txt".format(temp_directory, simulation)


            if temp_file1 not in os.listdir(temp_directory):
                with open(temp_file1, "w") as temp_output:
                    temp_output.write("id,exon,intron\n")
                    if simulation != "real":
                        np.random.seed()
                        for i, id in enumerate(exons):
                            if id in introns:
                                sequence = []
                                for exon_id in sorted(exons[id]):
                                    sequence.append(exons[id][exon_id])
                                    if exon_id in introns[id]:
                                        sequence.append(introns[id][exon_id])
                                sequence = "".join(sequence)
                                exon_coordinates = [[sequence.index(exons[id][i]), sequence.index(exons[id][i]) + len(exons[id][i])] for i in sorted(exons[id])]
                                intron_coordinates = [[sequence.index(introns[id][i]), sequence.index(introns[id][i]) + len(introns[id][i])] for i in sorted(introns[id])]
                                seq_nts = list(sequence)
                                np.random.shuffle(seq_nts)
                                random_seq = "".join(seq_nts)
                                new_exons = [random_seq[i[0]:i[1]] for i in exon_coordinates]
                                new_introns = [random_seq[i[0]:i[1]] for i in intron_coordinates]
                                exon_content = sequo.calc_purine_content(new_exons)
                                intron_content = sequo.calc_purine_content(new_introns)
                                temp_output.write("{0},{1},{2}\n".format(id, exon_content, intron_content))

                    else:
                        for i, id in enumerate(exons):
                            if id in introns:
                                exon_list = [exons[id][i] for i in exons[id]]
                                intron_list = [introns[id][i] for i in introns[id]]
                                exon_content = sequo.calc_purine_content(exon_list)
                                intron_content = sequo.calc_purine_content(intron_list)
                                temp_output.write("{0},{1},{2}\n".format(id, exon_content, intron_content))

            test = wilcox_test(temp_file1, "exon", "intron")

            outputs.append(temp_file2)
            with open(temp_file2, "w") as outfile:
                outfile.write("{0},{1},{2},{3},{4}\n".format(simulation, get_median(temp_file1, "exon"), get_median(temp_file1, "intron"), test.statistic, test.pvalue))
            gen.remove_file(temp_file1)

    return outputs


simulations = ['real']
simulations.extend(list(range(1000)))

temp_directory = "temp_purine_shuffle_dir"
gen.create_output_directories(temp_directory)

args = [exons, introns, temp_directory]
outputs = soc.run_simulation_function(simulations, args, get_sequences, sim_run = False)

with open(output_file, "w") as outfile:
    outfile.write("id,median_exon,median_intron,test_statistic,pvalue\n")

    for file in os.listdir(temp_directory):
        filepath = "{0}/{1}".format(temp_directory, file)
        data = gen.read_many_fields(filepath, "\t")[0]
        outfile.write("{0}\n".format(",".join(data)))

gen.remove_directory(temp_directory)
