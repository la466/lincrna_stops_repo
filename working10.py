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
from useful_motif_sets import stops, codon_map, nucleotides
from regex_patterns import codon_pattern
import containers as cont
import itertools as it
from Bio import SeqIO
import scipy.stats

nts = ["A", "C", "G", "T"]

# motifs = "source_data/motif_sets/int3.txt"
motifs = "source_data/motif_sets/RESCUE.txt"
# seqs_file = "clean_run/lincrna/lincRNA.multi_exon.exons.fasta"
# families_file = "clean_run/lincrna/lincRNA.multi_exon_families.bed"
# seqs_file = "clean_run/genome_sequences/lincrna/cabili/multi_exons.fasta"
# families_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt"
# seqs_file = "source_data/cabili_clean_alignments.fasta"

# seqs_file = "clean_run/genome_sequences/human/human.cds.clean_coding_exons.fasta"
families_file = "clean_run/genome_sequences/human/human.cds.families.bed"
seqs_file = "clean_run/genome_sequences/human.macaque.alignments.fasta"


motifs = sequo.read_motifs(motifs)

names, seqs = gen.read_fasta(seqs_file)
seq_list = collections.defaultdict(lambda: [])
[seq_list[name.split(".")[0]].append(seqs[i].split(",")) for i, name in enumerate(names)]

seq_list = sequo.pick_random_family_member(families_file, seq_list)


query = "".join([seq_list[i][0][0] for i in seq_list])
# print(query)


def get_nt_comp(sequence):
    nts = ["A", "C", "G", "T"]
    nt_comp = []
    for nt in nts:
        nt_comp.append(np.divide(sequence.count(nt), len(sequence)))
    return nt_comp

# def calc_rate(iterations, randomise = None):
#     outputs = []
#     if iterations:
#         np.random.seed()
#         for iter, iteration in enumerate(iterations):
#             gen.print_parallel_status(iter, iterations)
#
#             all_hit_human_motifs = []
#             all_hit_mac_motifs = []
#             all_non_hit_human_motifs = []
#             all_non_hit_mac_motifs = []
#             all_hit_stop_human_motifs = []
#             all_hit_stop_mac_motifs = []
#             all_hit_non_stop_human_motifs = []
#             all_hit_non_stop_mac_motifs = []
#             all_non_hit_stop_human_motifs = []
#             all_non_hit_stop_mac_motifs = []
#             all_non_hit_non_stop_human_motifs = []
#             all_non_hit_non_stop_mac_motifs = []
#
#             for count, id in enumerate(seq_list):
#                 if count:
#                     sequences = seq_list[id][0]
#                     human = sequences[0]
#                     mac = sequences[1]
#
#                     hits = sequo.sequence_overlap_indicies(human, motifs)
#                     chunked_hits = sequo.chunk_overlaps(hits)
#
#                     if randomise:
#                         hits = get_random_overlaps(human, chunked_hits)
#                         chunked_hits = sequo.chunk_overlaps(hits)
#
#                     non_hits = [i for i in range(len(human)) if i not in hits]
#                     chunked_non_hits = sequo.chunk_overlaps(non_hits)
#
#                     hit_human_motifs = ["".join([human[i] for i in chunk]) for chunk in chunked_hits]
#                     hit_mac_motifs = ["".join([mac[i] for i in chunk]) for chunk in chunked_hits]
#                     non_hit_human_motifs = ["".join([human[i] for i in chunk]) for chunk in chunked_non_hits]
#                     non_hit_mac_motifs = ["".join([mac[i] for i in chunk]) for chunk in chunked_non_hits]
#
#                     hit_stop_hits = [sequo.sequence_overlap_indicies(motif, stops) for motif in hit_human_motifs]
#                     all_hit_stop_human_motifs.append("".join(["".join([hit_human_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(hit_stop_hits)]))
#                     all_hit_stop_mac_motifs.append("".join(["".join([hit_mac_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(hit_stop_hits)]))
#
#                     hit_non_stop_hits = [[i for i in list(range(len(motif))) if i not in hit_stop_hits[no]] for no, motif in enumerate(hit_human_motifs)]
#                     all_hit_non_stop_human_motifs.append("".join(["".join([hit_human_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(hit_non_stop_hits)]))
#                     all_hit_non_stop_mac_motifs.append("".join(["".join([hit_mac_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(hit_non_stop_hits)]))
#
#
#
#                     non_hit_stop_hits = [sequo.sequence_overlap_indicies(motif, stops) for motif in non_hit_human_motifs]
#                     all_non_hit_stop_human_motifs.append("".join(["".join([non_hit_human_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(non_hit_stop_hits)]))
#                     all_non_hit_stop_mac_motifs.append("".join(["".join([non_hit_mac_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(non_hit_stop_hits)]))
#
#                     non_hit_non_stop_hits = [[i for i in list(range(len(motif))) if i not in non_hit_stop_hits[no]] for no, motif in enumerate(non_hit_human_motifs)]
#                     all_non_hit_non_stop_human_motifs.append("".join(["".join([non_hit_human_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(non_hit_non_stop_hits)]))
#                     all_non_hit_non_stop_mac_motifs.append("".join(["".join([non_hit_mac_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(non_hit_non_stop_hits)]))
#
#                     all_hit_human_motifs.append("".join(hit_human_motifs))
#                     all_hit_mac_motifs.append("".join(hit_mac_motifs))
#                     all_non_hit_human_motifs.append("".join(non_hit_human_motifs))
#                     all_non_hit_mac_motifs.append("".join(non_hit_mac_motifs))
#
#
#             all_hit_human_motifs = "".join(all_hit_human_motifs)
#             all_hit_mac_motifs = "".join(all_hit_mac_motifs)
#             all_non_hit_human_motifs = "".join(all_non_hit_human_motifs)
#             all_non_hit_mac_motifs = "".join(all_non_hit_mac_motifs)
#             all_hit_stop_human_motifs = "".join(all_hit_stop_human_motifs)
#             all_hit_stop_mac_motifs = "".join(all_hit_stop_mac_motifs)
#             all_hit_non_stop_human_motifs = "".join(all_hit_non_stop_human_motifs)
#             all_hit_non_stop_mac_motifs = "".join(all_hit_non_stop_mac_motifs)
#             all_non_hit_stop_human_motifs = "".join(all_non_hit_stop_human_motifs)
#             all_non_hit_stop_mac_motifs = "".join(all_non_hit_stop_mac_motifs)
#             all_non_hit_non_stop_human_motifs = "".join(all_non_hit_non_stop_human_motifs)
#             all_non_hit_non_stop_mac_motifs = "".join(all_non_hit_non_stop_mac_motifs)
#
#
#
#             hit_rate = get_sub_rate(all_hit_human_motifs, all_hit_mac_motifs)
#             non_hit_rate = get_sub_rate(all_non_hit_human_motifs, all_non_hit_mac_motifs)
#             hit_stop_rate = get_sub_rate(all_hit_stop_human_motifs, all_hit_stop_mac_motifs)
#             hit_non_stop_rate = get_sub_rate(all_hit_non_stop_human_motifs, all_hit_non_stop_mac_motifs)
#             non_hit_stop_rate = get_sub_rate(all_non_hit_stop_human_motifs, all_non_hit_stop_mac_motifs)
#             non_hit_non_stop_rate = get_sub_rate(all_non_hit_non_stop_human_motifs, all_non_hit_non_stop_mac_motifs)
#             output = [hit_rate, non_hit_rate, hit_stop_rate, hit_non_stop_rate, non_hit_stop_rate, non_hit_non_stop_rate]
#             outputs.append(output)
#
#     return outputs

def calc_dint_rate(iterations, randomise = None):
    outputs = []
    if iterations:
        np.random.seed()
        for iter, iteration in enumerate(iterations):
            gen.print_parallel_status(iter, iterations)

            all_hit_human_motifs = []
            all_hit_mac_motifs = []
            all_non_hit_human_motifs = []
            all_non_hit_mac_motifs = []

            hit_dint_count = collections.defaultdict(lambda: 0)
            hit_dint_sub = collections.defaultdict(lambda: 0)
            non_hit_dint_count = collections.defaultdict(lambda: 0)
            non_hit_dint_sub = collections.defaultdict(lambda: 0)

            for count, id in enumerate(seq_list):
                if count:
                    sequences = seq_list[id][0]
                    human = sequences[0]
                    mac = sequences[1]

                    hits = sequo.sequence_overlap_indicies(human, motifs)
                    chunked_hits = sequo.chunk_overlaps(hits)

                    if randomise:
                        hits = get_random_overlaps(human, chunked_hits)
                        chunked_hits = sequo.chunk_overlaps(hits)

                    non_hits = [i for i in range(len(human)) if i not in hits]
                    chunked_non_hits = sequo.chunk_overlaps(non_hits)

                    hit_human_motifs = ["".join([human[i] for i in chunk]) for chunk in chunked_hits]
                    hit_mac_motifs = ["".join([mac[i] for i in chunk]) for chunk in chunked_hits]
                    non_hit_human_motifs = ["".join([human[i] for i in chunk]) for chunk in chunked_non_hits]
                    non_hit_mac_motifs = ["".join([mac[i] for i in chunk]) for chunk in chunked_non_hits]

                    for i, hit in enumerate(hit_human_motifs):
                        human_motif = hit
                        mac_motif = hit_mac_motifs[i]

                        for pos in range(len(human_motif)-1):
                            if human_motif[pos] in nucleotides and mac_motif[pos] in nucleotides:
                                hit_dint_count[human_motif[pos:pos+2]] += 1
                                if human_motif[pos:pos+2] != mac_motif[pos:pos+2]:
                                    hit_dint_sub[human_motif[pos:pos+2]] += 1

                    for i, non_hit in enumerate(non_hit_human_motifs):
                        human_motif = non_hit
                        mac_motif = non_hit_mac_motifs[i]

                        for pos in range(len(human_motif)-1):
                            if human_motif[pos] in nucleotides and mac_motif[pos] in nucleotides:
                                non_hit_dint_count[human_motif[pos:pos+2]] += 1
                                if human_motif[pos:pos+2] != mac_motif[pos:pos+2]:
                                    non_hit_dint_sub[human_motif[pos:pos+2]] += 1

            dint_sub_rates = {}
            non_dint_sub_rates = {}
            for dint in sorted(hit_dint_count):
                dint_sub_rates[dint] = np.divide(hit_dint_sub[dint], hit_dint_count[dint])
                non_dint_sub_rates[dint] = np.divide(non_hit_dint_sub[dint], non_hit_dint_count[dint])

            stop_dints = ["TA", "TG", "GA", "AA"]
            non_stop_dints = [i for i in dint_sub_rates if i not in stop_dints]

            total_hit_dints = sum(hit_dint_count.values())
            total_hit_dint_subs = sum(hit_dint_sub.values())

            stop_dint_subs = sum([hit_dint_sub[i] for i in hit_dint_sub if i in stop_dints])
            non_stop_dint_subs = sum([hit_dint_sub[i] for i in hit_dint_sub if i in non_stop_dints])
            total_stop_dints = sum([hit_dint_count[i] for i in hit_dint_count if i in stop_dints])
            total_non_stop_dints = sum([hit_dint_count[i] for i in hit_dint_count if i in non_stop_dints])
            expected_stop_dints = np.divide(total_stop_dints, total_hit_dints)*total_hit_dint_subs
            expected_non_stop_dints = np.divide(total_non_stop_dints, total_hit_dints)*total_hit_dint_subs

            print(expected_stop_dints, stop_dint_subs, expected_non_stop_dints, non_stop_dint_subs)
            print(scipy.stats.chisquare([stop_dint_subs, non_stop_dint_subs], f_exp = [expected_stop_dints, expected_non_stop_dints]))


            total_hit_dints = sum(non_hit_dint_count.values())
            total_hit_dint_subs = sum(non_hit_dint_sub.values())
            stop_dint_subs = sum([non_hit_dint_sub[i] for i in non_hit_dint_sub if i in stop_dints])
            non_stop_dint_subs = sum([non_hit_dint_sub[i] for i in non_hit_dint_sub if i in non_stop_dints])
            total_stop_dints = sum([non_hit_dint_count[i] for i in non_hit_dint_count if i in stop_dints])
            total_non_stop_dints = sum([non_hit_dint_count[i] for i in non_hit_dint_count if i in non_stop_dints])
            expected_stop_dints = np.divide(total_stop_dints, total_hit_dints)*total_hit_dint_subs
            expected_non_stop_dints = np.divide(total_non_stop_dints, total_hit_dints)*total_hit_dint_subs

            print(expected_stop_dints, stop_dint_subs, expected_non_stop_dints, non_stop_dint_subs)
            print(scipy.stats.chisquare([stop_dint_subs, non_stop_dint_subs], f_exp = [expected_stop_dints, expected_non_stop_dints]))

            hit_stop_dint_rate = np.mean([dint_sub_rates[i] for i in stop_dints])
            hit_non_stop_dint_rate = np.mean([dint_sub_rates[i] for i in non_stop_dints])
            non_hit_stop_dint_rate = np.mean([non_dint_sub_rates[i] for i in stop_dints])
            non_hit_non_stop_dint_rate = np.mean([non_dint_sub_rates[i] for i in non_stop_dints])

            print(hit_stop_dint_rate)
            print(hit_non_stop_dint_rate)
            print(np.divide(hit_stop_dint_rate - hit_non_stop_dint_rate, hit_non_stop_dint_rate)*100)
            print(non_hit_stop_dint_rate)
            print(non_hit_non_stop_dint_rate)
            print(np.divide(non_hit_stop_dint_rate - non_hit_non_stop_dint_rate, non_hit_non_stop_dint_rate)*100)

            [print(i, dint_sub_rates[i]) for i in sorted(dint_sub_rates)]
            [print(i, dint_sub_rates[i]) for i in sorted(non_dint_sub_rates)]


    return outputs







# real = calc_rate([0])[0]
# outputs = soc.run_simulation_function(list(range(1000)), [], calc_rate, kwargs_dict = {"randomise": True}, sim_run = False)
#
# output_file = "temp_files/_rate1.csv"
# with open(output_file, "w") as outfile:
#     outfile.write("id,ese_rate,non_ese_rate,ese_stop_rate,ese_non_stop_rate,non_ese_stop_rate,non_ese_non_stop_rate\n")
#     outfile.write("{0},{1}\n".format("real",",".join(gen.stringify(real))))
#     [outfile.write("{0},{1}\n".format(i+1, ",".join(gen.stringify(output)))) for i, output in enumerate(outputs)]


real = calc_dint_rate([0])
# soc.run_simulation_function(list(range(5)), [], calc_dint_rate, kwargs_dict = {"randomise": True}, sim_run = False)
