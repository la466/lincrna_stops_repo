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

nts = ["A", "C", "G", "T"]

motifs = "source_data/motif_sets/int3.txt"
# motifs = "source_data/motif_sets/RESCUE.txt"
# seqs_file = "clean_run/lincrna/lincRNA.multi_exon.exons.fasta"
# families_file = "clean_run/lincrna/lincRNA.multi_exon_families.bed"
# seqs_file = "clean_run/genome_sequences/lincrna/cabili/multi_exons.fasta"
families_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt"
seqs_file = "source_data/cabili_clean_alignments.fasta"

# seqs_file = "clean_run/genome_sequences/human/human.cds.clean_coding_exons.fasta"
# families_file = "clean_run/genome_sequences/human/human.cds.families.bed"
# seqs_file = "clean_run/genome_sequences/human.macaque.alignments.fasta"


motifs = sequo.read_motifs(motifs)

names, seqs = gen.read_fasta(seqs_file)
seq_list = collections.defaultdict(lambda: [])
[seq_list[name.split(".")[0]].append(seqs[i].split(",")) for i, name in enumerate(names)]

seq_list = sequo.pick_random_family_member(families_file, seq_list)


query = "".join([seq_list[i][0][0] for i in seq_list])
# print(query)

def get_random_overlaps(sequence, chunked_hits):
    new_overlaps = []
    for chunk in chunked_hits:
        chunk_length = len(chunk)
        start = np.random.choice(list(range(len(sequence) - len(chunk))))
        new_overlaps.extend(list(range(start, start+len(chunk))))
    new_overlaps = list(set(new_overlaps))
    hits = new_overlaps
    return hits

def get_nt_comp(sequence):
    nts = ["A", "C", "G", "T"]
    nt_comp = []
    for nt in nts:
        nt_comp.append(np.divide(sequence.count(nt), len(sequence)))
    return nt_comp

# def get_random_overlap(motif):
#     real_nt_comp = get_nt_comp(motif)
#     motif_length = len(motif)
#     got = False
#     while not got:
#         start = np.random.choice(list(range(len(query) - motif_length)))
#         motif = query[start:start+motif_length]
#         if real_nt_comp == get_nt_comp(motif):
#             got = True
#     print(motif)
#     return motif

def get_sub_rate(seq1, seq2):
    nts = ["A", "C", "G", "T"]
    subs = len([i for i in range(len(seq1)) if seq1[i] in nts and seq2[i] in nts and seq1[i] != seq2[i]])
    nts = len([i for i in range(len(seq1)) if seq1[i] in nts and seq2[i] in nts])
    return np.divide(subs, nts)

def calc_rate(iterations, randomise = None):
    outputs = []
    if iterations:
        np.random.seed()
        for iter, iteration in enumerate(iterations):
            gen.print_parallel_status(iter, iterations)

            all_hit_human_motifs = []
            all_hit_mac_motifs = []
            all_non_hit_human_motifs = []
            all_non_hit_mac_motifs = []
            all_hit_stop_human_motifs = []
            all_hit_stop_mac_motifs = []
            all_hit_non_stop_human_motifs = []
            all_hit_non_stop_mac_motifs = []
            all_non_hit_stop_human_motifs = []
            all_non_hit_stop_mac_motifs = []
            all_non_hit_non_stop_human_motifs = []
            all_non_hit_non_stop_mac_motifs = []


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

                    hit_stop_hits = [sequo.sequence_overlap_indicies(motif, stops) for motif in hit_human_motifs]
                    all_hit_stop_human_motifs.append("".join(["".join([hit_human_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(hit_stop_hits)]))
                    all_hit_stop_mac_motifs.append("".join(["".join([hit_mac_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(hit_stop_hits)]))

                    hit_non_stop_hits = [[i for i in list(range(len(motif))) if i not in hit_stop_hits[no]] for no, motif in enumerate(hit_human_motifs)]
                    all_hit_non_stop_human_motifs.append("".join(["".join([hit_human_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(hit_non_stop_hits)]))
                    all_hit_non_stop_mac_motifs.append("".join(["".join([hit_mac_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(hit_non_stop_hits)]))



                    non_hit_stop_hits = [sequo.sequence_overlap_indicies(motif, stops) for motif in non_hit_human_motifs]
                    all_non_hit_stop_human_motifs.append("".join(["".join([non_hit_human_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(non_hit_stop_hits)]))
                    all_non_hit_stop_mac_motifs.append("".join(["".join([non_hit_mac_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(non_hit_stop_hits)]))

                    non_hit_non_stop_hits = [[i for i in list(range(len(motif))) if i not in non_hit_stop_hits[no]] for no, motif in enumerate(non_hit_human_motifs)]
                    all_non_hit_non_stop_human_motifs.append("".join(["".join([non_hit_human_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(non_hit_non_stop_hits)]))
                    all_non_hit_non_stop_mac_motifs.append("".join(["".join([non_hit_mac_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(non_hit_non_stop_hits)]))

                    all_hit_human_motifs.append("".join(hit_human_motifs))
                    all_hit_mac_motifs.append("".join(hit_mac_motifs))
                    all_non_hit_human_motifs.append("".join(non_hit_human_motifs))
                    all_non_hit_mac_motifs.append("".join(non_hit_mac_motifs))


            all_hit_human_motifs = "".join(all_hit_human_motifs)
            all_hit_mac_motifs = "".join(all_hit_mac_motifs)
            all_non_hit_human_motifs = "".join(all_non_hit_human_motifs)
            all_non_hit_mac_motifs = "".join(all_non_hit_mac_motifs)
            all_hit_stop_human_motifs = "".join(all_hit_stop_human_motifs)
            all_hit_stop_mac_motifs = "".join(all_hit_stop_mac_motifs)
            all_hit_non_stop_human_motifs = "".join(all_hit_non_stop_human_motifs)
            all_hit_non_stop_mac_motifs = "".join(all_hit_non_stop_mac_motifs)
            all_non_hit_stop_human_motifs = "".join(all_non_hit_stop_human_motifs)
            all_non_hit_stop_mac_motifs = "".join(all_non_hit_stop_mac_motifs)
            all_non_hit_non_stop_human_motifs = "".join(all_non_hit_non_stop_human_motifs)
            all_non_hit_non_stop_mac_motifs = "".join(all_non_hit_non_stop_mac_motifs)



            hit_rate = get_sub_rate(all_hit_human_motifs, all_hit_mac_motifs)
            non_hit_rate = get_sub_rate(all_non_hit_human_motifs, all_non_hit_mac_motifs)
            hit_stop_rate = get_sub_rate(all_hit_stop_human_motifs, all_hit_stop_mac_motifs)
            hit_non_stop_rate = get_sub_rate(all_hit_non_stop_human_motifs, all_hit_non_stop_mac_motifs)
            non_hit_stop_rate = get_sub_rate(all_non_hit_stop_human_motifs, all_non_hit_stop_mac_motifs)
            non_hit_non_stop_rate = get_sub_rate(all_non_hit_non_stop_human_motifs, all_non_hit_non_stop_mac_motifs)
            output = [hit_rate, non_hit_rate, hit_stop_rate, hit_non_stop_rate, non_hit_stop_rate, non_hit_non_stop_rate]
            outputs.append(output)







            # outputs.append([ese_rate, non_ese_rate, relative_rate, "", ese_stop_rate, non_ese_stop_rate, expected_ese_stop_rate, diff1 ,"", ese_stop_rate, ese_non_stop_rate])

    return outputs

real = calc_rate([0])[0]
outputs = soc.run_simulation_function(list(range(1000)), [], calc_rate, kwargs_dict = {"randomise": True}, sim_run = False)

output_file = "temp_files/_rate1.csv"
with open(output_file, "w") as outfile:
    outfile.write("id,ese_rate,non_ese_rate,ese_stop_rate,ese_non_stop_rate,non_ese_stop_rate,non_ese_non_stop_rate\n")
    outfile.write("{0},{1}\n".format("real",",".join(gen.stringify(real))))
    [outfile.write("{0},{1}\n".format(i+1, ",".join(gen.stringify(output)))) for i, output in enumerate(outputs)]
# outputs = soc.run_simulation_function(list(range(20)), [], calc_rate, kwargs_dict = {"randomise": True}, sim_run = False)
# output_file = "temp_files/__temp_rates.csv"
#
# with open(output_file, "w") as outfile:
#     outfile.write("ese_rate,non_ese_rate,relative_rate,,ese_stop_rate,non_ese_stop_rate,expected_ese_stop_rate,diff,,ese_stop_rate,ese_non_stop_rate\n")
#     outfile.write("{0}\n".format(",".join(gen.stringify(real))))
#     [outfile.write("{0}\n".format(",".join(gen.stringify(output)))) for output in outputs]
