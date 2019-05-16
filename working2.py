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

### masking eses for orf lengths



# file = "clean_run/genome_sequences/human/human.cds.clean_coding_exons.fasta"
# cds_file = "clean_run/genome_sequences/human/human.cds.multi_exons.fasta"
# families_file = "clean_run/genome_sequences/human/human.cds.families.bed"
file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_transcript_sequences.fasta"
families_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt"

test_file = "clean_run/tests/lincrna/orf_length_sim/cabili_sim_orf_lengths_zs_grouped.csv"

entries = gen.read_many_fields(test_file, ",")[1:]
entries = {i[0]: float(i[7]) for i in entries if float(i[7]) > 0}



motif_set = "int3"
motif_file = "source_data/motif_sets/{0}.txt".format(motif_set)
motifs = sequo.read_motifs(motif_file)


sim_path = "clean_run/dinucleotide_controls/{0}_dinucleotide_controls_matched_stops".format(motif_set)
sims = gen.get_filepaths(sim_path)[:100]


names, seqs = gen.read_fasta(file)
families = gen.read_many_fields(families_file, "\t")
seqs = {name.split(".")[0]: seqs[i] for i, name in enumerate(names)}

seqs = sequo.group_family_results(seqs, families)
seqs = {i: seqs[i] for i in seqs if i in entries}
# seqs = [seqs[i][0] for i in seqs]
# seqs = seqs[:0]


def remove_motifs(seqs, motifs):
    retained_seqs = []
    for seq in seqs:
        overlaps = sequo.sequence_overlap_indicies(seq, motifs)
        chunked_overlaps = sequo.chunk_overlaps(overlaps)

        # removed_indices = []
        # change_indices = []
        # for chunk in chunked_overlaps:
        #     if len(chunk) % 3 == 0:
        #         removed_indices.extend(chunk)
        #     else:
        #         required_to_maintain_frame = len(chunk) % 3
        #         change_indices.extend(chunk[:required_to_maintain_frame])
        #         removed_indices.extend(chunk[required_to_maintain_frame:])
        new_seq = []
        for pos, nt in enumerate(list(seq)):
            # if pos not in removed_indices:
            #     if pos in change_indices:
            #         new_seq.append("C")
            #     else:
            #         new_seq.append(nt)
            if pos in overlaps:
                new_seq.append("C")
            else:
                new_seq.append(nt)
        new_seq = "".join(new_seq)
        retained_seqs.append(new_seq)
    return retained_seqs



def randomise_seqs(seqs):
    new_seqs = []
    for seq in seqs:
        nts = list(seq)
        np.random.shuffle(nts)
        new_seqs.append("".join(nts))

    return new_seqs


# removed = remove_motifs(seqs, motifs)
# real_longest_orfs = seqo.get_longest_orfs(seqs)
# removed_longest_orfs = seqo.get_longest_orfs(removed)
#
# outputs = []
# for i, seq in enumerate(removed):
#     real = seqo.get_longest_orfs([seqs[i]])
#     ls = []
#     for i in range(100):
#         rseq = randomise_seqs([seq])
#         l = seqo.get_longest_orfs(rseq)
#         ls.append(l)
#     ls = gen.flatten(ls)
#     z = np.divide(real[0] - np.nanmean(ls), np.nanstd(ls))
#     outputs.append(z)
#
# print(len([i for i in outputs if i > 0]))
# longer = []
# for i, real in enumerate(real_longest_orfs):
#     if real > removed_longest_orfs[i]:
#         longer.append(i)
#
# print(len(longer))

def calc_overlaps(file):
    motifs = sequo.read_motifs(file)
    new_seqs = []
    for seq in seqs:
        overlaps = sequo.sequence_overlap_indicies(seq, motifs)
        new_seq = "".join([nt if i not in overlaps else "C" for i, nt in enumerate(list(seq))])
        new_seqs.append(new_seq)
    longest_orfs = seqo.get_longest_orfs(new_seqs)
    return longest_orfs




def mask_seq(seq, motifs):
    overlaps = sequo.sequence_overlap_indicies(seq, motifs)
    new_seq = "".join([nt if i not in overlaps else "C" for i, nt in enumerate(list(seq))])
    return new_seq


def calc_lengths(ids, seqs, motif_sets):
    sim_lengths = {}
    if len(seqs):
        for i, id in enumerate(ids):
            gen.print_parallel_status(i, ids)
            sim_lengths[id] = []
            seq = seqs[id][0]
            for motif_set in motif_sets:
                new_seq = mask_seq(seq, motif_set)
                longest_orf = seqo.get_longest_orfs([new_seq])[0]
                sim_lengths[id].append(longest_orf)

    return sim_lengths


motif_sets = [sequo.read_motifs(sim_file) for i, sim_file in enumerate(sims)]

real_masked = {i: mask_seq(seqs[i][0], motifs) for i in seqs}
real_longest_orfs = {i: seqo.get_longest_orfs([real_masked[i]])[0] for i in real_masked}

# print(real_longest_orfs)
sim_lengths = soc.run_simulation_function(list(seqs), [seqs, motif_sets], calc_lengths, sim_run = False)

zs = []
for id in real_longest_orfs:
    # print(id, real_longest_orfs[id], sim_lengths[id])
    z = np.divide(real_longest_orfs[id] - np.nanmean(sim_lengths[id]), np.nanmean(sim_lengths[id]))
    zs.append(z)

print(len([i for i in zs if i > 0]), len(zs))


# sim_orfs = []
# for i, file in enumerate(sims):
#     sim_orfs.append(calc_overlaps(file))
#
# zs = []
# for i, length in enumerate(real_longest_orfs):
#     sims = [sim[i] for sim in sim_orfs]
#     z = np.divide(length - np.nanmean(sims), np.nanstd(sims))
#     print(length, sims, np.nanmean(sims), z)
#     zs.append(z)
#
# print(real_longest_orfs)
#
# print(len([z for z in zs if z < 0]), len(zs))
# print(len([z for z in zs if z < 0 and z < -1.96]), len(zs))
