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

def randomise(seq):
    nts = list(seq)
    np.random.shuffle(nts)
    return "".join(nts)


motif_sets = ["int3.txt", "ESR.txt", "ke400.txt", "PESE.txt", "RESCUE.txt", "combined_eses.txt"]

for set in motif_sets:
    motifs = sequo.read_motifs("source_data/motif_sets/{0}".format(set))
    output = [set]
    [output.append(seqo.calc_motif_density(motifs, [stop])) for stop in sorted(stops)]
    print("\t".join(gen.stringify(output)))

motifs_file = "source_data/motif_sets/combined_eses.txt"

cabili_seqs = "clean_run/genome_sequences/lincrna/cabili/transcript_sequences.fasta"
# cabili_seqs = "clean_run/genome_sequences/lincrna/GENCODE_CLS/hsAllCompmerge.cage+polyASupported.full_transcripts.fasta"
families_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt"
# families_file = "clean_run/genome_sequences/lincrna/GENCODE_CLS/families.txt"
# families_file1 = "clean_run/genome_sequences/lincrna/GENCODE_CLS/families.txt"

# entries = gen.read_many_fields(families_file, "\t")
# with open(families_file1, "w") as outfile:
#     for entry in entries:
#         outfile.write("{0}\n".format("\t".join([".".join(i.split(":")) for i in entry])))


def sim_removed(simulations, sequence_list, motifs, inverse = None):
    """
    Calculate stop codon densities in lincRNA sequences and simulations

    Args:
        simulations (list): list containing simulation to iterate over
        sequence_list (dict): sequences
        output_directory (str): path to output directory

    Returns:
        outputs (list): list of output files
    """

    new_outputs = collections.defaultdict(lambda: [])

    if simulations:
        np.random.seed()
        for i, simulation_id in enumerate(simulations):
            gen.print_parallel_status(i, simulations)
            # check if the output file doesnt already exist

            hits = []
            for id in sequence_list:
                sequence = sequence_list[id]
                overlaps = sequo.sequence_overlap_indicies(sequence, motifs)
                non_overlaps = [i for i in list(range(len(sequence_list[id]))) if i not in overlaps]
                chunked_overlaps = sequo.chunk_overlaps(non_overlaps)
                overlap_motifs = ["".join([sequence_list[id][i] for i in chunk]) for chunk in chunked_overlaps]
                hits.extend(overlap_motifs)

            if simulation_id != "real":
                hit_lengths = [len(i) for i in hits]
                all_nts = list("".join(hits))
                np.random.shuffle(all_nts)
                nt_index = 0
                new_hits = []
                for length in hit_lengths:
                    new_hits.append("".join(all_nts[nt_index:nt_index+length]))
                    nt_index += length
                hits = new_hits

            for stop in sorted(stops):
                density = seqo.calc_motif_density(hits, [stop])
                new_outputs[simulation_id].append(density)

    outputs = {}
    for id in new_outputs:
        outputs[id] = new_outputs[id]

    return outputs


def randomise_seqs(seqs):
    sim_seqs = []
    for seq in seqs:
        nts = list(seq)
        np.random.shuffle(nts)
        sim_seqs.append("".join(nts))
    return sim_seqs


motifs = sequo.read_motifs(motifs_file)

names, seqs = gen.read_fasta(cabili_seqs)
seq_list = {id: seqs[i] for i, id in enumerate(names)}

seq_list = sequo.pick_random_family_member(families_file, seq_list)

seqs = [seq_list[id] for id in seq_list]
real_densities = {stop: seqo.calc_motif_density(seqs, [stop]) for stop in stops}


# sim_densities = collections.defaultdict(lambda: [])
# for i in range(50):
#     print(i+1)
#     sim_seqs = randomise_seqs(seqs)
#     for stop in stops:
#         sim_densities[stop].append(seqo.calc_motif_density(sim_seqs, [stop]))
#
# for stop in stops:
#     real = real_densities[stop]
#     sims = sim_densities[stop]
#     nd = np.divide(real - np.mean(sims), np.mean(sims))
#     print(stop, real, nd)

# simulations = ["real"] + list(range(10))
# outputs = soc.run_simulation_function(simulations, [seq_list, motifs], sim_removed, sim_run = False)
# print(outputs)
# # outputs = sim_removed(simulations, seq_list, motifs)
#
# for i, stop in enumerate(sorted(stops)):
#     real = outputs["real"][i]
#     sims = [outputs[it][i] for it in outputs if it != "real"]
#     nd = np.divide(real - np.mean(sims), np.mean(sims))
#     print(stop, real, nd)
