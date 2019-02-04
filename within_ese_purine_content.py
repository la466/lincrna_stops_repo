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
families_file = "clean_run/genome_sequences/human/human.cds.families.bed"


ese_file = "source_data/motif_sets/int3.txt"
eses = sequo.read_motifs(ese_file)

exon_names, exon_seqs = gen.read_fasta(exons_file)

exons = collections.defaultdict(lambda: [])
for i, name in enumerate(exon_names):
    if len(exon_seqs[i]) > 100:
        exons[name.split(".")[0]].append(exon_seqs[i])
exons = {i: exons[i] for i in exons}
exons = {i: exons[i] for i in exons if len(exons[i]) > 1}
families = gen.read_many_fields(families_file, "\t")


def calc_purine_content(ids, exons):

    outputs = {}

    if len(ids):
        for i, id in enumerate(ids):
            gen.print_parallel_status(i, ids)

            ese_sequence = []
            non_ese_sequence = []
            flank_sequence_ese = []
            flank_sequence_non_ese = []
            for sequence in exons[id]:
                overlaps = sequo.sequence_overlap_indicies(sequence, eses)
                ese_sequence.append("".join([sequence[i] for i in overlaps]))
                non_ese_sequence.append("".join([sequence[i] for i in range(len(sequence)) if i not in overlaps]))

                flank_sequence = "".join([sequence[:50], sequence[-50:]])
                flank_sequence_ese.append("".join([flank_sequence[i] for i in overlaps if i in range(len(flank_sequence))]))
                flank_sequence_non_ese.append("".join([flank_sequence[i] for i in range(len(flank_sequence)) if i not in overlaps]))

            outputs[id] = [sequo.calc_purine_content("".join(ese_sequence)), sequo.calc_purine_content("".join(non_ese_sequence)), sequo.calc_purine_content("".join(flank_sequence_ese)), sequo.calc_purine_content("".join(flank_sequence_non_ese))]

    return outputs

outputs = soc.run_simulation_function(list(exons), [exons], calc_purine_content, sim_run = False)

ese_purine = {}
non_ese_purine = {}
flank_sequence_ese = {}
flank_sequence_non_ese = {}

for id in outputs:
    ese_purine[id] = outputs[id][0]
    non_ese_purine[id] = outputs[id][1]
    flank_sequence_ese[id] = outputs[id][2]
    flank_sequence_non_ese[id] = outputs[id][3]

ese_purine = sequo.group_family_results(ese_purine, families)
non_ese_purine = sequo.group_family_results(non_ese_purine, families)
flank_sequence_ese = sequo.group_family_results(flank_sequence_ese, families)
flank_sequence_non_ese = sequo.group_family_results(flank_sequence_non_ese, families)


output_file = "clean_run/tests/purine_content/exon_ese_non_ese_purine_content.csv"
with open(output_file, "w") as outfile:
    outfile.write("id,ese_purine,non_ese_purine,flanking_50_ese_purine,flanking_50_non_ese_purine\n")
    [outfile.write("{0},{1},{2},{3},{4}\n".format(id, np.nanmedian(ese_purine[id]), np.nanmedian(non_ese_purine[id]), np.nanmedian(flank_sequence_ese[id]), np.nanmedian(flank_sequence_non_ese[id]))) for id in ese_purine]
