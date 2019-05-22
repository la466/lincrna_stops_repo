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

motifs = sequo.read_motifs("source_data/motif_sets/combined_eses.txt")

def remove_motifs(seqs, motifs):
    kept = {}
    for id in seqs:
        kept_seqs = []
        for seq in seqs[id]:
            overlaps = sequo.get_motifs_overlap_indices([seq], motifs)
            kept_seq = "".join([nt for i, nt in enumerate(list(seq)) if i not in overlaps])
            kept_seqs.append(kept_seq)
        kept[id] = kept_seqs
    return kept

exons = {i: [exons[i][j] for j in exons[i]] for i in exons}
introns = {i: [introns[i][j] for j in introns[i]] for i in introns}

exons = remove_motifs(exons, motifs)
introns = remove_motifs(introns, motifs)

exon_purine = {i: sequo.calc_purine_content(exons[i]) for i in exons}
intron_purine = {i: sequo.calc_purine_content(introns[i]) for i in introns}

families = gen.read_many_fields(families_file, "\t")
exon_purine = sequo.group_family_results(exon_purine, families)
intron_purine = sequo.group_family_results(intron_purine, families)
#
with open('clean_run/tests/purine_content/exon_intron_purine_no_eses.csv', "w") as outfile:
    outfile.write("id,exon_purine,intron_purine\n")
    [outfile.write("{0},{1},{2}\n".format(id, np.nanmedian(exon_purine[id]), np.nanmedian(intron_purine[id]))) for id in exon_purine]
