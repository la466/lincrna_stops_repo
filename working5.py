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


ese_file = "source_data/motif_sets/INT3.txt"
eses = sequo.read_motifs(ese_file)

output_directory = "clean_run/tests/ese_densities"
gen.create_output_directories(output_directory)

fileset = {
    "fileset_prefix": "pc",
    "exons_file": "clean_run/genome_sequences/human/human.cds.clean_coding_exons.fasta",
    "introns_file": "clean_run/genome_sequences/human/human.clean_introns.fasta",
    "families_file": "clean_run/genome_sequences/human/human.cds.families.bed" ,
}

lincrna_files = {
    "fileset_prefix": "lincrna",
    "exons_file": "clean_run/genome_sequences/lincrna/cabili/exons.fasta",
    "introns_file": "clean_run/genome_sequences/lincrna/cabili/introns.fasta",
    "families_file": "clean_run/genome_sequences/human/human.cds.families.bed" ,
}



exon_names, exon_seqs = gen.read_fasta(fileset["exons_file"])
exon_list = collections.defaultdict(lambda: collections.defaultdict())
for i, name in enumerate(exon_names[:1000]):
    if len(exon_seqs[i]) >= 100:
        exon_list[name.split(".")[0]][int(name.split(".")[-1].split("(")[0])] = exon_seqs[i]


intron_names, intron_seqs = gen.read_fasta(fileset["introns_file"])
intron_list = collections.defaultdict(lambda: collections.defaultdict())
for i, name in enumerate(intron_names):
    id = name.split(".")[0]
    flank1 = int(name.split(".")[-1].split("-")[0])
    if len(intron_seqs[i]) >= 400:
        intron_list[id][flank1] = intron_seqs[i]

temp_exon_list = collections.defaultdict(lambda: collections.defaultdict())
temp_intron_list = collections.defaultdict(lambda: collections.defaultdict())

for id in exon_list:
    for exon_id in sorted(exon_list[id]):
        if exon_id in intron_list[id] and exon_id - 1 in intron_list[id]:
            temp_exon_list[id][exon_id] = exon_list[id][exon_id]

exon_list = {i: {j: temp_exon_list[i][j] for j in temp_exon_list[i]} for i in temp_exon_list}

for id in exon_list:
    for exon_id in exon_list[id]:
        temp_intron_list[id][exon_id] = intron_list[id][exon_id]
        temp_intron_list[id][exon_id-1] = intron_list[id][exon_id-1]

intron_list = {i: {j: temp_intron_list[i][j] for j in temp_intron_list[i]} for i in temp_intron_list}

exons = []
[[exons.append(exon_list[id][exon_id]) for exon_id in exon_list[id]] for id in exon_list]
introns = []
[[introns.append(intron_list[id][intron_id]) for intron_id in intron_list[id]] for id in intron_list]

exon_five_prime_flanks = [seq[5:55] for seq in exons]
intron_three_prime_flanks = [seq[-55:-5] for seq in exons]

intron_five_prime_flanks = [seq[5:155] for seq in introns]
intron_three_prime_flanks = [seq[-170:-20] for seq in introns]


# motifs = ["AGAAGA", "GGAAGA", "GAAGAA", "TATTTT", "TTTTTT", "TTATTT"]
motifs = ["AGAAGA", "GGAAGA", "GAAGAA"]
# motifs = ["AGAAGA"]




def delta_calc(exon_no, intron_no, exon_frequency, intron_frequency):

    diff = exon_frequency - intron_frequency
    inv_exon_no = np.divide(1, exon_no)

    inv_intron_no = np.divide(1, intron_no)

    g = np.divide((exon_no * exon_frequency) + (intron_no + intron_frequency), (exon_no + intron_no))

    d = np.divide(diff, np.sqrt((inv_exon_no + inv_intron_no)*g))
    return(d)


def calc_delta(motif,exons,introns):

    exon_hits = sum([len([i for i in re.finditer("(?={0})".format(motif), seq)]) for seq in exons])
    intron_hits = sum([len([i for i in re.finditer("(?={0})".format(motif), seq)]) for seq in introns])

    # exon_overlaps = [sequo.return_overlap_motifs(i, [motif]) for i in exons]
    # intron_overlaps = [sequo.return_overlap_motifs(i, [motif]) for i in introns]

    exon_lengths = sum([len(i) for i in exons])
    intron_lengths = sum([len(i) for i in introns])

    # exon_overlap_nts =  sum([len(i) for i in gen.flatten(exon_overlaps)])
    # intron_overlap_nts =  sum([len(i) for i in gen.flatten(intron_overlaps)])

    exon_motif_frequency = np.divide(exon_hits, exon_lengths)
    intron_motif_frequency = np.divide(intron_hits, intron_lengths)

    exon_no = len(exons)
    intron_no = len(exons)

    d = delta_calc(exon_no, intron_no, exon_motif_frequency, intron_motif_frequency)
    return d



for motif in motifs:

    d = calc_delta(motif, exon_five_prime_flanks, intron_three_prime_flanks)

    print(motif, d)


# seq = "TAGAAGATAAAGAAGTAAGAAGAAGAAGATTTGAGA"
#
# test =
