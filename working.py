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
# import _rs_conservation as rsc
import os
import zipfile
from useful_motif_sets import stops
from regex_patterns import codon_pattern



def densities():

    prefix = "clean_run/genome_sequences/human/human"
    cds_file = "{0}.cds.clean.fasta".format(prefix)
    coding_exons_file = "{0}.cds.clean_coding_exons.fasta".format(prefix)
    coding_exons_introns_file = "{0}.clean_introns.fasta".format(prefix)
    families_file = "{0}.cds.families.bed".format(prefix)

    cds_list = gen.fasta_to_list(cds_file)
    coding_exons = gen.fasta_to_list(coding_exons_file, split = "(")
    coding_exons_list = collections.defaultdict(lambda: collections.defaultdict())

    for i in coding_exons:
        coding_exons_list[i.split(".")[0]][int(i.split(".")[1])] = coding_exons[i]


    # cds_joined_exons = {}
    # for i in coding_exons_list:
    #     exons = coding_exons_list[i]
    #     cds_joined_exons[i] = "".join([exons[i] for i in sorted(exons)])
    #
    # # exon_stop_densities = {i: seqo.calc_seqs_codon_set_density([cds_joined_exons[i]], codon_set = stops) for i in cds_joined_exons}
    # #
    # # print(exon_stop_densities)
    #

    # coding_exons = {i: coding_exons[i] for j, i in enumerate(coding_exons) if j < 50}

    exon_densities = collections.defaultdict(lambda: [])
    for i, exon in enumerate(coding_exons):
        print("{0}/{1}".format(i, len(coding_exons)))
        cds_seq = cds_list[exon.split(".")[0]]
        exon_seq = coding_exons[exon]
        exon_start_index = cds_seq.index(exon_seq)
        exon_start_frame = exon_start_index % 3
        if exon_start_frame == 1:
            exon_seq = "N" + exon_seq
        elif exon_start_frame == 2:
            exon_seq = "NN" + exon_seq
        exon_densities[exon.split(".")[0]].append(seqo.calc_seqs_codon_set_density([exon_seq], codon_set = stops, exclude_frames = [0]))

    introns = gen.fasta_to_list(coding_exons_introns_file, split = "(")

    intron_list = collections.defaultdict(lambda: collections.defaultdict())
    for intron in introns:
        splits = intron.split(".")
        intron_list[splits[0]][splits[1]] = introns[intron]

    intron_densities = collections.defaultdict(lambda: [])
    for i, id in enumerate(exon_densities):
        print("{0}/{1}".format(i, len(exon_densities)))
        introns = intron_list[id]
        for intron in introns:
            exon_id = "{0}.{1}".format(id, intron.split("-")[0])
            exon = coding_exons[exon_id]
            cds = cds_list[id]
            exon_end_frame = (cds.index(exon) + len(exon)) % 3

            intron_seq = introns[intron]
            if exon_end_frame == 0:
                intron_seq = "NN" + intron_seq
            elif exon_end_frame == 1:
                intron_seq = "N" + intron_seq

            density = seqo.calc_intron_seqs_stop_density([intron_seq])

            intron_densities[id].append(density)
            # intron_densities[id].append(np.mean([density, density1, density2]))

    exon_densities = {i: np.median(exon_densities[i]) for i in exon_densities}
    intron_densities = {i: np.median(intron_densities[i]) for i in intron_densities}

    if families_file:
        families = gen.read_many_fields(families_file, "\t")
        exon_densities = sequo.group_family_results(exon_densities, families)
        intron_densities = sequo.group_family_results(intron_densities, families)

    with open("temp_files/densities.csv", "w") as outfile:
        outfile.write("id,exon,intron\n")
        for id in exon_densities:
            if id in intron_densities:
                outfile.write("{0},{1},{2}\n".format(id, np.median(exon_densities[id]), np.median(intron_densities[id])))

def lincRNA_densities():

    prefix = "clean_run/genome_sequences/human/human"
    cds_file = "{0}.cds.clean.fasta".format(prefix)
    coding_exons_file = "{0}.cds.clean_coding_exons.fasta".format(prefix)
    coding_exons_introns_file = "{0}.clean_introns.fasta".format(prefix)
    families_file = "{0}.cds.families.bed".format(prefix)
    ese_file = "source_data/motif_sets/int3.txt"

    lncrna_exons = "clean_run/lincrna/lincRNA.multi_exon.exons.fasta"
    lncrna_introns = "clean_run/lincrna/lincRNA.multi_exon.introns.fasta"
    lncrna_families = "clean_run/lincrna/lincRNA.multi_exon_families.bed"

    # cds_list = gen.fasta_to_list(cds_file)
    # coding_exons = gen.fasta_to_list(coding_exons_file, split = "(")
    # coding_exons_list = collections.defaultdict(lambda: collections.defaultdict())
    #
    # for i in coding_exons:
    #     coding_exons_list[i.split(".")[0]][int(i.split(".")[1])] = coding_exons[i]


    exons = gen.fasta_to_list(lncrna_exons, split = "(")
    exon_list = collections.defaultdict(lambda: [])
    [exon_list[i.split(".")[0]].append(exons[i]) for i in exons]

    exon_densities = {}
    for id in exon_list:
        seqs = exon_list[id]
        id_densities = []
        for seq in seqs:
            id_densities.append(seqo.calc_seqs_codon_set_density([seq], codon_set = stops))
        exon_densities[id] = np.median(id_densities)

    introns = gen.fasta_to_list(lncrna_introns, split = "(")
    intron_list = collections.defaultdict(lambda: [])
    [intron_list[i.split(".")[0]].append(introns[i]) for i in introns]

    intron_densities = {}
    for id in exon_densities:
        seqs = intron_list[id]
        id_densities = []
        for seq in seqs:
            id_densities.append(seqo.calc_seqs_codon_set_density([seq], codon_set = stops))
        intron_densities[id] = np.median(id_densities)

    families = gen.read_many_fields(lncrna_families, "\t")
    exon_densities = sequo.group_family_results(exon_densities, families)
    intron_densities = sequo.group_family_results(intron_densities, families)

    with open("temp_files/lincRNA.csv", "w") as outfile:
        outfile.write("id,exon,intron\n")
        [outfile.write("{0},{1},{2}\n".format(id, np.median(exon_densities[id]), np.median(intron_densities[id]))) for id in exon_densities]


densities()
