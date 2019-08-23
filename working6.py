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

motif_file = "source_data/motif_sets/int3.txt"

cds_fasta = "clean_run/genome_sequences/human/human.cds.clean.fasta"
single_exon_cds_fasta = "clean_run/genome_sequences/human/human.cds.single_exons.fasta"

exons_bed = "clean_run/genome_sequences/human/human.exons.bed"
filtered_exons_bed = "clean_run/genome_sequences/human/human.all_exons_from_transcripts.bed"
all_exons_fasta = "clean_run/genome_sequences/human/human.all_exons_from_transcripts.fasta"
genome_fasta = "../source_data/genomes/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
families_file = "clean_run/genome_sequences/human/human.cds.families.bed"

single_exon_blast_file = "clean_run/genome_sequences/human/human.cds.blast_all_against_all_single.csv"
single_exon_families = "clean_run/genome_sequences/human/human.cds.families_single.bed"


full_mature_transcripts = "clean_run/genome_sequences/human/human.multi_exon_full_transcripts.fasta"
single_exon_mature_transcripts = "clean_run/genome_sequences/human/human.single_exon_full_transcripts.fasta"
multi_utr = "clean_run/genome_sequences/human/human.multi_exons_five_utr.fasta"
single_utr = "clean_run/genome_sequences/human/human.single_exon_five_utr.fasta"


# fo.fasta_from_intervals(exons_bed, all_exons_fasta, genome_fasta, names = True)



def get_full_transcripts(cds_fasta, exons_fasta, output_file):

    cds_names, cds_seqs = gen.read_fasta(cds_fasta)
    cds_list = {name: cds_seqs[i] for i, name in enumerate(cds_names)}

    exon_names, exon_seqs = gen.read_fasta(exons_fasta)
    exon_list = collections.defaultdict(lambda: collections.defaultdict(lambda: []))

    for i, name in enumerate(exon_names):
        id = name.split(".")[0]
        exon_id = int(name.split(".")[1].split("(")[0])
        exon_list[id][exon_id] = exon_seqs[i]

    full_spliced_transcripts = {}
    for id in exon_list:
        exons = []
        for exon_id in sorted(exon_list[id]):
            exons.append(exon_list[id][exon_id])
        full_spliced_transcripts[id] = "".join(exons)

    with open(output_file, "w") as outfile:
        for id in full_spliced_transcripts:
            outfile.write(">{0}\n{1}\n".format(id, full_spliced_transcripts[id]))

np.random.seed()
def randomise(seq):
    nts = list(seq)
    np.random.shuffle(nts)
    return "".join(nts)


def randomise_densities(ids, seq_list, iterations):

    random_densities = collections.defaultdict(lambda: [])
    if len(ids):

        for i, id in enumerate(ids):
            print("ID {0}/{1}".format(i+1, len(ids)))

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


def get_utrs(transcript_fasta, cds_fasta, output_file):

    cds_names, cds_seqs = gen.read_fasta(cds_fasta)
    cds_list = {name: cds_seqs[i] for i, name in enumerate(cds_names)}

    transcript_names, transcript_seqs = gen.read_fasta(transcript_fasta)
    transcript_list = {name: transcript_seqs[i] for i, name in enumerate(transcript_names)}

    with open(output_file, "w") as outfile:
        for id in cds_list:
            cds = cds_list[id]
            transcript = transcript_list[id]
            try:
                cds_start_index = transcript.index(cds)
                if cds_start_index > 0:
                    outfile.write(">{0}\n{1}\n".format(id, transcript[:cds_start_index]))
            except:
                pass


def get_output(entry):
    output = []
    for id in entry:
        output.extend(entry[id])
    return gen.stringify(output)

# cons.filter_families(single_exon_cds_fasta, single_exon_blast_file, single_exon_families, database_path = None, clean_run = None)

# get_full_transcripts(cds_fasta, all_exons_fasta, full_mature_transcripts)
# get_full_transcripts(single_exon_cds_fasta, all_exons_fasta, single_exon_mature_transcripts)

# get_utrs(full_mature_transcripts, cds_fasta, multi_utr)
# get_utrs(single_exon_mature_transcripts, single_exon_cds_fasta, single_utr)
#



single_cds_names, single_cds_seqs = gen.read_fasta(single_utr)
single_cds_list = {name: [single_cds_seqs[i]] for i, name in enumerate(single_cds_names) if len(single_cds_seqs[i]) > 50}
single_cds_list = sequo.pick_random_family_member(single_exon_families, single_cds_list)
# single_densities, single_nds = calc_values(single_cds_list)

# 473

cds_names, cds_seqs = gen.read_fasta(cds_fasta)
cds_list = {name: cds_seqs[i] for i, name in enumerate(cds_names)}

# print(list(cds_list))

multi_cds_names, multi_cds_seqs = gen.read_fasta(multi_utr)
multi_cds_list = {name: [multi_cds_seqs[i]] for i, name in enumerate(multi_cds_names) if len(multi_cds_seqs[i]) > 50 and name in cds_list}
multi_cds_list = sequo.pick_random_family_member(families_file, multi_cds_list)
#
# 6086
# multi_densities, multi_nds = calc_values(multi_cds_list)


motifs = sequo.read_motifs(motif_file)
multi_seqs = [multi_cds_list[i][0] for i in multi_cds_list]
single_seqs = [single_cds_list[i][0] for i in single_cds_list]
single_ese_desities = [seqo.calc_motif_density([i], motifs) for i in single_seqs]
multi_ese_densities = [seqo.calc_motif_density([i], motifs) for i in multi_seqs]

print(single_ese_desities)

output_file = "temp_data/utr_ese_densities.csv"
with open(output_file, "w") as outfile:
    outfile.write("single,{0}\n".format(",".join(gen.stringify(single_ese_desities))))
    outfile.write("multi,{0}\n".format(",".join(gen.stringify(multi_ese_densities))))


# output_file = "temp_data/utr_densities.csv"
#

#
# with open(output_file, "w") as outfile:
#
#     outfile.write("single_density,{0}\n".format(",".join(get_output(single_densities))))
#     outfile.write("single_nd,{0}\n".format(",".join(get_output(single_nds))))
#     outfile.write("multi_density,{0}\n".format(",".join(get_output(multi_densities))))
#     outfile.write("multi_nd,{0}\n".format(",".join(get_output(multi_nds))))
