import generic as gen
import containers as cont
import file_ops as fo
import seq_ops as seqo
import sequence_ops as sequo
import time
import random
import numpy as np
import collections

def coding_exons(input_file, families_file, output_directory):

    output_directory = "{0}/coding_exons".format(output_directory)
    gen.create_output_directories(output_directory)

    # get the coding exons
    coding_exon_names, coding_exon_seqs = gen.read_fasta(input_file)
    coding_exons = {name.split("(")[0]: coding_exon_seqs[i] for i, name in enumerate(coding_exon_names) if len(coding_exon_seqs[i]) >= 210}

    coding_exon_five_prime = {i: coding_exons[i][:69] for i in coding_exons}
    coding_exon_three_prime = {i: coding_exons[i][-69:] for i in coding_exons}

    coding_exon_cores = {}
    for i in coding_exons:
        middle = int(len(coding_exons[i]) / 2)
        coding_exon_cores[i] = coding_exons[i][middle-34:middle+35]

    five_gc = collections.defaultdict(lambda: [])
    three_gc = collections.defaultdict(lambda: [])
    core_gc = collections.defaultdict(lambda: [])

    [five_gc[i.split(".")[0]].append(seqo.calc_seq_gc(coding_exon_five_prime[i])) for i in coding_exon_five_prime]
    [three_gc[i.split(".")[0]].append(seqo.calc_seq_gc(coding_exon_three_prime[i])) for i in coding_exon_three_prime]
    [core_gc[i.split(".")[0]].append(seqo.calc_seq_gc(coding_exon_cores[i])) for i in coding_exon_cores]

    families = gen.read_many_fields(families_file, "\t")

    # group the densities by transcript id
    coding_exon_five_prime_list = collections.defaultdict(lambda: [])
    coding_exon_three_prime_list = collections.defaultdict(lambda: [])
    coding_exon_core_list = collections.defaultdict(lambda: [])

    [coding_exon_five_prime_list[i.split(".")[0]].append(seqo.calc_stop_density(coding_exon_five_prime[i])) for i in coding_exon_five_prime]
    [coding_exon_three_prime_list[i.split(".")[0]].append(seqo.calc_stop_density(coding_exon_three_prime[i])) for i in coding_exon_three_prime]
    [coding_exon_core_list[i.split(".")[0]].append(seqo.calc_stop_density(coding_exon_cores[i])) for i in coding_exon_cores]

    coding_exon_five_prime_list = {i: np.mean(coding_exon_five_prime_list[i]) for i in coding_exon_five_prime_list}
    coding_exon_three_prime_list = {i: np.mean(coding_exon_three_prime_list[i]) for i in coding_exon_three_prime_list}
    coding_exon_core_list = {i: np.mean(coding_exon_core_list[i]) for i in coding_exon_core_list}

    # group the densities by family
    coding_exon_five_prime = sequo.group_family_results(coding_exon_five_prime_list, families)
    coding_exon_three_prime = sequo.group_family_results(coding_exon_three_prime_list, families)
    coding_exon_cores = sequo.group_family_results(coding_exon_core_list, families)

    five_gc = {i: np.mean(five_gc[i]) for i in five_gc}
    three_gc = {i: np.mean(three_gc[i]) for i in three_gc}
    core_gc = {i: np.mean(core_gc[i]) for i in core_gc}



    five_gc = sequo.group_family_results(five_gc, families)
    three_gc = sequo.group_family_results(three_gc, families)
    core_gc = sequo.group_family_results(core_gc, families)

    output_file = "{0}/exon_region_densities.csv".format(output_directory)
    with open(output_file, "w") as outfile:
        outfile.write("id,start_density,core_density,end_density,start_gc,core_gc,end_gc\n")
        for i in sorted(coding_exon_five_prime):
            outfile.write("{0},{1},{2},{3},{4},{5},{6}\n".format(i, np.mean(coding_exon_five_prime[i]), np.mean(coding_exon_cores[i]), np.mean(coding_exon_three_prime[i]), np.mean(five_gc[i]), np.mean(core_gc[i]), np.mean(three_gc[i])))
