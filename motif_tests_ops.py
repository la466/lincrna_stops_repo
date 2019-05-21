import generic as gen
import containers as cont
import sim_ops_containers as simopc
import main_tests_ops as mto
import ops_containers as opsc
import seq_ops as seqo
import sequence_ops as sequo
import file_ops as fo
import time
from useful_motif_sets import stops
import os
import random
import collections
import numpy as np
import multiprocessing as mp



def motif_stop_codon_densities(motif_file, motif_controls_directory, required_simulations, output_file):

    filelist = {"real": motif_file}
    for i, file in enumerate(os.listdir(motif_controls_directory)[:required_simulations]):
        filelist[i] = "{0}/{1}".format(motif_controls_directory, file)
    file_ids = [i for i in filelist]

    temp_dir = "temp_motif_dir"
    gen.create_output_directories(temp_dir)

    args = [filelist, temp_dir]
    outputs = simopc.run_simulation_function(file_ids, args, calculate_stop_codon_densities, sim_run = False)

    with open(output_file, "w") as outfile:
        outfile.write("sim_id,stop_density\n")
        for file in outputs:
            outfile.write("{0}\n".format(",".join(gen.read_many_fields(file, "\t")[0])))



def calculate_stop_codon_densities(file_ids, filelist, output_directory):

    outputs = []

    if len(file_ids):
        for i, id in enumerate(file_ids):
            print("(W{0}) {1}/{2}".format(mp.current_process().name.split("-")[-1], i+1, len(file_ids)))
            motifs = [i[0] for i in gen.read_many_fields(filelist[id], "\t") if "#" not in i[0] and ">" not in i[0]]
            density = seqo.calc_seqs_codon_set_density(motifs, codon_set = stops)

            if "real" in str(id):
                temp_file = "{0}/real.txt".format(output_directory)
            else:
                temp_file = "{0}/{1}.txt".format(output_directory, random.random())
            outputs.append(temp_file)

            with open(temp_file, "w") as temp:
                temp.write("{0}\t{1}\n".format(id, density))

    return outputs


def motif_codon_densities(motif_file, codon_combinations_file, motif_controls_directory, required_simulations, output_file):

    filelist = {"real": motif_file}
    for i, file in enumerate(os.listdir(motif_controls_directory)[:required_simulations]):
        filelist[i] = "{0}/{1}".format(motif_controls_directory, file)
    file_ids = [i for i in filelist]


    codon_sets = gen.read_many_fields(codon_combinations_file, "\t")


    temp_dir = "temp_motif_dir"
    gen.create_output_directories(temp_dir)
    args = [filelist, codon_sets, temp_dir]
    outputs = simopc.run_simulation_function(file_ids, args, calculate_motif_densities, sim_run = False)

    real_density_list = {}
    sim_density_list = collections.defaultdict(lambda: [])

    for file in outputs:
        results = gen.read_many_fields(file, "\t")
        if "real" in file:
            for i in results:
                real_density_list[i[0]] = float(i[1])
        else:
            for i in results:
                sim_density_list[i[0]].append(float(i[1]))

    with open(output_file, "w") as outfile:
        outfile.write("codons,gc_content,purine_content,density,nd\n")
        for codon_set in sorted(real_density_list):
            nd = np.divide(real_density_list[codon_set] - np.mean(sim_density_list[codon_set]), np.mean(sim_density_list[codon_set]))
            outputs = [codon_set, seqo.calc_gc_seqs_combined(codon_set.split("_")), sequo.calc_purine_content(codon_set.split("_")), real_density_list[codon_set], nd]
            outfile.write("{0}\n".format(",".join(gen.stringify(outputs))))

    gen.remove_directory(temp_dir)


def calculate_motif_densities(file_ids, filelist, codon_sets, output_directory):

    outputs = []

    if len(file_ids):
        for i, id in enumerate(file_ids):
            print("(W{0}) {1}/{2}".format(mp.current_process().name.split("-")[-1], i+1, len(file_ids)))
            motifs = [i[0] for i in gen.read_many_fields(filelist[id], "\t") if "#" not in i[0] and ">" not in i[0]]
            density_list = {}
            for set in codon_sets:
                density = seqo.calc_seqs_codon_set_density(motifs, codon_set = set)
                density_list["_".join(sorted(set))] = density

            if "real" in str(id):
                temp_file = "{0}/real.txt".format(output_directory)
            else:
                temp_file = "{0}/{1}.txt".format(output_directory, random.random())
            outputs.append(temp_file)

            with open(temp_file, "w") as temp:
                [temp.write("{0}\t{1}\n".format(i, density_list[i])) for i in sorted(density_list)]

    return outputs


def calc_stop_densities(input_file):
    """
    Calculate the stop density in the set of motifs and other hexamers of same length.

    Args:
        input_file (str): path to input file
    """

    motifs = sequo.read_motifs(input_file)
    motif_lengths = list(set([len(i) for i in motifs]))

    all_motifs = []
    for length in motif_lengths:
        all_motifs.extend(["".join(j) for j in it.product(nucleotides, repeat = length)])

    non_motifs = [i for i in all_motifs if i not in motifs]

    print("Stop density in motifs: {0}".format(seqo.calc_gc_seqs_combined(motifs)))
    print("Stop density in non motifs {0}".format(seqo.calc_gc_seqs_combined(non_motif)))
