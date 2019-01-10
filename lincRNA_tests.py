import generic as gen
import sim_ops as simo
import sim_ops_containers as simoc
import seq_ops as seqo
import sequence_ops as sequo
import os
import collections
from useful_motif_sets import stops
import numpy as np
from progressbar import ProgressBar
import random
# TODO: simulation to get the normalised densities of lincRNA stop codons
# TODO: compare stop codon density in exons of lincRNA

def stop_density_test(input_fasta, output_file, required_simulations, families_file = None):
    simoc.simulate_sequence_stop_density(input_fasta, output_file, required_simulations, families_file = families_file)

def generate_dinucleotide_matched_seqs(simulations, seqs, names, dinucleotide_content, nucleotide_content, output_directory, seeds=None, seq_seeds=None, match_density = None):
    dinucleotide_choices = [dn for dn in sorted(dinucleotide_content)]
    dinucleotide_probabilities = [dinucleotide_content[dn] for dn in sorted(dinucleotide_content)]
    nucleotide_choices = [n for n in sorted(nucleotide_content)]
    nucleotide_probabilities = [nucleotide_content[n] for n in sorted(nucleotide_content)]

    output_files = []

    pbar = ProgressBar()

    if len(simulations):
        for sim_no, simulation in enumerate(pbar(simulations)):
            # set the seed
            simo.set_seed(seeds, simulation)

            # print the simulation number out
            # print("(W{0}) {1}/{2}: {3}".format(mp.current_process().name.split("-")[-1], sim_no+1, len(simulations), id))

            output_file = "{0}/{1}.txt".format(output_directory, random.random())
            output_files.append(output_file)
            with open(output_file, "w") as outfile:

                simulated_seqs = []
                # Generate a list of nucleotide matched sequences
                for i, seq in enumerate(seqs):
                    generated_seq = False
                    seq_seed = simo.get_seq_seed(seq_seeds, sim_no, i)
                    while not generated_seq:
                        sim_seq = seqo.generate_nt_matched_seq(seq, dinucleotide_choices, dinucleotide_probabilities, nucleotide_choices, nucleotide_probabilities, seed=seq_seed)
                        if sim_seq not in simulated_seqs:
                            generated_seq = True
                            simulated_seqs.append(sim_seq)

                            outfile.write(">{0}\n{1}\n".format(names[i],sim_seq))
    return output_files


def get_dinucleotide_matched_contols(input_file, simulations, output_directory):

    names, seqs = gen.read_fasta(input_file)
    nucleotide_content = seqo.get_nucleotide_content(seqs)
    dinucleotide_content = seqo.get_dinucleotide_content(seqs)
    generate_dinucleotide_matched_seqs(simulations, seqs, names, dinucleotide_content, nucleotide_content, output_directory)

def calculate_densities(filepath, families = None):

    names, seqs = gen.read_fasta(filepath)
    seq_list = collections.defaultdict(lambda: [])
    [seq_list[name.split(".")[0]].append(seqs[i]) for i, name in enumerate(names)]

    density = {i: seqo.calc_motif_density(seq_list[i], stops) for i in seq_list}

    if families:
        density = sequo.group_family_results(density, families)
    return density


def calc_gcs(filepath, families = None):
    names, seqs = gen.read_fasta(filepath)
    seq_list = collections.defaultdict(lambda: [])
    [seq_list[name.split(".")[0]].append(seqs[i]) for i, name in enumerate(names)]

    gcs = {i: seqo.calc_gc_seqs_combined(seq_list[i]) for i in seq_list}

    if families:
        gcs = sequo.group_family_results(gcs, families)
    return gcs


def density_simulation(exons_fasta, introns_fasta, simulations, families_file = None):

    # required = list(range(simulations))

    exon_output_directory = "clean_run/lincrna/exon_dinucleotide_controls"
    gen.create_output_directories(exon_output_directory)
    if len(os.listdir(exon_output_directory)) < simulations:
        required = list(range(simulations - len(os.listdir(exon_output_directory))))
        get_dinucleotide_matched_contols(exons_fasta, required, exon_output_directory)

    intron_output_directory = "clean_run/lincrna/intron_dinucleotide_controls"
    gen.create_output_directories(intron_output_directory)
    if len(os.listdir(intron_output_directory)) < simulations:
        required = list(range(simulations - len(os.listdir(intron_output_directory))))
        get_dinucleotide_matched_contols(introns_fasta, required, intron_output_directory)


    families = gen.read_many_fields(families_file, "\t")
    exons_density = calculate_densities(exons_fasta, families)
    introns_density = calculate_densities(introns_fasta, families)
    exons_gcs = calc_gcs(exons_fasta, families)
    introns_gcs = calc_gcs(introns_fasta, families)

    output_file = "clean_run/lincrna/tests/stop_densities.csv"
    with open(output_file, "w") as outfile:
        outfile.write("id,exon,intron,exon_gc,intron_gc\n")
        for i in exons_density:
            outfile.write("{0},{1},{2},{3},{4}\n".format(i, np.median(exons_density[i]), np.median(introns_density[i]), np.median(exons_gcs[i]), np.median(introns_gcs[i])))

    exon_sims = os.listdir(exon_output_directory)[:simulations]
    intron_sims = os.listdir(intron_output_directory)[:simulations]


    exons_density = calculate_densities(exons_fasta)
    introns_density = calculate_densities(introns_fasta)

    sim_exon_densities = collections.defaultdict(lambda: [])
    sim_intron_densities = collections.defaultdict(lambda: [])
    sim_exon_gcs = {}
    sim_intron_gcs = {}

    for i, file in enumerate(exon_sims):

        exon_filepath = "{0}/{1}".format(exon_output_directory, file)
        intron_filepath = "{0}/{1}".format(intron_output_directory, intron_sims[i])

        sim_exon_density = calculate_densities(exon_filepath)
        sim_intron_density = calculate_densities(intron_filepath)

        [sim_exon_densities[i].append(sim_exon_density[i]) for i in sim_exon_density]
        [sim_intron_densities[i].append(sim_intron_density[i]) for i in sim_intron_density]

    nd_exon = {}
    nd_intron = {}
    for id in exons_density:
        exon_nd = np.divide(exons_density[id] - np.mean(sim_exon_densities[id]), np.mean(sim_exon_densities[id]))
        intron_nd = np.divide(introns_density[id] - np.mean(sim_intron_densities[id]), np.mean(sim_intron_densities[id]))
        nd_exon[id] = exon_nd
        nd_intron[id] = intron_nd

    nd_exon = sequo.group_family_results(nd_exon, families)
    nd_intron = sequo.group_family_results(nd_intron, families)

    output_file = "clean_run/lincrna/tests/stop_densities_nd.csv"
    with open(output_file, "w") as outfile:
        outfile.write("id,exon,intron\n")
        for i in nd_exon:
            outfile.write("{0},{1},{2}\n".format(id, np.median(nd_exon[i]), np.median(nd_intron[i])))
