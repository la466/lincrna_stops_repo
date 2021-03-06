import generic as gen
import sim_ops as simo
import sim_ops_containers as simoc
import seq_ops as seqo
import sequence_ops as sequo
import os
import collections
from useful_motif_sets import stops, codon_map, dinucleotides
import numpy as np
from progressbar import ProgressBar
import random
import pandas as pd
import scipy.stats
import scipy
import shutil

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


def process_length_sim(input_file, output_file, families_file = None):

    data = pd.read_csv(input_file)
    ids = list(data)[1:]
    nds = {}
    z_scores = {}
    gcs = {}
    empirical_ps = {}
    reals = {}
    means = {}
    medians = {}

    for id in ids:
        gc = data[id][0]
        gcs[id] = gc
        real = data[id][1]
        reals[id] = real
        sims = data[id][2:]
        nd = np.divide(real - np.nanmean(sims), np.nanmean(sims))
        nds[id] = nd
        empirical_ps[id] = np.divide(len([i for i in sims if i >= real]) + 1, len(sims) + 1)
        z = np.divide(real - np.nanmean(sims), np.nanstd(sims))
        z_scores[id] = z
        means[id] = np.nanmean(sims)
        medians[id] = np.nanmedian(sims)

    with open(output_file, "w") as outfile:
        outfile.write("id,gc,real,mean_sims,median_sims,nd,empirical_p,z,p\n")

        if families_file:
            families = gen.read_many_fields(families_file, "\t")
            reals = sequo.group_family_results(reals, families)
            nds, nd_groups = sequo.group_family_results(nds, families, return_groups = True)
            z_scores = sequo.group_family_results(z_scores, families)
            gcs = sequo.group_family_results(gcs, families)
            means = sequo.group_family_results(means, families)
            medians = sequo.group_family_results(medians, families)

            for id in z_scores:
                real = np.nanmedian(reals[id])
                gc = np.nanmedian(gcs[id])
                z = np.nanmedian(z_scores[id])
                nd = np.nanmedian(nds[id])
                mean = np.nanmedian(means[id])
                median = np.nanmedian(medians[id])

                if not np.isnan(nd):
                    if nd in nds[id]:
                        index = nds[id].index(nd)
                    else:
                        mid_percentile = np.nanpercentile(nds[id],50,interpolation='nearest')
                        index = nds[id].index(mid_percentile)
                    loc_id = nd_groups[id][index]
                    empirical_p = empirical_ps[loc_id]
                else:
                    empirical_p = "nan"

                p = scipy.stats.norm.sf(abs(z))*2
                outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(id, gc, real, mean, median, nd, empirical_p, z, p))
        else:
            for id in z_scores:
                real = reals[id]
                gc = gcs[id]
                z = z_scores[id]
                nd = nds[id]
                mean = means[id]
                median = medians[id]
                empirical_p = empirical_ps[id]
                p = scipy.stats.norm.sf(abs(z))*2
                outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(id, gc, real, mean, median, nd, empirical_p, z, p))


def calculate_lengths(exons_fasta, output_file, families_file = None):

    names, seqs = gen.read_fasta(exons_fasta)
    lincRNA = {name: seqs[i] for i, name in enumerate(names)}

    if families_file:
        families = gen.read_many_fields(families_file, "\t")
        lincRNA = sequo.group_family_results(lincRNA, families)

    gcs = {id: [seqo.calc_seq_gc(i) for i in lincRNA[id]] for id in lincRNA}
    lengths = {id: [len(i) for i in lincRNA[id]] for id in lincRNA}

    with open(output_file, "w") as outfile:
        outfile.write("id,gc,length\n")
        for id in lengths:
            outfile.write("{0},{1},{2}\n".format(id, np.median(gcs[id]), np.median(lengths[id])))


def sim_stop_density(input_fasta, output_file, simulations = None, families_file = None, introns = None, input_fasta2 = None):
    """
    Wrapper to calculate and simulate the stop codon density in lincRNA sequences.

    Args:
        input_fasta (str): path to input fasta containing sequences
        output_file (str): path to output file
        simulations (int): if set, the number of simulations to run
        families_file (str): if set, the path to the file containing the paralogous families
        introns (bool): if true, look at introns (also need input_fasta2)
        input_fasta2 (str): if set, the path to the exons file whe looking at introns
    """
    # get the sequences
    names, sequences = gen.read_fasta(input_fasta)

    if introns:
        exon_names, exon_seqs = gen.read_fasta(input_fasta2)
        exon_list = collections.defaultdict(lambda: [])
        [exon_list[name.split(".")[0]].append(exon_seqs[i]) for i, name in enumerate(exon_names)]
        # exon_list = sequo.pick_random_family_member(families_file, exon_list)

        sequence_list = collections.defaultdict(lambda: [])
        [sequence_list[name.split(".")[0]].append(sequences[i]) for i, name in enumerate(names) if name.split(".")[0] in exon_list]
        sequence_list = {i: "".join(sequence_list[i]) for i in sequence_list}
    else:
        sequence_list = {name.split(".")[0]: sequences[i] for i, name in enumerate(names)}

    if families_file:
        sequence_list = sequo.pick_random_family_member(families_file, sequence_list)



    # create a temporary output directory
    temp_dir = "temp_lincrna_sim"
    gen.create_output_directories(temp_dir)
    # set up the simulation
    simulation_list = ["real"]
    simulation_list.extend(list(range(simulations)))
    # get a list of the output files that might have already been created
    output_filelist = {}
    for file in os.listdir(temp_dir):
        if "real" not in file:
            output_filelist[int(file.split(".")[0].split("_")[-1])] = "{0}/{1}".format(temp_dir, file)
        else:
            output_filelist[file.split(".")[0].split("_")[-1]] = "{0}/{1}".format(temp_dir, file)

    args = [sequence_list, temp_dir, output_filelist]
    # run the simulation
    outputs = simoc.run_simulation_function(simulation_list, args, simo.simulate_lincrna_stop_codon_density, sim_run = False)
    # join all outputs
    outputs = {**output_filelist, **outputs}
    real_output = outputs["real"]
    sim_outputs = {i: outputs[i] for i in outputs if i != "real"}

    with open(output_file, "w") as outfile:
        outfile.write("id,seq_count,gc,stop_codon_density\n")
        # real data
        outfile.write("{0}\n".format(gen.read_many_fields(real_output, "\t")[0][0]))
        # simulation data
        for sim_id in sorted(sim_outputs):
            data = gen.read_many_fields(sim_outputs[sim_id], "\t")[0][0]
            outfile.write("{0},{1}\n".format(sim_id + 1, ",".join(data.split(",")[1:])))
    # remove the temp directory
    gen.remove_directory(temp_dir)

def sim_stop_density_removed_motifs(input_fasta, output_file, motif_file, simulations = None, families_file = None):
    """
    Wrapper to calculate and simulate the stop codon density in lincRNA sequences.

    Args:
        input_fasta (str): path to input fasta containing sequences
        output_file (str): path to output file
        motif_file (str): path to file containing motifs
        simulations (int): if set, the number of simulations to run
        families_file (str): if set, the path to the file containing the paralogous families
    """
    # get the motifs
    motifs = sequo.read_motifs(motif_file)
    # get the sequences
    names, sequences = gen.read_fasta(input_fasta)
    sequence_list = {name.split(".")[0]: sequences[i] for i, name in enumerate(names)}

    if families_file:
        sequence_list = sequo.pick_random_family_member(families_file, sequence_list)

    # create a temporary output directory
    temp_dir = "temp_lincrna_sim_remove"
    gen.remove_directory(temp_dir)
    gen.create_output_directories(temp_dir)
    # set up the simulation
    simulation_list = ["real"]
    simulation_list.extend(list(range(simulations)))
    # get a list of the output files that might have already been created
    output_filelist = {}
    for file in os.listdir(temp_dir):
        if "real" not in file:
            output_filelist[int(file.split(".")[0].split("_")[-1])] = "{0}/{1}".format(temp_dir, file)
        else:
            output_filelist[file.split(".")[0].split("_")[-1]] = "{0}/{1}".format(temp_dir, file)

    args = [sequence_list, temp_dir, output_filelist, motifs]
    # run the simulation
    outputs = simoc.run_simulation_function(simulation_list, args, simo.simulate_lincrna_stop_codon_density_removed_motifs, sim_run = False)
    # join all outputs
    outputs = {**output_filelist, **outputs}
    real_output = outputs["real"]
    sim_outputs = {i: outputs[i] for i in outputs if i != "real"}

    with open(output_file, "w") as outfile:
        outfile.write("id,seq_count,gc,stop_codon_density\n")
        # real data
        outfile.write("{0}\n".format(gen.read_many_fields(real_output, "\t")[0][0]))
        # simulation data
        for sim_id in sorted(sim_outputs):
            data = gen.read_many_fields(sim_outputs[sim_id], "\t")[0][0]
            outfile.write("{0},{1}\n".format(sim_id + 1, ",".join(data.split(",")[1:])))
    # remove the temp directory
    gen.remove_directory(temp_dir)

def sim_stop_density_removed_motifs_seq_sim(input_fasta, output_file, motif_file, sim_dir, simulations = None, families_file = None):
    """
    Wrapper to calculate and simulate the stop codon density in lincRNA sequences.

    Args:
        input_fasta (str): path to input fasta containing sequences
        output_file (str): path to output file
        motif_file (str): path to file containing motifs
        sim_dir (str): path to simulation directory
        simulations (int): if set, the number of simulations to run
        families_file (str): if set, the path to the file containing the paralogous families
    """
    # get the motifs
    motifs = sequo.read_motifs(motif_file)
    # get the sequences
    names, sequences = gen.read_fasta(input_fasta)
    sequence_list = {name.split(".")[0]: sequences[i] for i, name in enumerate(names)}

    if families_file:
        sequence_list = sequo.pick_random_family_member(families_file, sequence_list)

    # create a temporary output directory
    temp_dir = "temp_lincrna_sim_remove"
    gen.remove_directory(temp_dir)
    gen.create_output_directories(temp_dir)
    # set up the simulation
    motif_files = {"real": motif_file}
    random_motif_files = np.random.choice(list(range(len(os.listdir(sim_dir)))), simulations)
    for i, file in enumerate(random_motif_files):
        motif_files["sim_{0}".format(i+1)] = "{0}/{1}".format(sim_dir, os.listdir(sim_dir)[file])
    simulation_list = list(motif_files)

    args = [sequence_list, temp_dir, motif_files]
    # run the simulation
    outputs = simoc.run_simulation_function(simulation_list, args, simo.simulate_lincrna_stop_codon_density_removed_motifs_sim_seq, sim_run = False)
    # join all outputs
    # outputs = {**output_filelist, **outputs}
    real_output = outputs["real"]
    sim_outputs = {i: outputs[i] for i in outputs if i != "real"}

    with open(output_file, "w") as outfile:
        outfile.write("id,seq_count,gc,stop_codon_density\n")
        # real data
        outfile.write("{0}\n".format(gen.read_many_fields(real_output, "\t")[0][0]))
        # simulation data
        for sim_id in sorted(sim_outputs):
            data = gen.read_many_fields(sim_outputs[sim_id], "\t")[0][0]
            outfile.write("{0}\n".format(",".join(data.split(","))))
    # remove the temp directory
    gen.remove_directory(temp_dir)


def process_sim_stop_density_outputs(output_dir, output_file, test_col = None, reverse = None):
    """
    Wrapper to calculate p value for lincRNA stop codon density after shuffling
    of sequences
    """

    if not test_col:
        test_col = "stop_codon_density"

    with open(output_file, "w") as outfile:
        outfile.write("run_number,seq_count,gc,density,median_simulated_density,normalised_density,p_value\n")
        files = os.scandir(output_dir)
        for i, file in sorted(enumerate(files)):
            data = pd.read_csv(file.path)
            real = data[data['id'] == "real"]
            sims = data[data['id'] != "real"]
            if reverse:
                count_list = sims[test_col] >= real[test_col][0]
            else:
                count_list = sims[test_col] <= real[test_col][0]
            p = np.divide(count_list.sum() + 1, len(sims) + 1)
            median_sims = sims[test_col].median()
            nd = np.divide(real[test_col][0] - sims[test_col].mean(), sims[test_col].mean())
            outfile.write("{0},{1},{2},{3},{4},{5},{6}\n".format(i+1, real["seq_count"][0], real["gc"][0], real[test_col][0], median_sims, nd, p))


def sim_stop_density_within_genes(input_fasta, output_file, simulations = None, families_file = None):
    """
    Wrapper to calculate and simulate the stop codon density in lincRNA sequences.

    Args:
        input_fasta (str): path to input fasta containing sequences
        output_file (str): path to output file
        simulations (int): if set, the number of simulations to run
        families_file (str): if set, the path to the file containing the paralogous families
    """
    # get the sequences
    names, sequences = gen.read_fasta(input_fasta)
    sequence_list = {name.split(".")[0]: sequences[i] for i, name in enumerate(names)}

    if families_file:
        sequence_list = sequo.pick_random_family_member(families_file, sequence_list)

    # create a temporary output directory
    temp_dir = "temp_lincrna_sim_within_genes"
    gen.create_output_directories(temp_dir)
    # get a list of the output files that might have already been created
    output_filelist = {}
    for file in os.listdir(temp_dir):
        output_filelist[file.split(".")[0]] = "{0}/{1}".format(temp_dir, file)
    simulation_list = [i for i in sequence_list if i not in output_filelist]

    args = [sequence_list, temp_dir, output_filelist, simulations]
    # run the simulation
    outputs = simoc.run_simulation_function(simulation_list, args, simo.simulate_lincrna_stop_codon_density_within_genes, sim_run = False)
    # join all outputs
    outputs = {**output_filelist, **outputs}

    with open(output_file, "w") as outfile:
        outfile.write("id,stop_density,median_simulated_stop_density,normalised_density,p_value\n")
        for id in sorted(outputs):
            data = [float(i) for i in gen.read_many_fields(outputs[id], ",")[0]]
            real = data[0]
            sims = data[1:]
            p = np.divide(len([i for i in sims if i <= real]) + 1, len(sims) + 1)
            nd = np.divide(real - np.mean(sims), np.mean(sims))
            outfile.write("{0},{1},{2},{3},{4}\n".format(id,real,np.median(sims), nd, p))
    # remove the temp directory
    gen.remove_directory(temp_dir)

def process_sim_stop_density_within_gene_outputs(output_dir, output_file):
    """
    Wrapper to calculate binomial p value
    """

    with open(output_file, "w") as outfile:
        outfile.write("run_number,data_points,depletions,depletions_binomial_p,significant_depletions,significant_depletions_binomial_p\n")
        files = os.scandir(output_dir)
        for i, file in sorted(enumerate(files)):
            data = pd.read_csv(file.path)
            depletions = data["normalised_density"] < 0
            significant = data["p_value"] < 0.05
            rows = data.shape[0]
            print(file, rows)
            binom_test_depletions = scipy.stats.binom_test(depletions.sum(), rows, p = 0.05, alternative = "greater")
            binom_test_significant_depletions = scipy.stats.binom_test(significant.sum(), rows, p = 0.05, alternative = "greater")
            outfile.write("{0},{1},{2},{3},{4},{5}\n".format(i+1, rows, depletions.sum(), binom_test_depletions, significant.sum(), binom_test_significant_depletions))


def calculcate_motif_nd(input_fasta, motif_file, output_file, simulations = None, families_file = None):

    # get the sequences
    names, sequences = gen.read_fasta(input_fasta)
    sequence_list = {name: sequences[i] for i, name in enumerate(names)}

    if families_file:
        sequence_list = sequo.pick_random_family_member(families_file, sequence_list)

    # create a temporary output directory
    temp_dir = "temp_lincrna_motif_nd"
    gen.create_output_directories(temp_dir)
    # get a list of the output files that might have already been created
    # set up the simulation
    simulation_list = ["real"]
    simulation_list.extend(list(range(simulations)))
    # get a list of the output files that might have already been created
    output_filelist = {}
    for file in os.listdir(temp_dir):
        if "real" not in file:
            output_filelist[int(file.split(".")[0].split("_")[-1])] = "{0}/{1}".format(temp_dir, file)
        else:
            output_filelist[file.split(".")[0].split("_")[-1]] = "{0}/{1}".format(temp_dir, file)

    args = [sequence_list, motif_file, temp_dir, output_filelist]
    # run the simulation
    outputs = simoc.run_simulation_function(simulation_list, args, simo.lincrna_motif_nd, sim_run = False)
    # join all outputs
    outputs = {**output_filelist, **outputs}
    real_output = outputs["real"]
    sim_outputs = {i: outputs[i] for i in outputs if i != "real"}

    real_densities = {i[0]: float(i[1]) for i in gen.read_many_fields(real_output, ",")}
    sim_densities = collections.defaultdict(lambda: [])
    for sim_id in sorted(sim_outputs):
        data = gen.read_many_fields(sim_outputs[sim_id], ",")
        [sim_densities[i[0]].append(float(i[1])) for i in data]

    nds = {motif: np.divide(real_densities[motif] - np.mean(sim_densities[motif]), np.std(sim_densities[motif])) for motif in real_densities}

    with open(output_file, "w") as outfile:
        outfile.write("motif,density,nd\n")
        [outfile.write("{0},{1},{2}\n".format(motif, real_densities[motif], nds[motif])) for motif in sorted(nds)]

    # remove the temp directory
    gen.remove_directory(temp_dir)



def excess_test(input_fasta, motif_file, output_file, simulations = None, families_file = None):
    """
    """

    if not simulations:
        simulations = 1000


    # read in the motifs
    motifs = sequo.read_motifs(motif_file)
    # read in the sequences
    names, seqs = gen.read_fasta(input_fasta)
    seq_list = collections.defaultdict(lambda: [])
    [seq_list[name.split(".")[0]].append(seqs[i]) for i, name in enumerate(names)]

    if families_file:
        seq_list = sequo.pick_random_family_member(families_file, seq_list)

    seqs = []
    [seqs.extend(seq_list[i]) for i in seq_list]

    # create output directory
    temp_dir = "temp_sim_seqs"
    gen.create_output_directories(temp_dir)
    # write the set of real sequences to a temp file
    with open("{0}/real.txt".format(temp_dir), "w") as temp_output:
        temp_output.write("{0}\n".format(",".join(seqs)))
    # simulate the sequences
    simoc.run_simulation_function(list(range(simulations)), [seqs, temp_dir], simo.shuffle_sequences, sim_run = False)
    # get the temp files
    filelist = ["{0}/{1}".format(temp_dir, i) for i in os.listdir(temp_dir)]
    # calculate the excesses
    outputs = simoc.run_simulation_function(filelist, [motifs], sequo.motif_excess_test, sim_run = False)
    # write to file
    with open(output_file, "w") as outfile:
        outfile.write("id,core_density,flank_density,excess\n")
        [outfile.write("{0}\n".format(",".join(gen.stringify(output)))) for output in outputs]
    # delete the temp directory
    gen.remove_directory(temp_dir)


    # now calculate the densities

def calc_substitution_rates(input_fasta, motif_file, required_simulations, output_file, families_file = None):
    """
    Given a set of motifs and alignments, calculate the substitution rates between
    for nucleotides that are part of stop codons and those that aren't.

    Args:
        input_fasta (str): path to input fasta
        motif_file (str): path to motif file
        required_simulations (int): number of simulations
        output_file (str): path to output file
        families_file (str): if set, path to families file
    """
    if not required_simulations:
        required_simulations = 1000
    # read in motifs
    motifs = sequo.read_motifs(motif_file)
    # read in alignments
    names, seqs = gen.read_fasta(input_fasta)
    seq_list = collections.defaultdict(lambda: [])
    [seq_list[name.split(".")[0]].append(seqs[i].split(",")) for i, name in enumerate(names)]
    # pick a random family member
    seq_list = sequo.pick_random_family_member(families_file, seq_list)

    temp_dir = "temp_sub_rates"
    gen.create_output_directories(temp_dir)

    # set the arguments
    sim_randomisations = {"real": False}
    for i in range(required_simulations):
        sim_randomisations[i+1] = True
    args = [motifs, seq_list, sim_randomisations, temp_dir]

    outputs = simoc.run_simulation_function(list(sim_randomisations), args, simoc.calculate_substitution_rates, sim_run = False)
    real_file = [i for i in outputs if "real" in i][0]
    sim_files = [i for i in outputs if i != real_file]

    # write to file
    with open(output_file, "w") as outfile:
        outfile.write("id,ese_rate,non_ese_rate,,ese_stop_rate,ese_non_stop_rate,,non_ese_stop_rate,non_ese_non_stop_rate,,relative_ese_rate,relative_non_ese_rate,relative_difference\n")
        outfile.write("{0}\n".format(gen.read_many_fields(real_file, "\t")[0][0]))
        [outfile.write("{0}\n".format(gen.read_many_fields(sim_file, "\t")[0][0])) for sim_file in sim_files]

    # remove the temp dir
    gen.remove_directory(temp_dir)

def calc_substitution_rates_motif(input_fasta, motif_file, required_simulations, controls_dir, output_file, families_file = None):
    """
    Given a set of motifs and alignments, calculate the substitution rates between
    for nucleotides that are part of stop codons and those that aren't.

    Args:
        input_fasta (str): path to input fasta
        motif_file (str): path to motif file
        required_simulations (int): number of simulations
        output_file (str): path to output file
        families_file (str): if set, path to families file
    """

    # read in motifs
    motifs = sequo.read_motifs(motif_file)
    # read in alignments
    names, seqs = gen.read_fasta(input_fasta)
    seq_list = collections.defaultdict(lambda: [])
    [seq_list[name.split(".")[0]].append(seqs[i].split(",")) for i, name in enumerate(names)]
    # pick a random family member
    seq_list = sequo.pick_random_family_member(families_file, seq_list)
    seq_list = {i: seq_list[i] for i in seq_list}

    temp_dir = "temp_sub_rates_motifs"
    gen.create_output_directories(temp_dir)

    # set the arguments
    # sim_randomisations = {"real": False}
    motif_files = {"real": motif_file}
    if required_simulations:
        for i, file in enumerate(np.random.choice(os.listdir(controls_dir), required_simulations, replace = False)):
            motif_files[i+1] = "{0}/{1}".format(controls_dir, file)
    args = [motif_file, motif_files, seq_list, temp_dir]

    outputs = simoc.run_simulation_function(list(motif_files), args, simoc.calculate_substitution_rates_motifs, sim_run = False, parallel = True)
    real_file = [i for i in outputs if "real" in i][0]
    sim_files = [i for i in outputs if i != real_file]

    # write to file
    with open(output_file, "w") as outfile:
        outfile.write("id,ese_rate,non_ese_rate,ese_stop_rate,ese_non_stop_rate,non_ese_stop_rate,non_ese_non_stop_rate,relative_ese_diff,relative_non_ese_diff,log_diff_ratio\n")
        outfile.write("{0}\n".format(gen.read_many_fields(real_file, "\t")[0][0]))
        [outfile.write("{0}\n".format(gen.read_many_fields(sim_file, "\t")[0][0])) for sim_file in sim_files]

    # remove the temp dir
    gen.remove_directory(temp_dir)


def calc_dinucleotide_substitution_rates(input_fasta, motif_file, required_simulations, output_file, families_file = None):
    """
    Given a set of motifs and alignments, calculate the substitution rates for
    each dinucleotide. Then also compare the expected number of those that
    contribute to stop codons with those that dont

    Args:
        input_fasta (str): path to input fasta
        motif_file (str): path to motif file
        required_simulations (int): number of simulations
        output_file (str): path to output file
        families_file (str): if set, path to families file
    """

    # read in motifs
    motifs = sequo.read_motifs(motif_file)
    # read in alignments
    names, seqs = gen.read_fasta(input_fasta)
    seq_list = collections.defaultdict(lambda: [])
    [seq_list[name.split(".")[0]].append(seqs[i].split(",")) for i, name in enumerate(names)]
    # pick a random family member
    seq_list = sequo.pick_random_family_member(families_file, seq_list)

    args = [seq_list, motifs]
    outputs = simoc.run_simulation_function(list(seq_list), args, sequo.get_dinucleotide_substitutions, sim_run = False)

    motif_dint_counts = collections.defaultdict(lambda: 0)
    motif_dint_subs = collections.defaultdict(lambda: 0)
    non_motif_dint_counts = collections.defaultdict(lambda: 0)
    non_motif_dint_subs = collections.defaultdict(lambda: 0)

    # gather the results
    for i, output in enumerate(outputs):
        if i+4 % 4 == 0:
            for dint in output:
                motif_dint_counts[dint] += output[dint]
        if i+4 % 4 == 1:
            for dint in output:
                motif_dint_subs[dint] += output[dint]
        if i+4 % 4 == 2:
            for dint in output:
                non_motif_dint_counts[dint] += output[dint]
        if i+4 % 4 == 3:
            for dint in output:
                non_motif_dint_subs[dint] += output[dint]


    # get the substitution rates in motifs and not in motifs
    motif_sub_rates = sequo.calc_dinucleotide_substitution_rates(motif_dint_subs, motif_dint_counts)
    non_motif_sub_rates = sequo.calc_dinucleotide_substitution_rates(non_motif_dint_subs, non_motif_dint_counts)
    # define the dinucleotides in stops and not in stops
    stop_dints = ["TA", "TG", "GA", "AA", "AG"]
    non_stop_dints = [i for i in dinucleotides if i not in stop_dints]

    # do a chisquare test to determine whether the stop dints sub are underrepresented
    motif_chisquare_outputs = sequo.calc_dint_chisquare(motif_dint_subs, motif_dint_counts, stop_dints, non_stop_dints)
    non_motif_chisquare_outputs = sequo.calc_dint_chisquare(non_motif_dint_subs, non_motif_dint_counts, stop_dints, non_stop_dints)

    # get the mean rates and percentage difference
    motif_stop_average_rate = np.mean([motif_sub_rates[i] for i in stop_dints])
    motif_non_stop_average_rate = np.mean([motif_sub_rates[i] for i in non_stop_dints])
    non_motif_stop_average_rate = np.mean([non_motif_sub_rates[i] for i in stop_dints])
    non_motif_non_stop_average_rate = np.mean([non_motif_sub_rates[i] for i in non_stop_dints])
    # diffs
    motif_diff = np.divide(motif_stop_average_rate - motif_non_stop_average_rate, motif_non_stop_average_rate)*100
    non_motif_diff = np.divide(non_motif_stop_average_rate - non_motif_non_stop_average_rate, non_motif_non_stop_average_rate)*100

    with open(output_file, "w") as outfile:
        outfile.write("within_motifs\n")
        outfile.write("\n")
        outfile.write(",counts,subs,expected_subs\n")
        outfile.write("stop_dints,{0},{1},{2}\n".format(motif_chisquare_outputs[0], motif_chisquare_outputs[1], motif_chisquare_outputs[2]))
        outfile.write("non_stop_dints,{0},{1},{2}\n".format(motif_chisquare_outputs[3], motif_chisquare_outputs[4], motif_chisquare_outputs[5]))
        outfile.write("totals,{0},{1},{2}".format(motif_chisquare_outputs[6], motif_chisquare_outputs[7], motif_chisquare_outputs[6] + motif_chisquare_outputs[7]))
        outfile.write("\n")
        outfile.write("\nchiquare_test\n")
        outfile.write("statistic:,{0}\n".format(motif_chisquare_outputs[-1].statistic))
        outfile.write("p_value:,{0}\n".format(motif_chisquare_outputs[-1].pvalue))
        outfile.write("\n\noutside_motifs\n")
        outfile.write("\n")
        outfile.write(",counts,subs,expected_subs\n")
        outfile.write("stop_dints,{0},{1},{2}\n".format(non_motif_chisquare_outputs[0], non_motif_chisquare_outputs[1], non_motif_chisquare_outputs[2]))
        outfile.write("non_stop_dints,{0},{1},{2}\n".format(non_motif_chisquare_outputs[3], non_motif_chisquare_outputs[4], non_motif_chisquare_outputs[5]))
        outfile.write("totals,{0},{1},{2}".format(non_motif_chisquare_outputs[6], non_motif_chisquare_outputs[7], non_motif_chisquare_outputs[6] + non_motif_chisquare_outputs[7]))
        outfile.write("\n")
        outfile.write("\nchiquare_test\n")
        outfile.write("statistic:,{0}\n".format(non_motif_chisquare_outputs[-1].statistic))
        outfile.write("p_value:,{0}\n".format(non_motif_chisquare_outputs[-1].pvalue))

        outfile.write("\n\naverage_rates\n")
        outfile.write(",stop_dints,non_stop_dints,diff\n")
        outfile.write("within_motifs,{0},{1},{2}\n".format(motif_stop_average_rate, motif_non_stop_average_rate, motif_diff))
        outfile.write("outside_motifs,{0},{1},{2}\n".format(non_motif_stop_average_rate, non_motif_non_stop_average_rate, non_motif_diff))

        outfile.write("\n\nindividual_rates\n")
        outfile.write("dint,within_motifs,outside_motifs\n")
        for dint in sorted(motif_sub_rates):
            outfile.write("{0},{1},{2}\n".format(dint, motif_sub_rates[dint], non_motif_sub_rates[dint]))


def calc_gc(input_fasta, output_file, families_file = None):

    # get sequences
    names, seqs = gen.read_fasta(input_fasta)
    seq_list = {name.split(".")[0]: seqs[i] for i, name in enumerate(names)}

    with open(output_file, "w") as outfile:
        outfile.write("run,count,median_gc\n")
        if families_file:
            for i in range(10):
                chosen_seqs = sequo.pick_random_family_member(families_file, seq_list)
                gcs = [seqo.calc_seq_gc(chosen_seqs[id]) for id in chosen_seqs]
                median_gc = np.median(gcs)
                outfile.write("{0},{1},{2}\n".format(i+1, len(chosen_seqs), median_gc))
        else:
            gcs = [seqo.calc_seq_gc(seq_list[id]) for id in seq_list]
            median_gc = np.median(gcs)
            outfile.write("all_seqs,{0},{1}\n".format(len(seq_list), median_gc))

def motif_overlap_test(input_fasta, motif_file, output_file, runs = 1, families_file = None):
    """
    Wrapper for the test to ask whether stop containing motifs overlap than
    non-stop containing motifs.

    Args:
        input_fasta (str): path to input fasta file
        motif_file (str): path to motif file
        output_file (str): path to output file
        runs (int): the number of times to run the test
        families_file (str): if set, path to the families_file
    """

    exon_names, exon_seqs = gen.read_fasta(input_fasta)
    exon_list = collections.defaultdict(lambda: [])
    for i, name in enumerate(exon_names):
        exon_list[name.split(".")[0]].append(exon_seqs[i])
    exon_list = {i: exon_list[i] for i in exon_list}

    motifs = sequo.read_motifs(motif_file)

    outputs = soc.run_simulation_function(list(range(runs)), [motifs, exon_list, families_file], sequo.calc_motif_overlaps, sim_run = False)
    pvals = [i[-1].pvalue for i in outputs]
    adj_pvals = [i*runs for i in pvals]

    with open(output_file, "w") as outfile:
        outfile.write("id,expected_stop_overlaps,observed_stop_overlaps,expected_non_stop_overlaps,observed_non_stop_overlaps,mean_stop_overlap_length,mean_non_stop_overlap_length,,chi_statistic,p_value,adj_p\n")
        for i, output in enumerate(outputs):
            outfile.write("run_{0},{1},,{2},{3},{4}\n".format(i+1, ",".join(gen.stringify(output[:-1])), output[-1].statistic, output[-1].pvalue, adj_pvals[i]))


def motif_overlap_density_test(input_fasta, motif_file, output_file, runs = 1, families_file = None):
    """
    Wrapper for the test to ask whether overlapping motif hits have a higher density
    of stop codons than those that do not overlap.

    Args:
        input_fasta (str): path to input fasta file
        motif_file (str): path to motif file
        output_file (str): path to output file
        runs (int): the number of times to run the test
        families_file (str): if set, path to the families_file
    """

    exon_names, exon_seqs = gen.read_fasta(input_fasta)
    exon_list = collections.defaultdict(lambda: [])
    for i, name in enumerate(exon_names):
        exon_list[name.split(".")[0]].append(exon_seqs[i])
    exon_list = {i: exon_list[i] for i in exon_list}

    motifs = sequo.read_motifs(motif_file)

    real_output = soc.run_simulation_function(list(range(runs)), [motifs, exon_list, families_file], sequo.calc_motif_overlap_density, sim_run = False)[0]
    outputs = soc.run_simulation_function(list(range(runs)), [motifs, exon_list, families_file], sequo.calc_motif_overlap_density, kwargs_dict = {"sim": True}, sim_run = False)

    with open(output_file, "w") as outfile:
        outfile.write("sim_id,non_overlap_stop_density,overlap_stop_density\n")
        outfile.write("real,{0},{1}\n".format(real_output[0], real_output[1]))
        for i, output in enumerate(outputs):
            outfile.write("sim_{0},{1},{2}\n".format(i+1,output[0], output[1]))



def clean_alignments(input_bed, alignments_file, output_exon_file, output_intron_file):
    """
    Extract the exons and introns from a full alignments file. Note this uses alignements
    extracted from https://usegalaxy.org/ which have spaces below each fasta entry.

    Args:
        input_bed (str): path to bed file
        alignments_file (str): path to alignments file
        output_exon_file (str): path to file containing exon alignment outputs
        output_introns_file (str): path to file containing intron alignment outputs
    """

    # get the sequence coordinates to extract the parts
    bed_entries = gen.read_many_fields(input_bed, "\t")
    entries = collections.defaultdict(lambda: collections.defaultdict())
    for entry in bed_entries:
        entries[int(entry[1])][int(entry[2])] = entry

    # now extract the parts
    with open(output_exon_file, "w") as exon_outfile:
        with open(output_intron_file, "w") as intron_outfile:
            with open(alignments_file, "r") as input_file:
                lines = input_file.readlines()
                lines = [i.strip("\n") for i in lines]
                for entry_id in range(0, len(lines), 5):
                    if entry_id:
                        human_name = lines[entry_id]
                        human_seq = lines[entry_id+1]
                        mac_name = lines[entry_id+2]
                        mac_seq = lines[entry_id+3]

                        if "N" not in human_seq and "N" not in mac_seq:
                            # ensure that there isnt too many indels
                            if np.divide(len([i for i in human_seq if i == "-"]), len(human_seq)) < 0.15 and np.divide(len([i for i in mac_seq if i == "-"]), len(mac_seq)) < 0.15:
                                coordinates = [int(i) for i in human_name.split(":")[-1].split("-")]
                                bed_entry = entries[coordinates[0]][coordinates[1]]
                                id = bed_entry[3]
                                exon_sizes = [int(i) for i in bed_entry[-2].split(",") if len(i)]
                                exon_starts = [int(i) for i in bed_entry[-1].split(",") if len(i)]
                                exon_coordinates = [[start, start + exon_sizes[i]] for i, start in enumerate(exon_starts)]
                                intron_coordinates = [[start + exon_sizes[i], exon_starts[i+1]] for i, start in enumerate(exon_starts) if i+1 < len(exon_starts)]

                                human_exon_sequences = [human_seq[i[0]:i[1]].upper() for i in exon_coordinates]
                                mac_exon_sequences = [mac_seq[i[0]:i[1]].upper() for i in exon_coordinates]
                                [exon_outfile.write(">{0}.{1}\n{2},{3}\n".format(id, i+1, human_exon_sequences[i], mac_exon_sequences[i])) for i in range(len(human_exon_sequences))]

                                human_intron_sequences = [human_seq[i[0]:i[1]].upper() for i in intron_coordinates]
                                mac_intron_sequences = [mac_seq[i[0]:i[1]].upper() for i in intron_coordinates]
                                [intron_outfile.write(">{0}.{1}\n{2},{3}\n".format(id, i+1, human_intron_sequences[i], mac_intron_sequences[i])) for i in range(len(human_intron_sequences))]


def sim_stop_density_diff(input_fasta, output_file, motif_file, sim_dir, simulations = None, families_file = None):
    """
    Wrapper to calculate the differnece in stop density between ESE hits and the
    remaining bits of sequence.

    Args:
        input_fasta (str): path to input fasta containing sequences
        output_file (str): path to output file
        simulations (int): if set, the number of simulations to run
        families_file (str): if set, the path to the file containing the paralogous families
    """
    # get the sequences
    names, sequences = gen.read_fasta(input_fasta)
    sequence_list = {name.split(".")[0]: sequences[i] for i, name in enumerate(names)}

    if families_file:
        sequence_list = sequo.pick_random_family_member(families_file, sequence_list)

    # create a temporary output directory
    temp_dir = "temp_lincrna_sim_diff"
    gen.remove_directory(temp_dir)
    gen.create_output_directories(temp_dir)
    # set up the simulation
    motif_files = {"real": motif_file}
    random_motif_files = np.random.choice(list(range(len(os.listdir(sim_dir)))), simulations)
    for i, file in enumerate(random_motif_files):
        motif_files["sim_{0}".format(i+1)] = "{0}/{1}".format(sim_dir, os.listdir(sim_dir)[file])
    simulation_list = list(motif_files)

    args = [sequence_list, temp_dir, motif_files]
    # run the simulation
    outputs = simoc.run_simulation_function(simulation_list, args, simo.simulate_lincrna_stop_codon_density_diff, sim_run = False)
    # join all outputs
    # outputs = {**output_filelist, **outputs}
    real_output = outputs["real"]
    sim_outputs = {i: outputs[i] for i in outputs if i != "real"}

    with open(output_file, "w") as outfile:
        outfile.write("id,seq_count,hit_stop_density,non_hit_stop_density,diff\n")
        # real data
        outfile.write("{0}\n".format(gen.read_many_fields(real_output, "\t")[0][0]))
        # simulation data
        for sim_id in sorted(sim_outputs):
            data = gen.read_many_fields(sim_outputs[sim_id], "\t")[0][0]
            outfile.write("{0}\n".format(",".join(data.split(","))))
    # remove the temp directory
    gen.remove_directory(temp_dir)

def process_sim_stop_density_diffs(output_dir, output_file, greater_than = True):
    """
    Wrapper to process stop density diff outputs
    """

    test_col = "diff"
    with open(output_file, "w") as outfile:
        outfile.write("run_number,seq_count,hit_stop_density,non_hit_stop_density,diff,median_sim_diff,normalised_diff,p_value\n")
        files = os.scandir(output_dir)
        for i, file in sorted(enumerate(files)):
            data = pd.read_csv(file.path)
            real = data[data['id'] == "real"]
            sims = data[data['id'] != "real"]
            if greater_than:
                count_list = sims[test_col] >= real[test_col][0]
            else:
                count_list = sims[test_col] <= real[test_col][0]
            p = np.divide(count_list.sum() + 1, len(sims) + 1)
            median_sims = sims[test_col].median()
            nd = np.divide(real[test_col][0] - sims[test_col].mean(), sims[test_col].mean())
            outfile.write("{0},{1},{2},{3},{4},{5},{6},{7}\n".format(i+1, real["seq_count"][0], real["hit_stop_density"][0], real["non_hit_stop_density"][0], real[test_col][0], median_sims, nd, p))


def upstream_atg(input_fasta, output_file, simulations = None, families_file = None):
    """
    Calculate the stop codon density before the first ATG.

    Args:
        input_fasta (str): path to input file
        output_file (str): path to output file
        families_file (str): if set, path to families file
    """

    # read in the sequences
    names, seqs = gen.read_fasta(input_fasta)
    all_seqs = dict(zip([name.split(".")[0] for name in names], seqs))
    # pick random family member if required
    if families_file:
        all_seqs = sequo.pick_random_family_member(families_file, all_seqs)

    length_outputs = []

    max_length = 100
    step = 10
    for length in list(range(step, max_length+step, step)):
        print(length)
        # keep only the bits of sequences before the atg
        kept = []
        for id in all_seqs:
            seq = all_seqs[id]
            try:
                atg = seq.index("ATG")
                if atg > length:
                    kept.append(seq[:atg])
            except:
                pass

        real_density = seqo.calc_motif_density(kept, stops)

        sim_densities  = []
        sims = list(range(simulations))
        sim_densities = simoc.run_simulation_function(sims, [kept], sequo.calc_upstream_sim_densities, sim_run = False)
        # calculate the number of simulants with density less than the real
        less_than = [i for i in sim_densities if i <= real_density]
        output = [length, real_density, np.median(sim_densities), np.divide(len(less_than) + 1, len(sims) + 1)]
        length_outputs.append(output)

    # write outputs to file
    with open(output_file, "w") as outfile:
        outfile.write("minimum_upstream_threshold,real_density,median_sim_density,p_value\n")
        [outfile.write("{0}\n".format(",".join(gen.stringify(output)))) for output in length_outputs]

def exon_intron_stop_density(exon_fasta, intron_fasta, output_file, families_file = None):
    """
    Get the stop density in corresponding exons and introns

    Args:
        exon_fasta (str): path to exon fasta
        intron_fasta (str): path to intron_fasta
        output_file (str): path to intron_fasta
        families_file (str): if set, path to families file
    """
    # get exons
    names, seqs = gen.read_fasta(exon_fasta)
    exon_list = collections.defaultdict(lambda: [])
    [exon_list[name.split(".")[0]].append(seqs[i]) for i, name in enumerate(names)]
    if families_file:
        exon_list = sequo.pick_random_family_member(families_file, exon_list)

    # get introns
    intron_names, intron_seqs = gen.read_fasta(intron_fasta)
    intron_list = collections.defaultdict(lambda: [])
    [intron_list[name.split(".")[0]].append(intron_seqs[i]) for i, name in enumerate(intron_names) if name.split(".")[0] in exon_list]

    with open(output_file, "w") as outfile:
        outfile.write("id,exon_density,intron_density\n")
        for id in exon_list:
            exon_density = seqo.calc_motif_density(exons[id], stops)
            intron_density = seqo.calc_motif_density(introns[id], stops)
            outfile.write("{0},{1},{2}\n".format(id, exon_density, intron_density))


def orf_exceed_length_threshold(input_fasta, output_directory, required_simulations = None, families_file = None):
    """
    Calculate the number of ORFs that exceed given thresholds compared
    with simulated sequences.

    Args:
        input_fasta (str): path to input transcript file
        output_directory (str): path to directory to hold outputs
        required_simulations (int): if set, the number of simulations to run
        families_file (str): if set, path to families file
    """

    names, seqs = gen.read_fasta(input_fasta)
    seqs = {name.split(".")[0]: seqs[i] for i, name in enumerate(names)}
    # join all seqeunces together
    seqs_string = "".join([seqs[i] for i in seqs])
    # pick only one sequence per family for tests
    seqs = sequo.pick_random_family_member(families_file, seqs)
    # get the lengths of sequences
    seqs = [seqs[i] for i in seqs]
    seqs_lengths = [len(i) for i in seqs]

    real_longest_orfs = seqo.get_longest_orfs(seqs)

    simulant_longest_orfs = simoc.run_simulation_function(list(range(required_simulations)), [seqs_string, seqs_lengths], calc_sim_orfs, sim_run = False)

    # now for each threshold, see how many exceed
    for threshold in list(range(200,610,10)):
        real_greater = len([i for i in real_longest_orfs if i >= threshold])
        output_file = "{0}/{1}.csv".format(output_dir, threshold)
        with open(output_file, "w") as outfile:
            outfile.write("id,greater,max_orf,total\n")
            outfile.write("real,{0},{1},{2}\n".format(real_greater, max(real_longest_orfs), len(real_longest_orfs)))
            for i in simulant_longest_orfs:
                sim_output = simulant_longest_orfs[i]
                outfile.write("sim_{0},{1},{2},{3}\n".format(i, len([j for j in sim_output if j >= threshold]), max(sim_output), len(sim_output)))


def calc_sim_orfs(sims, seqs_string, seqs_lengths):
    """
    For each simulation, shuffle sequence and pick a simulant set of sequences
    with the same lengths as the real ones. Then calculate the longest orfs.

    Args:
        sims (list): list of simulations to iterate over
        seqs_string (str): sequence to take simulants from
        seqs_lengths (list): list of the lengths of the real sequeneces

    Returns:
        outputs (dict): longest orfs for each simulation
    """

    outputs = {}
    if len(sims):
        np.random.seed()
        for i, sim in enumerate(sims):
            gen.print_parallel_status(i, sims)
            sim_nts = list(seqs_string)
            np.random.shuffle(sim_nts)
            sim_string = "".join(sim_nts)
            new_seqs = []
            index = 0
            for length in seqs_lengths:
                new_seqs.append(sim_string[index:index+length])
            sim_orfs = seqo.get_longest_orfs(new_seqs)
            outputs[sim] = sim_orfs
    return outputs


def randomise_seqs(seqs):
    new_seqs = []
    for seq in seqs:
        nts = list(seq)
        np.random.shuffle(nts)
        new_seqs.append("".join(nts))
    return new_seqs


def randomise_regions(iterations, seq_list):
    output = {}
    if len(iterations):
        np.random.seed()
        for i, iteration in enumerate(iterations):
            gen.print_parallel_status(i, iterations)
            output[iteration] = {region: randomise_seqs(seq_list[region]) for region in seq_list}
    return output


def density_regions(input_fasta, motif_file, output_file, output_file1, output_file2, required_simulations = None, families_file = None):
    """
    """
    # get the motifs
    motifs = sequo.read_motifs(motif_file)
    # read in the seqs only keeping those longer than 207 nts
    names, seqs = gen.read_fasta(input_fasta)
    seqs_list = collections.defaultdict(lambda: [])
    [seqs_list[name.split(".")[0]].append(seqs[i]) for i, name in enumerate(names) if len(seqs[i]) >= 207]
    if families_file:
        seqs_list = sequo.pick_random_family_member(families_file, seqs_list)

    # get the different parts of the sequences
    seq_parts = collections.defaultdict(lambda: [])
    for id in seqs_list:
        for seq in seqs_list[id]:
            seq_parts["five"].append(seq[2:69])
            seq_parts["three"].append(seq[-69:-2])
            midpoint = int(len(seq)/2)
            seq_parts["core"].append(seq[midpoint-33:midpoint+34])
    seq_parts = {i: seq_parts[i] for i in seq_parts}


    # do chisquare test on hit counts
    region_hit_counts = {region: seqo.calc_motif_counts(seq_parts[region], stops) for region in seq_parts}
    total_hits = sum([region_hit_counts[i] for i in region_hit_counts])
    sequence_total_lengths = {region: sum([len(i) for i in seq_parts[region]]) for region in seq_parts}
    total_lengths = sum([sequence_total_lengths[i] for i in sequence_total_lengths])
    expecteds = {i: np.divide(sequence_total_lengths[i], total_lengths)*total_hits for i in sorted(sequence_total_lengths)}
    observeds = {i: region_hit_counts[i] for i in sorted(sequence_total_lengths)}
    expected = [expecteds[i] for i in expecteds]
    observed = [observeds[i] for i in observeds]
    chisq = scipy.stats.chisquare(observed, f_exp = expected)
    oe = {i: np.divide(observeds[i], expecteds[i]) for i in expecteds}

    ese_densities = {region: seqo.calc_motif_density(seq_parts[region], motifs) for region in seq_parts}
    stop_densities = {region: seqo.calc_motif_density(seq_parts[region], stops) for region in seq_parts}

    ese_densities_per_seq = collections.defaultdict(lambda: [])
    stop_densities_per_seq = collections.defaultdict(lambda: [])
    for region in seq_parts:
        for seq in seq_parts[region]:
            ese_densities_per_seq[region].append(seqo.calc_motif_density([seq], motifs))
            stop_densities_per_seq[region].append(seqo.calc_motif_density([seq], stops))

    with open(output_file, "w") as outfile:
        outfile.write("region,nucleotides,prop_nucleotides,stop_hits,prop_stop_hits,expected,o/e\n")
        for region in ["five", "core", "three"]:
            output = [region, sequence_total_lengths[region], np.divide(sequence_total_lengths[region], total_lengths), region_hit_counts[region], np.divide(region_hit_counts[region], total_hits), expecteds[region], oe[region]]
            outfile.write("{0}\n".format(",".join(gen.stringify(output))))
        outfile.write("total,{0},,{1},,{2}".format(total_lengths, total_hits, sum(expected)))
        outfile.write("\n")
        outfile.write("\nChisquare\n")
        outfile.write("Statistic:,{0}\n".format(chisq.statistic))
        outfile.write("P value:,{0}\n".format(chisq.pvalue))

    sims = simoc.run_simulation_function(list(range(1, required_simulations + 1)), [seq_parts], randomise_regions, sim_run = False)
    sim_stop_densities = collections.defaultdict(lambda: [])
    sim_seqs = collections.defaultdict(lambda: collections.defaultdict(lambda: []))
    for it in sims:
        for region in sims[it]:
            stop_density = seqo.calc_motif_density(sims[it][region], stops)
            sim_stop_densities[region].append(stop_density)

            for i, seq in enumerate(sims[it][region]):
                sim_seqs[region][i].append(seq)

    nds = {}
    for region in stop_densities:
        nds[region] = np.divide(stop_densities[region] - np.mean(sim_stop_densities[region]), np.mean(sim_stop_densities[region]))

    for region in per_seq_nds:
        print(region, per_seq_nds[region])
    with open(output_file1, "w") as outfile:
        outfile.write("region,ese_density,stop_density,mean_sim_stop_density,nd\n")
        for region in ["five", "core", "three"]:
            outfile.write("{0},{1},{2},{3},{4}\n".format(region, ese_densities[region], stop_densities[region], np.mean(sim_stop_densities[region]), nds[region]))


    with open(output_file2, "w") as outfile:
        for region in ["five", "core", "three"]:
            outfile.write("{0}_ese,{1}\n".format(region, ",".join(gen.stringify(ese_densities_per_seq[region]))))
            outfile.write("{0}_stop,{1}\n".format(region, ",".join(gen.stringify(stop_densities_per_seq[region]))))
