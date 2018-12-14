import generic as gen
import seq_ops as seqo
import ops
import numpy as np
import random
import re
import collections
from useful_motif_sets import dinucleotides, nucleotides, stops
import multiprocessing as mp
from progressbar import ProgressBar

def generate_dint_controls(input_ids, input_seqs, dinucleotide_content, nucleotide_content, output_directory):

    required_simulations = 1000

    temp_output_directory = "temp_dint_sims"
    gen.create_output_directories(temp_output_directory)

    outputs = []

    for i, id in enumerate(input_ids):
        print("(W{0}) {1}/{2}: {3}".format(mp.current_process().name.split("-")[-1], i+1, len(input_ids), id))
        simulations = list(range(1))
        inputs = [input_seqs[id]] * required_simulations
        temp_file = generate_dinucleotide_matched_seqs(simulations, inputs, dinucleotide_content, nucleotide_content, temp_output_directory)[0]

        output_file = "{0}/{1}.txt".format(output_directory, id)
        outputs.append(output_file)
        seqs = gen.read_fasta(temp_file)[1]
        with open(output_file, "w") as outfile:
            outfile.write("{0}\n".format(",".join(seqs)))
        gen.remove_file(temp_file)

    return []

def get_gc_matched_seqs(sequence_ids, sequence_list, untranscribed_sequence, output_directory, required_simulations):

    temp_dir = "temp_files"
    gen.create_output_directories(temp_dir)

    for i, id in enumerate(sequence_ids):
        print("(W{0}) {1}/{2}: {3}".format(mp.current_process().name.split("-")[-1], i+1, len(sequence_ids), id))
        output_file = "{0}/{1}.txt".format(output_directory, id)
        temp_file = "{0}/{1}.txt".format(temp_dir, random.random())
        seq = sequence_list[id]
        # set the number of simulation for each seq to do
        input_seqs = {"{0}_{1}".format(id, i+1): seq for i in range(required_simulations)}
        seqo.get_gc_matched_seqs(input_seqs, untranscribed_sequence, 0.05, temp_file)
        with open(output_file, "w") as outfile:
            outfile.write("{0}".format(",".join(gen.read_fasta(temp_file)[1])))
        gen.remove_file(temp_file)

    return []


def get_seq_seed(seq_seeds, simulation_number, sequence_number):
    """
    Return a predefined seed if set for testing

    Args:
        seq_seeds (list): a list of list containing seeds
        simulation_number (int): the simulation number
        sequence_number (int): sequence i in the simulation

    Returns:
        seq_seed: the predefined seed
    """

    if seq_seeds:
        seq_seed = seq_seeds[simulation_number][sequence_number]
    else:
        seq_seed = None
    return seq_seed


def set_seed(seeds, simulation_number):
    """
    Optionally set seed if required

    Args:
        seeds (list): a list of seeds
        simulation_number (int): the simulation number
    """

    if seeds:
        seed = seeds[simulation_number]
    else:
        seed = None
    np.random.seed(seed)


def sim_exon_flank_counts(simulations, seqs, seq_frames, dinucleotide_content, nucleotide_content, window_start, window_end, temp_dir, seeds=None, seq_seeds=None):

    temp_files = []

    dinucleotide_choices = [dn for dn in sorted(dinucleotide_content)]
    dinucleotide_probabilities = [dinucleotide_content[dn] for dn in sorted(dinucleotide_content)]
    nucleotide_choices = [n for n in sorted(nucleotide_content)]
    nucleotide_probabilities = [nucleotide_content[n] for n in sorted(nucleotide_content)]

    exons_to_exclude = []

    for sim_no, simulation in enumerate(simulations):

        # set the seed
        set_seed(seeds, simulation)

        # print the simulation number out
        gen.print_simulation(sim_no+1, simulations)

        temp_file = "{0}/{1}.{2}.txt".format(temp_dir, random.random(), simulation+1)
        temp_files.append(temp_file)


        with open(temp_file, "w") as outfile:
            # Generate a list of nucleotide matched sequences
            for i, name in enumerate(seqs):
                # this will prevent simulating sequences that have already been excluded
                seqs_list = seqs[name]
                seq_count = []
                for index, seq in enumerate(seqs_list):
                    query_name = "{0}.{1}".format(name, index)
                    if query_name not in exons_to_exclude:
                        seq = seqs_list[index]
                        frame = seq_frames[name][index]
                        generated = False
                        seq_seed = get_seq_seed(seq_seeds, sim_no, i)

                        while not generated:
                            sim_seq, failed = seqo.generate_nt_matched_seq_frame(seq, dinucleotide_choices, dinucleotide_probabilities, nucleotide_choices, nucleotide_probabilities, seq_frame=frame, seed=seq_seed)
                            if failed:
                                generated = True
                                exons_to_exclude.append(query_name)
                                seq_count.append(np.nan)
                            else:
                                generated = True
                                sim_seqs = [sim_seq]
                                counts = seqo.get_stop_counts(sim_seqs)
                                seq_count.append(counts[0])
                    else:
                        seq_count.append(np.nan)

                outfile.write(">{0}\n{1},{2}\n".format(name,seq_count[0],seq_count[1]))

    return temp_files, exons_to_exclude


def sim_exon_region_counts(simulations, seqs, seq_frames, dinucleotide_content, nucleotide_content, window_start, window_end, temp_dir, seeds=None, seq_seeds=None):

    temp_files = []

    dinucleotide_choices = [dn for dn in sorted(dinucleotide_content)]
    dinucleotide_probabilities = [dinucleotide_content[dn] for dn in sorted(dinucleotide_content)]
    nucleotide_choices = [n for n in sorted(nucleotide_content)]
    nucleotide_probabilities = [nucleotide_content[n] for n in sorted(nucleotide_content)]

    exons_to_exclude = []

    for sim_no, simulation in enumerate(simulations):

        # set the seed
        set_seed(seeds, simulation)

        # print the simulation number out
        gen.print_simulation(sim_no+1, simulations)

        temp_file = "{0}/{1}.{2}.txt".format(temp_dir, random.random(), simulation+1)
        temp_files.append(temp_file)

        with open(temp_file, "w") as outfile:
            # Generate a list of nucleotide matched sequences
            for i, name in enumerate(seqs):
                # this will prevent simulating sequences that have already been excluded
                if name not in exons_to_exclude:
                    seq = seqs[name]
                    frame = seq_frames[name]
                    generated = False

                    seq_seed = get_seq_seed(seq_seeds, sim_no, i)

                    while not generated:
                        sim_seq, failed = seqo.generate_nt_matched_seq_frame(seq, dinucleotide_choices, dinucleotide_probabilities, nucleotide_choices, nucleotide_probabilities, seq_frame=frame, seed=seq_seed)
                        if failed:
                            generated = True
                            exons_to_exclude.append(name)
                        else:
                            generated = True
                            sim_seqs = {name: sim_seq}
                            counts = ops.get_region_stop_counts(sim_seqs, window_start, window_end)
                            outfile.write(">{0}\n{1},{2}\n".format(name, counts[0],counts[1]))

    return temp_files, exons_to_exclude


def sim_orf_lengths(simulations, seqs, dinucleotide_content, nucleotide_content, temp_dir, seeds=None, seq_seeds=None):
    """
    Simulation to look at the lengths of ORFs in sequences
    compared with simulated null sequences

    Args:
        simulations (list): list of simluations to iterate over
        seqs (dict): a dictionary of sequences
        dinucleotide_content (dict): dictionary containing dinucleotide proportions
        temp_dir (str): a temporary directory to hold the outputs of the simulations
        seeds (list): list of seeds to be used for the randomisations

    Returns:
        temp_files (list): list containing file paths to simulation outputs
    """

    temp_files = []

    dinucleotide_choices = [dn for dn in sorted(dinucleotide_content)]
    dinucleotide_probabilities = [dinucleotide_content[dn] for dn in sorted(dinucleotide_content)]
    nucleotide_choices = [n for n in sorted(nucleotide_content)]
    nucleotide_probabilities = [nucleotide_content[n] for n in sorted(nucleotide_content)]

    for sim_no, simulation in enumerate(simulations):

        # set the seed
        set_seed(seeds, simulation)

        # print the simulation number out
        gen.print_simulation(sim_no+1, simulations)

        simulated_seqs = {}
        # Generate a list of nucleotide matched sequences
        for i, name in enumerate(seqs):
            seq = seqs[name]
            generated = False

            seq_seed = get_seq_seed(seq_seeds, sim_no, i)

            while not generated:
                sim_seq = seqo.generate_nt_matched_seq(seq, dinucleotide_choices, dinucleotide_probabilities, nucleotide_choices, nucleotide_probabilities, seed=seq_seed)
                if sim_seq not in simulated_seqs:
                    generated = True
                    simulated_seqs[name] = sim_seq

        longest_orfs = seqo.get_longest_orfs(simulated_seqs)
        temp_file = "{0}/{1}.{2}.txt".format(temp_dir, random.random(), simulation+1)
        temp_files.append(temp_file)
        with open(temp_file, "w") as temp:
            for name in longest_orfs:
                temp.write(">{0}\n{1}\n".format(name, longest_orfs[name]))

    return temp_files


def sim_stop_counts(simulations, seqs, dinucleotide_content, nucleotide_content, temp_dir, seeds=None, seq_seeds=None):
    """
    Simulation to count the number of stop codons in any reading frame
    compared with simulated null sequences

    Args:
        simulations (list): list of simluations to iterate over
        seqs (dict): a dictionary of sequences
        dinucleotide_content (dict): dictionary containing dinucleotide proportions
        temp_dir (str): a temporary directory to hold the outputs of the simulations
        seeds (list): list of seeds to be used for the randomisations

    Returns:
        temp_files (list): list containing file paths to simulation outputs
    """

    temp_files = []

    dinucleotide_choices = [dn for dn in sorted(dinucleotide_content)]
    dinucleotide_probabilities = [dinucleotide_content[dn] for dn in sorted(dinucleotide_content)]
    nucleotide_choices = [n for n in sorted(nucleotide_content)]
    nucleotide_probabilities = [nucleotide_content[n] for n in sorted(nucleotide_content)]

    for sim_no, simulation in enumerate(simulations):

        # set the seed
        set_seed(seeds, simulation)

        # print the simulation number out
        gen.print_simulation(sim_no+1, simulations)

        simulated_seqs = {}
        # Generate a list of nucleotide matched sequences
        for i, name in enumerate(seqs):
            seq = seqs[name]
            generated = False

            seq_seed = get_seq_seed(seq_seeds, sim_no, i)

            while not generated:
                sim_seq = seqo.generate_nt_matched_seq(seq, dinucleotide_choices, dinucleotide_probabilities, nucleotide_choices, nucleotide_probabilities, seed=seq_seed)
                if sim_seq not in simulated_seqs:
                    generated = True
                    simulated_seqs[name] = sim_seq

        stop_counts = seqo.get_stop_counts(simulated_seqs)
        temp_file = "{0}/{1}.{2}.txt".format(temp_dir, random.random(), simulation+1)
        temp_files.append(temp_file)
        with open(temp_file, "w") as temp:
            for name in stop_counts:
                temp.write(">{0}\n{1}\n".format(name, stop_counts[name]))

    return temp_files


def sim_motif_stop_counts(simulations, seqs, dinucleotide_content, nucleotide_content, temp_dir, seeds=None, seq_seeds=None):
    """
    Simulation to count the number of stop codons in any reading frame
    compared with simulated null sequences

    Args:
        simulations (list): list of simluations to iterate over
        seqs (dict): a dictionary of sequences
        dinucleotide_content (dict): dictionary containing dinucleotide proportions
        temp_dir (str): a temporary directory to hold the outputs of the simulations
        seeds (list): list of seeds to be used for the randomisations

    Returns:
        temp_files (list): list containing file paths to simulation outputs
    """

    temp_files = []

    dinucleotide_choices = [dn for dn in sorted(dinucleotide_content)]
    dinucleotide_probabilities = [dinucleotide_content[dn] for dn in sorted(dinucleotide_content)]
    nucleotide_choices = [n for n in sorted(nucleotide_content)]
    nucleotide_probabilities = [nucleotide_content[n] for n in sorted(nucleotide_content)]

    for sim_no, simulation in enumerate(simulations):

        # set the seed
        set_seed(seeds, simulation)

        # print the simulation number out
        gen.print_simulation(sim_no+1, simulations)

        simulated_seqs = []
        # Generate a list of nucleotide matched sequences
        for i, seq in enumerate(seqs):
            generated = False
            seq_seed = get_seq_seed(seq_seeds, sim_no, i)

            while not generated:
                sim_seq = seqo.generate_nt_matched_seq(seq, dinucleotide_choices, dinucleotide_probabilities, nucleotide_choices, nucleotide_probabilities, seed=seq_seed)
                if sim_seq not in simulated_seqs:
                    generated = True
                    simulated_seqs.append(sim_seq)

        stop_counts = seqo.get_stop_counts(simulated_seqs)
        temp_file = "{0}/{1}.{2}.txt".format(temp_dir, random.random(), simulation+1)
        temp_files.append(temp_file)
        with open(temp_file, "w") as temp:
            temp.write(">{0}\n{1}\n".format(simulation+1, sum(stop_counts)))

    return temp_files


def sim_motif_codon_densities(simulations, seqs, dinucleotide_content, nucleotide_content, query_set, seeds=None, seq_seeds=None):

    sim_densities = []

    dinucleotide_choices = [dn for dn in sorted(dinucleotide_content)]
    dinucleotide_probabilities = [dinucleotide_content[dn] for dn in sorted(dinucleotide_content)]
    nucleotide_choices = [n for n in sorted(nucleotide_content)]
    nucleotide_probabilities = [nucleotide_content[n] for n in sorted(nucleotide_content)]

    for sim_no, simulation in enumerate(simulations):
        # set the seed
        set_seed(seeds, simulation)

        # # print the simulation number out
        # gen.print_simulation(sim_no+1, simulations)

        simulated_seqs = []
        # Generate a list of nucleotide matched sequences
        for i, seq in enumerate(seqs):
            generated = False
            seq_seed = get_seq_seed(seq_seeds, sim_no, i)
            while not generated:
                sim_seq = seqo.generate_nt_matched_seq(seq, dinucleotide_choices, dinucleotide_probabilities, nucleotide_choices, nucleotide_probabilities, seed=seq_seed)
                if sim_seq not in simulated_seqs:
                    generated = True
                    simulated_seqs.append(sim_seq)

        sim_density = seqo.calc_codon_density_in_seqs(query_set, simulated_seqs)
        sim_densities.append(sim_density)

    return sim_densities

def generate_dinucleotide_matched_seqs(simulations, seqs, dinucleotide_content, nucleotide_content, output_directory, seeds=None, seq_seeds=None, match_density = None):

    dinucleotide_choices = [dn for dn in sorted(dinucleotide_content)]
    dinucleotide_probabilities = [dinucleotide_content[dn] for dn in sorted(dinucleotide_content)]
    nucleotide_choices = [n for n in sorted(nucleotide_content)]
    nucleotide_probabilities = [nucleotide_content[n] for n in sorted(nucleotide_content)]

    output_files = []

    pbar = ProgressBar()

    if len(simulations):
        for sim_no, simulation in enumerate(pbar(simulations)):
            generated_set = False
            while not generated_set:
                # set the seed
                set_seed(seeds, simulation)

                # print the simulation number out
                # print("(W{0}) {1}/{2}: {3}".format(mp.current_process().name.split("-")[-1], sim_no+1, len(simulations), id))

                simulated_seqs = []
                # Generate a list of nucleotide matched sequences
                for i, seq in enumerate(seqs):
                    generated_seq = False
                    seq_seed = get_seq_seed(seq_seeds, sim_no, i)

                    if match_density:
                        stop_density = seqo.calc_seqs_codon_set_density([seq], codon_set = stops)

                    while not generated_seq:
                        sim_seq = seqo.generate_nt_matched_seq(seq, dinucleotide_choices, dinucleotide_probabilities, nucleotide_choices, nucleotide_probabilities, seed=seq_seed)
                        if match_density:
                            sim_density = seqo.calc_seqs_codon_set_density([sim_seq], codon_set = stops)

                        if sim_seq not in simulated_seqs:
                            if not match_density or match_density and stop_density == sim_density:
                                generated_seq = True
                                simulated_seqs.append(sim_seq)
                if sorted(seqs) != sorted(simulated_seqs):
                    generated_set = True

            output_file = "{0}/{1}.txt".format(output_directory, random.random())
            output_files.append(output_file)
            with open(output_file, "w") as outfile:
                for i, seq in enumerate(simulated_seqs):
                    outfile.write(">{0}\n{1}\n".format(i+1,seq))

    return output_files


def generate_matched_gc_seqs(simulations, seqs_fasta, non_features_seqs, threshold, output_dir, seeds=None):

    temp_files = []

    if len(simulations):
        # get the sequences
        input_names, input_seqs = gen.read_fasta(seqs_fasta)
        input_seqs_list = {name: input_seqs[i] for i, name in enumerate(input_names)}
        genome_seqs = non_features_seqs
        # get the non-transcribed region of the genome as a string to query
        genome_query_string = "".join(genome_seqs)

        for sim_no, simulation in enumerate(simulations):
            # set the seed
            set_seed(seeds, simulation)
            # print the simulation number out
            gen.print_simulation(sim_no+1, simulations)
            temp_file = "{0}/{1}.{2}.txt".format(output_dir, random.random(), simulation+1)
            seqo.get_gc_matched_seqs(input_seqs_list, genome_query_string, threshold, temp_file)
            temp_files.append(temp_file)

    return temp_files


def randomise_multiarray(multiarray, axis=1):
    """
    Shuffle a numpy array along a certain axis

    Args:
        multiarray (np.array): an array of lists, i.e. list of lists
        axis (int): the axis to randomise along (0 = down, 1 = across)

    Returns:
        randomised (np.array): the randomised array
    """

    randomised_positions = np.random.random(multiarray.shape)
    ids = np.argsort(randomised_positions, axis=axis)
    randomised = multiarray[np.arange(multiarray.shape[0])[:, None], ids]
    return randomised


# def sim_cds_seqs(codons_spaced, starts, stops, seed=None):
#     """
#     Given a list of sequences, shuffle the order of codons
#
#     Args:
#         seq_list (list): list containing sequences as strings
#         seed (list): if set, set the randomisation seed
#
#     Returns:
#         randomised_seqs (list): list of sequences with codons shuffled
#     """
#
#     if seed:
#         # set the seed
#         np.random.seed(seed)
#     # now randomise the set of codons at once
#     randomised_codons = randomise_multiarray(codons_spaced)
#     # join the sequences back together
#     randomised_seqs = ["{0}{1}{2}".format(starts[i], "".join(codon_set), stops[i]) for i, codon_set in enumerate(randomised_codons)]
#
#     return randomised_seqs


def sim_cds_seqs(codon_list, starts, stops, seed=None):
    """
    Given a list of internal codons, starts and stops, generate shuffled CDS sequences.

    Args:
        codon_list (list): list containing internal codons
        starts (list): list of start codons
        stops (list): list of stop codons
        seed (int): if set, set the randomisation seed

    Returns:
        randomised_seqs (list): list of sequences with codons shuffled
    """

    if seed:
        np.random.seed(seed)

    if isinstance(codon_list, list):
        randomised_seqs = []
        for i, codon_set in enumerate(codon_list):
            np.random.shuffle(codon_set)
            randomised_seqs.append("{0}{1}{2}".format(starts[i], "".join(codon_set), stops[i]))
        return randomised_seqs

    if isinstance(codon_list, dict):
        randomised_seqs = {}
        for name in codon_list:
            np.random.shuffle(codon_list[name])
            randomised_seqs[name] = "{0}{1}{2}".format(starts[name], "".join(codon_list[name]), stops[name])
        return randomised_seqs



def sim_cds_seqs_stop_counts(simulations, seqs, temp_dir, seeds=None):
    """
    Shuffle coding sequences and count the number of stop codons present.

    Args:
        simulations (list): list of simluations to iterate over
        seqs (dict): a dictionary of sequences
        temp_dir (str): a temporary directory to hold the outputs of the simulations
        seeds (list): list of seeds to be used for the randomisations

    Returns:
        temp_files (list): list containing file paths to simulation outputs
    """

    temp_files = []
    if len(simulations):

        codon_regex = re.compile(".{3}")
        codon_list = [re.findall(codon_regex, seq) for seq in seqs]
        # get a list of start and stop codons
        starts = [i[0] for i in codon_list]
        stops = [i[-1] for i in codon_list]
        # get a list of internal codons in the list
        codon_list = [codon_set[1:-1] for codon_set in codon_list]

        for sim_no, simulation in enumerate(simulations):
            # set the seed
            set_seed(seeds, simulation)
            # print the simulation number out
            gen.print_simulation(sim_no+1, simulations)
            # get the randomised seqs
            randomised_seqs = sim_cds_seqs(codon_list, starts, stops)
            sim_stop_counts = seqo.get_stop_counts([seq[:-3] for seq in randomised_seqs])
            temp_file = "{0}/{1}.{2}.txt".format(temp_dir, random.random(), simulation+1)
            temp_files.append(temp_file)
            with open(temp_file, "w") as outfile:
                outfile.write("{0}".format(",".join(gen.stringify([np.divide(count, len(randomised_seqs[i])) for i, count in enumerate(sim_stop_counts)]))))

    return temp_files


def get_exon_info(bed_file):

    exon_info = collections.defaultdict(lambda: collections.defaultdict(lambda: []))
    lines = gen.read_many_fields(bed_file, "\t")
    for line in lines:
        transcript = line[0]
        exons = [int(i) for i in line[1].split(",")]
        starts = [int(i) for i in line[2].split(",")]
        lengths = [int(i) for i in line[3].split(",")]
        frames = [int(i) for i in line[4].split(",")]
        for i, exon_id in enumerate(exons):
            exon_info[transcript][exon_id] = [starts[i], lengths[i], frames[i]]
    return exon_info


def get_exon_stop_densities(seqs, exon_info):


    stops = ["TAA", "TAG", "TGA"]

    densities = collections.defaultdict(lambda: collections.defaultdict(lambda: []))
    core_densities = []

    for name in sorted(seqs):
        cds_seq = seqs[name]

        exon_densities = collections.defaultdict(lambda: collections.defaultdict(lambda: []))
        exon_core_densities = []

        for exon_id in exon_info[name]:
            info = exon_info[name][exon_id]
            start = info[0]
            length = info[1]
            frame = info[2]

            exon_seq = cds_seq[start:start+length]
            half = int(len(exon_seq) / 2)

            region_size = 20

            regions = list(range(0, min(half, 140), region_size))

            for region in regions:
                # get the regions that arent at the ends of the exon
                # and not at the core
                if region != 0 and region < 120:
                    five_prime_seq = exon_seq[region:region+region_size]
                    three_prime_seq = exon_seq[length-region:length-region+region_size]
                    five_prime_count = 0
                    three_prime_count = 0
                    # get just the normal stop codons not at ends for regions
                    for i in range(0, len(five_prime_seq)):
                        codon = five_prime_seq[i:i+3]
                        if codon in stops:
                            five_prime_count += 3
                    for i in range(0, len(three_prime_seq)):
                        codon = three_prime_seq[i:i+3]
                        if codon in stops:
                            three_prime_count += 3
                    # # now get the counts at the 5' end of the region
                    # five_prime_region_five_prime_end = exon_seq[region-2:region+2]
                    # three_prime_region_five_prime_end = exon_seq[length-region-2:length-region+2]
                    # if five_prime_region_five_prime_end[:3] in stops:
                    #     five_prime_count += 1
                    # if five_prime_region_five_prime_end[1:] in stops:
                    #     five_prime_count += 2
                    # if three_prime_region_five_prime_end[:3] in stops:
                    #     three_prime_count += 1
                    # if three_prime_region_five_prime_end[1:] in stops:
                    #     three_prime_count += 2
                    # # now get the counts at the 3' end of the region
                    # five_prime_region_three_prime_end = exon_seq[region+region_size-2:region+region_size+2]
                    # three_prime_region_three_prime_end = exon_seq[length-region+region_size-2:length-region+region_size+2]
                    # if five_prime_region_three_prime_end[:3] in stops:
                    #     five_prime_count += 2
                    # if five_prime_region_five_prime_end[1:] in stops:
                    #     five_prime_count += 1
                    # if three_prime_region_five_prime_end[:3] in stops:
                    #     three_prime_count += 2
                    # if three_prime_region_five_prime_end[1:] in stops:
                    #     three_prime_count += 1

                    five_prime_density = np.divide(five_prime_count, region_size)
                    three_prime_density = np.divide(three_prime_count, region_size)
                    exon_densities[5][region].append(five_prime_density)
                    exon_densities[3][region].append(three_prime_density)

                # get the core region
                if region >= 120:
                    core_seq = exon_seq[region:length-region]
                    core_count = 0
                    for i in range(0, len(core_seq), 3):
                        codon = core_seq[i:i+3]
                        if codon in stops:
                            core_count += 3
                    # five_prime_core = exon_seq[region-2:region+2]
                    # if five_prime_core[:3] in stops:
                    #     core_count += 1
                    # if five_prime_core[1:] in stops:
                    #     core_count += 2
                    # three_prime_core = exon_seq[length-region-2:length-region+2]
                    # if three_prime_core[:3] in stops:
                    #     core_count += 2
                    # if three_prime_core[1:] in stops:
                    #     core_count += 1
                    core_density = np.divide(core_count, len(core_seq))
                    exon_core_densities.append(core_density)

                # get the five prime exon end
                five_prime_end = exon_seq[:region_size]
                five_prime_end_count = 0
                for i in range(0, len(five_prime_end), 3):
                    codon = five_prime_end[i:i+3]
                    if codon in stops:
                        five_prime_end_count += 3
                # five_prime_start = cds_seq[start-2:start+2]
                # if five_prime_start[:3] in stops:
                #     five_prime_end_count += 1
                # if five_prime_start[1:] in stops:
                #     five_prime_end_count += 2
                # five_prime_end = exon_seq[region_size-2:region_size+2]
                # if five_prime_end[:3] in stops:
                #     five_prime_end_count += 2
                # if five_prime_end[1:] in stops:
                #     five_prime_end_count += 1
                five_prime_end_density = np.divide(five_prime_end_count, region_size)
                exon_densities[5][0].append(five_prime_end_density)

                # get the three prime exon end
                three_prime_end = exon_seq[start+length-region_size:start+length]
                three_prime_end_count = 0
                for i in range(0, len(three_prime_end), 3):
                    codon = three_prime_end[i:i+3]
                    if codon in stops:
                        three_prime_end_count += 3
                # three_prime_end_start = exon_seq[start+length-region_size-2:start+length-region_size+2]
                # if three_prime_end_start[:3] in stops:
                #     three_prime_end_count += 1
                # if three_prime_end_start[1:] in stops:
                #     three_prime_end_count += 2
                # three_prime_end_end = cds_seq[start+length-2:start+length+2]
                # if three_prime_end_end[:3] in stops:
                #     three_prime_end_count += 2
                # if three_prime_end_end[1:] in stops:
                #     three_prime_end_count += 1
                three_prime_end_density = np.divide(three_prime_end_count, region_size)
                exon_densities[3][0].append(three_prime_end_density)

        for end in exon_densities:
            for region in exon_densities[end]:
                densities[end][region].append(np.mean(exon_densities[end][region]))
        core_densities.append(np.mean(exon_core_densities))


    return densities, core_densities

def sim_exon_stop_density(simulations, seqs, exon_positions_bed, temp_dir, seeds=None):
    """
    Shuffle coding sequences and count the number of stop codons present.

    Args:
        simulations (list): list of simluations to iterate over
        seqs (dict): a dictionary of sequences
        exons_positions_bed (dict): a dictionary of sequences
        temp_dir (str): a temporary directory to hold the outputs of the simulations
        seeds (list): list of seeds to be used for the randomisations

    Returns:
        temp_files (list): list containing file paths to simulation outputs
    """

    temp_files = []
    if len(simulations):

        exon_info = get_exon_info(exon_positions_bed)

        codon_regex = re.compile(".{3}")
        codon_list = {name: re.findall(codon_regex, seqs[name]) for name in seqs}
        # get a list of start and stop codons
        starts = {name: codon_list[name][0] for name in codon_list}
        stops = {name: codon_list[name][-1] for name in codon_list}
        # get a list of internal codons in the list
        codon_list = {name: codon_list[name][1:-1] for name in codon_list}

        for sim_no, simulation in enumerate(simulations):
            # set the seed
            set_seed(seeds, simulation)
            # print the simulation number out
            gen.print_simulation(sim_no+1, simulations)
            # get the randomised seqs
            randomised_seqs = sim_cds_seqs(codon_list, starts, stops)
            sim_densities, sim_core_densities = get_exon_stop_densities(randomised_seqs, exon_info)
            temp_file = "{0}/{1}.{2}.txt".format(temp_dir, random.random(), simulation+1)
            temp_files.append(temp_file)
            with open(temp_file, "w") as outfile:
                for end in sim_densities:
                    for region in sim_densities[end]:
                        outfile.write("{0},{1},{2}\n".format(end, region, ",".join(gen.stringify(sim_densities[end][region]))))
                outfile.write("core,{0}\n".format(",".join(gen.stringify(sim_core_densities))))

    return temp_files


def sim_introns(seqs):

    randomised_seqs = []
    for seq in seqs:
        np.random.shuffle(seq)
        randomised_seqs.append("".join(seq))
    return randomised_seqs


def sim_intron_seqs_stop_density(simulations, seqs, temp_dir, seeds=None):
    """
    Shuffle intron sequences and count the number of stop codons present.

    Args:
        simulations (list): list of simluations to iterate over
        seqs (dict): a dictionary of sequences
        temp_dir (str): a temporary directory to hold the outputs of the simulations
        seeds (list): list of seeds to be used for the randomisations

    Returns:
        temp_files (list): list containing file paths to simulation outputs
    """

    temp_files = []
    if len(simulations):

        nt_list = [list(seq) for seq in seqs]

        for sim_no, simulation in enumerate(simulations):
            # set the seed
            set_seed(seeds, simulation)
            # print the simulation number out
            gen.print_simulation(sim_no+1, simulations)
            # get the randomised seqs
            randomised_seqs = sim_introns(nt_list)
            sim_stop_counts = seqo.get_stop_counts(randomised_seqs)
            temp_file = "{0}/{1}.{2}.txt".format(temp_dir, random.random(), simulation+1)
            temp_files.append(temp_file)
            with open(temp_file, "w") as outfile:
                outfile.write("{0}".format(",".join(gen.stringify([np.divide(count, len(randomised_seqs[i])) for i, count in enumerate(sim_stop_counts)]))))

    return temp_files


def simulate_sequence_stop_density(simulations, sequences, nucleotide_probabilities, dinucleotide_probabilities, output_dir, seeds = None, seq_seeds = None):

    outputs = []
    seq_ids = [i for i in sorted(sequences)]
    seq_list = [sequences[i] for i in sorted(sequences)]

    for sim_no, simulation in enumerate(simulations):
        # set the seed
        set_seed(seeds, simulation)
        # print the simulation number out
        gen.print_simulation(sim_no+1, simulations)
        simulated_seqs = []
        # Generate a list of nucleotide matched sequences
        for i, seq in enumerate(seq_list):
            generated = False
            seq_seed = get_seq_seed(seq_seeds, sim_no, i)
            while not generated:
                sim_seq = seqo.generate_nt_matched_seq(seq, dinucleotides, dinucleotide_probabilities, nucleotides, nucleotide_probabilities, seed = seq_seed)
                if sim_seq not in simulated_seqs:
                    generated = True
                    simulated_seqs.append(sim_seq)

        sim_densities = [seqo.calc_stop_density(i) for i in simulated_seqs]
        temp_file = "{0}/{1}.bed".format(output_dir, random.random())
        with open(temp_file, "w") as outfile:
            [outfile.write("{0}\t{1}\n".format(seq_ids[i], sim_density)) for i, sim_density in enumerate(sim_densities)]
        outputs.append(temp_file)

    return outputs


    # dinucleotide_choices = [dn for dn in sorted(dinucleotide_content)]
    # dinucleotide_probabilities = [dinucleotide_content[dn] for dn in sorted(dinucleotide_content)]
    # nucleotide_choices = [n for n in sorted(nucleotide_content)]
    # nucleotide_probabilities = [nucleotide_content[n] for n in sorted(nucleotide_content)]
    #
    # output_files = []
    #
    # if len(simulations):
    #     for sim_no, simulation in enumerate(simulations):
    #         # set the seed
    #         set_seed(seeds, simulation)
    #
    #         # print the simulation number out
    #         gen.print_simulation(sim_no+1, simulations)
    #
    #         simulated_seqs = []
    #         # Generate a list of nucleotide matched sequences
    #         for i, seq in enumerate(seqs):
    #             generated = False
    #             seq_seed = get_seq_seed(seq_seeds, sim_no, i)
    #             while not generated:
    #                 sim_seq = seqo.generate_nt_matched_seq(seq, dinucleotide_choices, dinucleotide_probabilities, nucleotide_choices, nucleotide_probabilities, seed=seq_seed)
    #                 if sim_seq not in simulated_seqs:
    #                     generated = True
    #                     simulated_seqs.append(sim_seq)
    #
    #         output_file = "{0}/{1}.txt".format(output_directory, random.random())
    #         output_files.append(output_file)
    #         with open(output_file, "w") as outfile:
    #             for i, seq in enumerate(simulated_seqs):
    #                 outfile.write(">{0}\n{1}\n".format(i+1,seq))
    #
    # return output_files
