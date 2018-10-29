import generic as gen
import seq_ops as seqo
import ops
import numpy as np
import random
import re

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
    dicnucleotide_probabilities = [dinucleotide_content[dn] for dn in sorted(dinucleotide_content)]
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
                            sim_seq, failed = seqo.generate_nt_matched_seq_frame(seq, dinucleotide_choices, dicnucleotide_probabilities, nucleotide_choices, nucleotide_probabilities, seq_frame=frame, seed=seq_seed)
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
    dicnucleotide_probabilities = [dinucleotide_content[dn] for dn in sorted(dinucleotide_content)]
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
                        sim_seq, failed = seqo.generate_nt_matched_seq_frame(seq, dinucleotide_choices, dicnucleotide_probabilities, nucleotide_choices, nucleotide_probabilities, seq_frame=frame, seed=seq_seed)
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
        dicnucleotide_content (dict): dictionary containing dinucleotide proportions
        temp_dir (str): a temporary directory to hold the outputs of the simulations
        seeds (list): list of seeds to be used for the randomisations

    Returns:
        temp_files (list): list containing file paths to simulation outputs
    """

    temp_files = []

    dinucleotide_choices = [dn for dn in sorted(dinucleotide_content)]
    dicnucleotide_probabilities = [dinucleotide_content[dn] for dn in sorted(dinucleotide_content)]
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
                sim_seq = seqo.generate_nt_matched_seq(seq, dinucleotide_choices, dicnucleotide_probabilities, nucleotide_choices, nucleotide_probabilities, seed=seq_seed)
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
        dicnucleotide_content (dict): dictionary containing dinucleotide proportions
        temp_dir (str): a temporary directory to hold the outputs of the simulations
        seeds (list): list of seeds to be used for the randomisations

    Returns:
        temp_files (list): list containing file paths to simulation outputs
    """

    temp_files = []

    dinucleotide_choices = [dn for dn in sorted(dinucleotide_content)]
    dicnucleotide_probabilities = [dinucleotide_content[dn] for dn in sorted(dinucleotide_content)]
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
                sim_seq = seqo.generate_nt_matched_seq(seq, dinucleotide_choices, dicnucleotide_probabilities, nucleotide_choices, nucleotide_probabilities, seed=seq_seed)
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
        dicnucleotide_content (dict): dictionary containing dinucleotide proportions
        temp_dir (str): a temporary directory to hold the outputs of the simulations
        seeds (list): list of seeds to be used for the randomisations

    Returns:
        temp_files (list): list containing file paths to simulation outputs
    """

    temp_files = []

    dinucleotide_choices = [dn for dn in sorted(dinucleotide_content)]
    dicnucleotide_probabilities = [dinucleotide_content[dn] for dn in sorted(dinucleotide_content)]
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
                sim_seq = seqo.generate_nt_matched_seq(seq, dinucleotide_choices, dicnucleotide_probabilities, nucleotide_choices, nucleotide_probabilities, seed=seq_seed)
                if sim_seq not in simulated_seqs:
                    generated = True
                    simulated_seqs.append(sim_seq)

        stop_counts = seqo.get_stop_counts(simulated_seqs)
        temp_file = "{0}/{1}.{2}.txt".format(temp_dir, random.random(), simulation+1)
        temp_files.append(temp_file)
        with open(temp_file, "w") as temp:
            temp.write(">{0}\n{1}\n".format(simulation+1, sum(stop_counts)))

    return temp_files


def generate_matched_gc_seqs(simulations, seqs_fasta, non_features_fasta, threshold, output_dir, seeds=None):

    temp_files = []

    if len(simulations):
        # get the sequences
        input_names, input_seqs = gen.read_fasta(seqs_fasta)
        input_seqs_list = {name: input_seqs[i] for i, name in enumerate(input_names)}
        genome_seqs = gen.read_fasta(non_features_fasta)[1]
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

    randomised_seqs = []
    for i, codon_set in enumerate(codon_list):
        np.random.shuffle(codon_set)
        randomised_seqs.append("{0}{1}{2}".format(starts[i], "".join(codon_set), stops[i]))
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
