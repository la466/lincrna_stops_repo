import generic as gen
import seq_ops as seqo
import ops
import numpy as np
import random

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


def sim_stop_counts_mm(simulations, seqs, temp_dir, seeds=None, seq_seeds=None):
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

    nts = ["A", "C", "G", "T"]
    # get the dicnucleotide and nucleotide content of the sequences
    # have to get here because multiprocessing cant pickle beforehand
    unique_seqs = [seqs[i] for i in seqs]
    dinucleotide_probs = seqo.get_dinucleotide_probabilities_markov(unique_seqs)
    # nucleotide_probs = seqo.get_nucleotide_probabilities_markov(unique_seqs)

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
                sim_seq = seqo.generate_nt_matched_seq_mm(seq, nts, dinucleotide_probs, seed=seq_seed)
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
