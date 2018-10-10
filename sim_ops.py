import generic as gen
import seq_ops as seqo
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

            # print(sim_no, name, sim_seq)
        longest_orfs = seqo.get_longest_orfs(simulated_seqs)
        temp_file = "{0}/{1}.{2}.txt".format(temp_dir, random.random(), simulation+1)
        temp_files.append(temp_file)
        with open(temp_file, "w") as temp:
            for name in longest_orfs:
                temp.write(">{0}\n{1}\n".format(name, longest_orfs[name]))

    return temp_files
