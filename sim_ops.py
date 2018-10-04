import generic as gen
import seq_ops as seqo


def sim_orf_length(seq_fasta, required_simulations, output_file):
    """
    Simulation to look at the lengths of ORFs in sequences
    compared with simulated null sequences

    Args:
        seq_fasta (str): path to fasta file containing sequences
        required_simulations (int): number of simulations to do
        output_file (str): path to output file
    """

    names, seqs = gen.read_fasta(seq_fasta)
    shortest_orfs = seqo.get_shortest_orfs(seqs)
    # print(seqs[0])
    print(shortest_orfs)
