import generic as gen
import sim_ops_containers as simoc
import sequence_ops as sequo


# TODO: simulation to get the normalised densities of lincRNA stop codons
# TODO: compare stop codon density in exons of lincRNA

def stop_density_test(input_fasta, output_file, required_simulations, families_file = None):
    simoc.simulate_sequence_stop_density(input_fasta, output_file, required_simulations, families_file = families_file)
