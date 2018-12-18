import generic as gen
import containers as cont
import sim_ops_containers as simopc
import main_tests_ops as mto
import ops_containers as opsc
import seq_ops as seqo
import sequence_ops as sequo
import motif_tests_ops as mtop
import file_ops as fo
import time
from useful_motif_sets import stops
import os

def main():

    arguments = ["output_directory", "required_simulations", "motif_codon_densities"]

    description = ""
    args = gen.parse_arguments(description, arguments, flags = [2], ints=[1])
    output_directory, required_simulations, motif_codon_densities = args.output_directory, args.required_simulations, args.motif_codon_densities

    motif_file = "source_data/motif_sets/int3.txt"
    motif_controls_directory = "{0}/int3_dinucleotide_controls".format(output_directory)
    gen.create_output_directories(motif_controls_directory)

    test_output_directory = "{0}/motif_tests".format(output_directory)
    gen.create_output_directories(test_output_directory)

    output_file = "{0}/int3_densities.csv".format(test_output_directory)
    codon_combinations_file = "{0}/codon_combinations.txt".format(output_directory)
    if motif_codon_densities:

        if not os.path.isfile(codon_combinations_file):
            seqo.generate_all_motif_combinations(stops, codon_combinations_file)

        if required_simulations > len(os.listdir(motif_controls_directory)):
            gen.remove_directory(motif_controls_directory)
            simopc.generate_motif_dinucleotide_controls(motif_file, required_simulations, motif_controls_directory, match_density = False)

        mtop.motif_codon_densities(motif_file, codon_combinations_file, motif_controls_directory, required_simulations, output_file)



if __name__ == "__main__":
    main()
