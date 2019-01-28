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

    arguments = ["output_directory", "motif_file", "required_simulations", "clean_run", "exons_fasta", "motif_stop_codon_densities", "motif_codon_densities", "motif_densities_exon_dinucleotides"]

    description = ""
    args = gen.parse_arguments(description, arguments, flags = [3,5,6,7], ints=[2], opt_flags=[4])
    output_directory, motif_file, required_simulations, clean_run, exons_fasta, motif_stop_codon_densities, motif_codon_densities, motif_densities_exon_dinucleotides = args.output_directory, args.motif_file, args.required_simulations, args.clean_run, args.exons_fasta, args.motif_stop_codon_densities, args.motif_codon_densities, args.motif_densities_exon_dinucleotides

    test_output_directory = "{0}/motif_tests".format(output_directory)
    gen.create_output_directories(test_output_directory)


    motif_controls_directory = "{0}/dinucleotide_controls/{1}_dinucleotide_controls".format(output_directory, motif_file.split("/")[-1].split(".")[0])
    output_file = "{0}/{1}_stop_codon_densities.csv".format(test_output_directory, motif_file.split("/")[-1].split(".")[0])

    if motif_stop_codon_densities:
        gen.create_output_directories(motif_controls_directory)

        if required_simulations > len(os.listdir(motif_controls_directory)):
            gen.remove_directory(motif_controls_directory)
            simopc.generate_motif_dinucleotide_controls(motif_file, required_simulations, motif_controls_directory, match_density = False)
        # # calculate densities
        mtop.motif_stop_codon_densities(motif_file, motif_controls_directory, required_simulations, output_file)

    output_file = "{0}/{1}_stop_codon_densities_exon_dinucleotides.csv".format(test_output_directory, motif_file.split("/")[-1].split(".")[0])
    exon_controls_directory = "{0}/dinucleotide_controls/exon_dinucleotide_controls".format(output_directory)
    gen.create_output_directories(exon_controls_directory)
    if motif_densities_exon_dinucleotides:
        if clean_run:
            gen.remove_directory(exon_controls_directory)
        gen.create_output_directories(exon_controls_directory)

        if required_simulations > len(os.listdir(exon_controls_directory)):
            simopc.generate_exon_dinucleotide_controls(motif_file, exons_fasta, required_simulations, exon_controls_directory, match_density = False)
        # calculate densities
        mtop.motif_stop_codon_densities(motif_file, exon_controls_directory, required_simulations, output_file)


    output_file = "{0}/{1}_densities.csv".format(test_output_directory, motif_file.split("/")[-1].split(".")[0])
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
