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

    arguments = ["output_directory", "motif_file", "simulations", "controls_directory", "exons_fasta", "motif_stop_codon_densities", "motif_codon_densities", "motif_densities_exon_dinucleotides", "generate_motif_controls", "match_density", "match_subs"]

    description = ""
    args = gen.parse_arguments(description, arguments, opt_flags=[2,3,4], flags = [5,6,7,8,9,10])
    output_directory, motif_file, simulations, controls_directory, exons_fasta, motif_stop_codon_densities, motif_codon_densities, motif_densities_exon_dinucleotides, generate_motif_controls, match_density, match_subs = args.output_directory, args.motif_file, args.simulations, args.controls_directory, args.exons_fasta, args.motif_stop_codon_densities, args.motif_codon_densities, args.motif_densities_exon_dinucleotides, args.generate_motif_controls, args.match_density, args.match_subs

    # interger the simulations
    if simulations:
        simulations = int(simulations)

    # create the global output directory
    global_output_directory = "{0}/motif_tests".format(output_directory)
    gen.create_output_directories(global_output_directory)

    # if we want to generate the controls
    if generate_motif_controls:
        simopc.generate_motif_dinucleotide_controls(motif_file, simulations, output_directory, match_density = match_density, match_subs = match_subs)

    if motif_stop_codon_densities:
        # create a local output directory
        local_output_directory = "{0}/motif_stop_density_simulations".format(global_output_directory)
        gen.create_output_directories(local_output_directory)
        # output filepath
        output_file = "{0}/{1}_stop_codon_densities.csv".format(local_output_directory, motif_file.split("/")[-1].split(".")[0])
        # run if we need some more controls
        if simulations > len(os.listdir(controls_directory)):
            print("Please create more simulants...")
            raise Exception
        # # calculate densities
        mtop.motif_stop_codon_densities(motif_file, controls_directory, simulations, output_file)

    # output_file = "{0}/{1}_stop_codon_densities_exon_dinucleotides.csv".format(global_output_directory, motif_file.split("/")[-1].split(".")[0])
    # exon_controls_directory = "{0}/dinucleotide_controls/exon_dinucleotide_controls".format(output_directory)
    # gen.create_output_directories(exon_controls_directory)
    # if motif_densities_exon_dinucleotides:
    #     if clean_run:
    #         gen.remove_directory(exon_controls_directory)
    #     gen.create_output_directories(exon_controls_directory)
    #
    #     if simulations > len(os.listdir(exon_controls_directory)):
    #         simopc.generate_exon_dinucleotide_controls(motif_file, exons_fasta, simulations, exon_controls_directory, match_density = False)
    #     # calculate densities
    #     mtop.motif_stop_codon_densities(motif_file, exon_controls_directory, simulations, output_file)

    if motif_codon_densities:
        local_output_directory = "{0}/codon_combination_densities".format(global_output_directory)
        gen.create_output_directories(local_output_directory)
        # get all the possible sets of 3 unique codon combinations
        codon_combinations_file = "{0}/codon_combinations.txt".format(local_output_directory)
        if not os.path.isfile(codon_combinations_file):
            seqo.generate_all_motif_combinations(stops, codon_combinations_file)

        output_file = "{0}/{1}_codon_combination_densities.csv".format(local_output_directory, motif_file.split("/")[-1].split(".")[0])
        if simulations > len(os.listdir(controls_directory)):
            gen.remove_directory(controls_directory)
            simopc.generate_motif_controls(motif_file, simulations, controls_directory, match_density = False)
        mtop.motif_codon_densities(motif_file, codon_combinations_file, controls_directory, simulations, output_file)



if __name__ == "__main__":
    main()
