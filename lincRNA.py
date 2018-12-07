import generic as gen
import containers as cont
import lincRNA_tests as ltests
import sim_ops_containers as simopc
import time
import os

def main():

    arguments = ["input_bed", "input_fasta", "output_directory", "required_simulations", "extract_sequences", "clean_run", "get_exon_dint_controls", "get_intron_dint_controls"]
    description = "Container for analysis on lincRNAs"
    args = gen.parse_arguments(description, arguments, flags = [4,5,6,7], ints=[3])
    input_bed, input_fasta, output_directory, required_simulations, extract_sequences, clean_run, get_exon_dint_controls, get_intron_dint_controls = args.input_bed, args.input_fasta, args.output_directory, args.required_simulations, args.extract_sequences, args.clean_run, args.get_exon_dint_controls, args.get_intron_dint_controls

    # set a start time
    start = time.time()

    # create the output_directory if it doenst already exist
    gen.create_output_directories(output_directory)

    lincRNA_single_exon_bed = "{0}/lincRNA.single_exon.bed".format(output_directory)
    lincRNA_single_exon_fasta = "{0}/lincRNA.single_exon.fasta".format(output_directory)
    lincRNA_single_exon_families = "{0}/lincRNA.single_exon_families.bed".format(output_directory)
    lincRNA_multi_exon_bed = "{0}/lincRNA.multi_exon.bed".format(output_directory)
    lincRNA_multi_exon_intron_bed = "{0}/lincRNA.multi_exon.introns.bed".format(output_directory)
    lincRNA_multi_exon_fasta = "{0}/lincRNA.multi_exon.fasta".format(output_directory)
    lincRNA_multi_exon_intron_fasta = "{0}/lincRNA.multi_exon.introns.fasta".format(output_directory)
    lincRNA_multi_exon_families = "{0}/lincRNA.multi_exon_families.bed".format(output_directory)
    # get the sequences
    if extract_sequences or not os.path.isfile(lincRNA_single_exon_fasta) or not os.path.isfile(lincRNA_multi_exon_fasta):
        cont.extract_lincRNA_sequences(input_bed, input_fasta, lincRNA_single_exon_bed, lincRNA_multi_exon_bed, lincRNA_single_exon_fasta, lincRNA_multi_exon_fasta, lincRNA_multi_exon_intron_bed, lincRNA_multi_exon_intron_fasta, lincRNA_single_exon_families, lincRNA_multi_exon_families, clean_run = clean_run)

    # multi_exon_sequence_density_simulation = "{0}/tests/lincRNA.multi_exon.full_sequence.stop_density_simulation.csv".format(output_directory)
    # if stop_density_test:
    #     ltests.stop_density_test(lincRNA_multi_exon_fasta, multi_exon_sequence_density_simulation, required_simulations, families_file = lincRNA_multi_exon_families)
    # single_exon_sequence_density_simulation = "{0}/tests/lincRNA.single_exon.full_sequence.stop_density_simulation.csv".format(output_directory)
    # if stop_density_test:
    #     ltests.stop_density_test(lincRNA_single_exon_fasta, single_exon_sequence_density_simulation, required_simulations, families_file = lincRNA_single_exon_families)

    exon_dint_control_directory = "{0}/exon_dinucleotide_controls".format(output_directory)
    if get_exon_dint_controls:
        simopc.generate_dint_controls(lincRNA_multi_exon_fasta, exon_dint_control_directory)

    intron_dint_control_directory = "{0}/intron_dinucleotide_controls".format(output_directory)
    if get_intron_dint_controls:
        simopc.generate_dint_intron_controls(lincRNA_multi_exon_intron_fasta, intron_dint_control_directory)

if __name__ == "__main__":
    main()
