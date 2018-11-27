import generic as gen
import containers as cont
import lincRNA_tests as ltests
import time
import os

def main():

    arguments = ["input_file1", "input_file2", "output_directory", "required_simulations", "extract_sequences", "clean_run", "stop_density_test"]
    description = "Container for analysis on lincRNAs"
    args = gen.parse_arguments(description, arguments, flags = [4,5,6], ints=[3])
    input_file1, input_file2, output_directory, required_simulations, extract_sequences, clean_run, stop_density_test = args.input_file1, args.input_file2, args.output_directory, args.required_simulations, args.extract_sequences, args.clean_run, args.stop_density_test

    # set a start time
    start = time.time()

    # create the output_directory if it doenst already exist
    gen.create_output_directories(output_directory)

    lincRNA_single_exon_bed = "{0}/lincRNA.single_exon.bed".format(output_directory)
    lincRNA_single_exon_fasta = "{0}/lincRNA.single_exon.fasta".format(output_directory)
    lincRNA_single_exon_families = "{0}/lincRNA.single_exon_families.bed".format(output_directory)
    lincRNA_multi_exon_bed = "{0}/lincRNA.multi_exon.bed".format(output_directory)
    lincRNA_multi_exon_fasta = "{0}/lincRNA.multi_exon.fasta".format(output_directory)
    lincRNA_multi_exon_families = "{0}/lincRNA.multi_exon_families.bed".format(output_directory)
    # get the sequences
    if extract_sequences or not os.path.isfile(lincRNA_single_exon_fasta) or not os.path.isfile(lincRNA_multi_exon_fasta):
        cont.extract_lincRNA_sequences(input_file1, input_file2, lincRNA_single_exon_bed, lincRNA_multi_exon_bed, lincRNA_single_exon_fasta, lincRNA_multi_exon_fasta, lincRNA_single_exon_families, lincRNA_multi_exon_families, clean_run = clean_run)

    # multi_exon_sequence_density_simulation = "{0}/tests/lincRNA.multi_exon.full_sequence.stop_density_simulation.csv".format(output_directory)
    # if stop_density_test:
    #     ltests.stop_density_test(lincRNA_multi_exon_fasta, multi_exon_sequence_density_simulation, required_simulations, families_file = lincRNA_multi_exon_families)
    single_exon_sequence_density_simulation = "{0}/tests/lincRNA.single_exon.full_sequence.stop_density_simulation.csv".format(output_directory)
    if stop_density_test:
        ltests.stop_density_test(lincRNA_single_exon_fasta, single_exon_sequence_density_simulation, required_simulations, families_file = lincRNA_single_exon_families)

if __name__ == "__main__":
    main()
