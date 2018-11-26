import generic as gen
import containers as cont
import lincRNA_tests as ltests
import time
import os

def main():

    arguments = ["input_file1", "input_file2", "output_directory", "extract_sequences", "clean_run", "stop_density_test"]
    description = "Container for analysis on lincRNAs"
    args = gen.parse_arguments(description, arguments, flags = [3,4,5], ints=[])
    input_file1, input_file2, output_directory, extract_sequences, clean_run, stop_density_test = args.input_file1, args.input_file2, args.output_directory, args.extract_sequences, args.clean_run, args.stop_density_test

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

    if stop_density_test:
        ltests.stop_density_test(lincRNA_multi_exon_fasta, lincRNA_multi_exon_families)

if __name__ == "__main__":
    main()
