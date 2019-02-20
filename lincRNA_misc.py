import generic as gen
import lincRNA_misc_ops as lmo
import time
import os

def main():

    arguments = ["input_directory", "output_directory", "input_bed", "input_fasta", "extract_exon_intron_fasta", "extract_sequences"]
    description = "Wrapper for miscellaneous operations on lincRNA"
    args = gen.parse_arguments(description, arguments, flags = [4,5], opt_flags = [2,3])
    input_directory, output_directory, input_bed, input_fasta, extract_exon_intron_fasta, extract_sequences = args.input_directory, args.output_directory, args.input_bed, args.input_fasta, args.extract_exon_intron_fasta, args.extract_sequences

    # create the directories
    gen.create_output_directories(input_directory)
    gen.create_output_directories(output_directory)

    # file paths
    exons_fasta = "{0}/exons.fasta".format(input_directory)
    introns_fasta = "{0}/introns.fasta".format(input_directory)
    sequences_fasta = "{0}/sequences.fasta".format(input_directory)

    # create the exons and introns files
    if extract_exon_intron_fasta:
        # copy the main file to the folder
        gen.copy_file(input_fasta, "{0}/{1}".format(output_directory, input_fasta.split("/")[-1]))
        # extract the features
        lmo.extract_exon_intron_fasta(input_fasta, exons_fasta, introns_fasta)

    # create the full length sequences from the exons
    if extract_sequences:
        lmo.extract_sequences(exons_fasta, sequences_fasta)




if __name__ == "__main__":
    main()
