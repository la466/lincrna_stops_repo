import generic as gen
import containers as cont
import sim_ops_containers as simopc
import sequence_ops as sequo
import ops_containers as opsc
import file_ops as fo
import random
import time
import os

def main():

    arguments = ["genome_gtf", "genome_fasta", "ortholog_gtf", "ortholog_fasta", "ensembl_links", "output_directory", "extract_sequences", "extract_exons", "extract_introns", "extract_coding_exons", "clean_run"]

    description = ""
    args = gen.parse_arguments(description, arguments, flags = [6,7,8,9,10], ints=[])
    genome_gtf, genome_fasta, ortholog_gtf, ortholog_fasta, ensembl_links, output_directory, extract_sequences, extract_exons, extract_introns, extract_coding_exons, clean_run = args.genome_gtf, args.genome_fasta, args.ortholog_gtf, args.ortholog_fasta, args.ensembl_links, args.output_directory, args.extract_sequences, args.extract_exons, args.extract_introns, args.extract_coding_exons, args.clean_run

    # set a start time
    start = time.time()

    # create the output_directory if it doenst already exist
    gen.create_output_directories(output_directory)

    # get the sequences
    if extract_sequences:
        # input_file1 = gtf genome 1, genome_fasta = genome fasta 1, ortholog_gtf = gtf genome 2, ortholog_fasta = genome fasta 2, ensembl_links = orthlogs file
        cont.extract_clean_sequences(genome_gtf, genome_fasta, ortholog_gtf, ortholog_fasta, ensembl_links, output_directory, clean_run = clean_run)

    full_exon_file = "{0}/genome_sequences/human/human.exons.bed".format(output_directory)
    if extract_exons:
        cont.extract_exons(genome_gtf, genome_fasta, output_directory, full_exon_file, clean_run = clean_run)
        if not os.path.isfile(full_exon_file) or clean_run:
            sequo.clean_feature_file(full_exon_file)

    exons_bed = "{0}/genome_sequences/{1}/{1}.cds.clean_filtered_exons.bed".format(output_directory, "human")
    coding_exons_bed = "{0}/genome_sequences/{1}/{1}.cds.clean_coding_exons.bed".format(output_directory, "human")
    if extract_coding_exons:
        sequo.get_coding_exons(full_exon_file, exons_bed, coding_exons_bed)

    non_coding_exons_bed = "{0}/genome_sequences/{1}/{1}.cds.clean_non_coding_exons.bed".format(output_directory, "human")
    if extract_non_coding_exons:
        sequo.get_non_coding_exons(full_exon, exons_bed, non_coding_exons_bed)

    intron_bed = "{0}/genome_sequences/human/human.introns.bed".format(output_directory)
    intron_fasta = "{0}/genome_sequences/human/human.introns.fasta".format(output_directory)
    if extract_introns:
        sequo.extract_introns(coding_exons_bed, intron_bed, intron_fasta, genome_fasta, clean_run = clean_run)


if __name__ == "__main__":
    main()
