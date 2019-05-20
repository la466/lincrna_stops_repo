import generic as gen
import containers as cont
import sim_ops_containers as simopc
import sequence_ops as sequo
import seq_ops as seqo
import ops_containers as opsc
import file_ops as fo
import random
import time
import os

def main():

    arguments = [
        "output_directory", "genome_gtf", "genome_fasta", "ortholog_gtf", "ortholog_fasta", "input_file", "genome_fasta", "mapping_file", "codes_file", "ensembl_links", "extract_protein_coding", "extract_exons", "extract_introns", "extract_coding_exons", "extract_non_coding_exons", "extract_non_transcribed_regions", "extract_lincrna_seqs", "clean_run"
    ]

    description = ""
    args = gen.parse_arguments(description, arguments, opt_flags = [1,2,3,4,5,6,7,8,9], flags = [10,11,12,13,14,15,16,17])
    output_directory,  genome_gtf,  genome_fasta,  ortholog_gtf,  ortholog_fasta,  input_file,  genome_fasta,  mapping_file,  codes_file,  ensembl_links,  extract_protein_coding,  extract_exons,  extract_introns,  extract_coding_exons,  extract_non_coding_exons,  extract_non_transcribed_regions,  extract_lincrna_seqs,  clean_run = args.output_directory, args.genome_gtf, args.genome_fasta, args.ortholog_gtf, args.ortholog_fasta, args.input_file, args.genome_fasta, args.mapping_file, args.codes_file, args.ensembl_links, args.extract_protein_coding, args.extract_exons, args.extract_introns, args.extract_coding_exons, args.extract_non_coding_exons, args.extract_non_transcribed_regions, args.extract_lincrna_seqs, args.clean_run

    # set a start time
    start = time.time()

    # create the output_directory if it doenst already exist
    gen.create_output_directories(output_directory)

    # get the sequences
    if extract_protein_coding:
        # input_file1 = gtf genome 1, genome_fasta = genome fasta 1, ortholog_gtf = gtf genome 2, ortholog_fasta = genome fasta 2, ensembl_links = orthlogs file
        cont.extract_clean_sequences(genome_gtf, genome_fasta, ortholog_gtf, ortholog_fasta, ensembl_links, output_directory, clean_run = clean_run)

    full_exon_file = "{0}/genome_sequences/human/human.exons.bed".format(output_directory)
    if extract_exons:
        cont.extract_exons(genome_gtf, genome_fasta, output_directory, full_exon_file, clean_run = clean_run)
        sequo.clean_feature_file(full_exon_file)

    exons_bed = "{0}/genome_sequences/{1}/{1}.cds.clean_filtered_exons.bed".format(output_directory, "human")
    coding_exons_bed = "{0}/genome_sequences/{1}/{1}.cds.clean_coding_exons.bed".format(output_directory, "human")
    coding_exons_fasta = "{0}/genome_sequences/{1}/{1}.cds.clean_coding_exons.fasta".format(output_directory, "human")
    if extract_coding_exons:
        sequo.get_coding_exon_coordinates(full_exon_file, exons_bed, coding_exons_bed)
        fo.fasta_from_intervals(coding_exons_bed, coding_exons_fasta, genome_fasta, names = True)

    if extract_non_coding_exons:
        non_coding_exons_bed = "{0}/genome_sequences/{1}/{1}.cds.clean_non_coding_exons.bed".format(output_directory, "human")
        non_coding_exons_fasta = "{0}/genome_sequences/{1}/{1}.cds.clean_non_coding_exons.fasta".format(output_directory, "human")
        sequo.get_non_coding_exon_coordinates(full_exon_file, exons_bed, non_coding_exons_bed)
        fo.fasta_from_intervals(non_coding_exons_bed, non_coding_exons_fasta, genome_fasta, names = True)


    if extract_introns:
        intron_bed = "{0}/genome_sequences/human/human.clean_introns.bed".format(output_directory)
        intron_fasta = "{0}/genome_sequences/human/human.clean_introns.fasta".format(output_directory)
        sequo.get_intron_coordinates(coding_exons_bed, intron_bed)
        fo.fasta_from_intervals(intron_bed, intron_fasta, genome_fasta, names = True)

    if extract_non_transcribed_regions:
        all_features_bed = "{0}/genome_sequences/human/human.all_features.bed".format(output_directory)
        non_transcribed_bed = "{0}/genome_sequences/human/human.non_transcribed.bed".format(output_directory)
        non_transcribed_fasta = "{0}/genome_sequences/human/human.non_transcribed.fasta".format(output_directory)
        seqo.get_non_transcribed_regions(genome_gtf, genome_fasta, all_features_bed, non_transcribed_bed, non_transcribed_fasta, output_directory)


    # extract sequences from source file
    if extract_lincrna_seqs:
        # set up the output fasta to contain the exon seqs
        lincrna_exons_bed = "{0}/lincRNA_exons.bed".format(output_directory)
        lincrna_exons_fasta = "{0}/lincRNA_exons.fasta".format(output_directory)
        lincrna_seqs_fasta = "{0}/lincRNA_seqs.fasta".format(output_directory)
        print("Extracting lincRNA seqs...")
        fo.extract_seqs(input_file, genome_fasta, lincrna_exons_bed, lincrna_exons_fasta, lincrna_seqs_fasta, mapping_file, codes_file, exclude_XY=True, hg38=hg38, NONCODE=NONCODE)

if __name__ == "__main__":
    main()
