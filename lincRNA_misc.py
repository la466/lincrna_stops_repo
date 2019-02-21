import generic as gen
import lincRNA_misc_ops as lmo
import conservation as cons
import file_ops as fo
import time
import os

import lincRNA_tests as ltests
import sim_ops_containers as simopc

def main():

    arguments = ["working_directory", "output_directory", "genome_path", "input_bed", "input_fasta", "clean_run", "extract_exon_intron_bed", "extract_exons", "extract_introns", "sort_by_exon_number", "build_transcripts", "extract_families", "orf_length_sim"]
    description = "Wrapper for miscellaneous operations on lincRNA"
    args = gen.parse_arguments(description, arguments, flags = [5,6,7,8,9,10,11,12,13], opt_flags = [2,3,4])
    working_directory, output_directory, genome_path, input_bed, input_fasta, clean_run, extract_exon_intron_bed, extract_exons, extract_introns, sort_by_exon_number, build_transcripts, extract_families, orf_length_sim = args.working_directory, args.output_directory, args.genome_path, args.input_bed, args.input_fasta, args.clean_run, args.extract_exon_intron_bed, args.extract_exons, args.extract_introns, args.sort_by_exon_number, args.build_transcripts, args.extract_families, args.orf_length_sim

    # create the directories
    gen.create_output_directories(working_directory)
    gen.create_output_directories(output_directory)

    # file paths
    exons_bed = "{0}/exons.bed".format(working_directory)
    single_exons_bed = "{0}/single_exons.bed".format(working_directory)
    multi_exons_bed = "{0}/multi_exons.bed".format(working_directory)
    exons_fasta = "{0}/exons.fasta".format(working_directory)
    single_exons_fasta = "{0}/single_exons.fasta".format(working_directory)
    multi_exons_fasta = "{0}/multi_exons.fasta".format(working_directory)
    introns_bed = "{0}/introns.bed".format(working_directory)
    introns_fasta = "{0}/introns.fasta".format(working_directory)
    transcript_sequences_fasta = "{0}/transcript_sequences.fasta".format(working_directory)
    multi_exon_transcript_sequences_fasta = "{0}/multi_exon_transcript_sequences.fasta".format(working_directory)
    multi_exon_blast_file = "{0}/multi_exons_blast_all_against_all.csv".format(working_directory)
    multi_exon_blast_database = "{0}/multi_exon_blast_all_against_all".format(working_directory)
    multi_exon_families_file = "{0}/multi_exon_families.txt".format(working_directory)

    # create the exons and introns files from bed
    if extract_exon_intron_bed:
        # copy the main file to the folder
        gen.copy_file(input_bed, "{0}/{1}".format(working_directory, input_bed.split("/")[-1]))
        # extract the features
        lmo.extract_bed_coordinates_block_format(input_bed, exons_bed, introns_bed)
    # get files for each
    if sort_by_exon_number:
        gen.check_files_exists([exons_bed])
        lmo.sort_by_exon_number(exons_bed, single_exons_bed, multi_exons_bed)

    # get exons
    if extract_exons:
        gen.check_files_exists([exons_bed])
        fo.fasta_from_intervals(exons_bed, exons_fasta, genome_path, names=True)
        # if the single exons bed file exists, get just the single exon sequences
        if os.path.isfile(single_exons_bed):
            lmo.sort_fasta_by_bed(single_exons_bed, exons_fasta, single_exons_fasta)
        # if the multi exons bed file exists, get just the multi exon sequences
        if os.path.isfile(multi_exons_bed):
            lmo.sort_fasta_by_bed(multi_exons_bed, exons_fasta, multi_exons_fasta)

    # get introns
    if extract_introns:
        gen.check_files_exists([introns_bed])
        fo.fasta_from_intervals(introns_bed, introns_fasta, genome_path, names=True)

    # build transcripts
    if build_transcripts:
        gen.check_files_exists([exons_fasta])
        lmo.build_transcripts(exons_fasta, transcript_sequences_fasta)
        # if the multi exons bed file exists, get just the multi exon sequences
        if os.path.isfile(multi_exons_bed):
            lmo.sort_fasta_by_bed(multi_exons_bed, transcript_sequences_fasta, multi_exon_transcript_sequences_fasta)

    # now group into paralagous families
    if extract_families:
        gen.check_files_exists([multi_exon_transcript_sequences_fasta])
        cons.filter_families(multi_exon_transcript_sequences_fasta, multi_exon_blast_file, multi_exon_families_file, database_path = multi_exon_blast_database, clean_run = clean_run)

    sim_orf_length_output_file = "{0}/orf_length_sim.csv".format(output_directory)
    sim_orf_length_z_file = "{0}/orf_length_sim_z.csv".format(output_directory)
    if orf_length_sim:
        simopc.sim_orf_length(multi_exon_transcript_sequences_fasta, 6, sim_orf_length_output_file)
        ltests.process_length_sim(sim_orf_length_output_file, sim_orf_length_z_file, families_file = multi_exon_families_file)
        # ltests.calculate_lengths(lincRNA_multi_exon_fasta, lincRNA_length_output_file, families_file = lincRNA_multi_exon_families)


if __name__ == "__main__":
    main()
