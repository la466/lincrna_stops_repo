import generic as gen
import ops as ops
import file_ops as fo
import seq_ops as seqo
import sim_ops_containers as simoc
import numpy as np
import itertools as it
import os
import re
import collections

def get_non_coding_exons(genome_fasta, gtf, coding_exons_bed, non_coding_exons_bed, coding_exons_fasta, non_coding_exons_fasta, output_directory, clean_run=None):
    """
    Extract the coding and non coding exons.

    Args:
        genome_fasta (str): path to genome fasta
        gtf (str): path to genome gtf
        coding_exons_bed (str): path to coding exons bed output file
        non_coding_exons_bed (str): path to non coding exons bed output file
        coding_exons_fasta (str): path to coding exons fasta output file
        non_coding_exons_fasta (str): path to non coding exons fasta output file
        output_directory (str): path to the output directory
        clean_run (bool): if true, run all the steps
    """

    # get the required features
    full_exons_bed = "{0}/full_exons.bed".format(output_directory)
    full_cds_bed = "{0}/full_cds.bed".format(output_directory)
    full_stops_bed = "{0}/full_stop_codons.bed".format(output_directory)
    if clean_run or not os.path.isfile(full_exons_bed) or not os.path.isfile(full_cds_bed) or not os.path.isfile(full_stops_bed):
        print("Extracting features from .gtf file...")
        ops.extract_gtf_features(gtf, ["exon"], full_exons_bed)
        ops.extract_gtf_features(gtf, ["CDS"], full_cds_bed)
        ops.extract_gtf_features(gtf, ["stop_codon"], full_stops_bed)

    # get the fasta sequences of the features
    full_exons_fasta = "{0}/full_exons.fasta".format(output_directory)
    full_cds_fasta = "{0}/full_cds.fasta".format(output_directory)
    full_stops_fasta = "{0}/full_stop_codons.fasta".format(output_directory)
    if clean_run or not os.path.isfile(full_exons_fasta) or not os.path.isfile(full_cds_fasta) or not os.path.isfile(full_stops_fasta):
        print("Converting bed to fasta...")
        fo.fasta_from_intervals(full_exons_bed, full_exons_fasta, genome_fasta, names=True)
        fo.fasta_from_intervals(full_cds_bed, full_cds_fasta, genome_fasta, names=True)
        fo.fasta_from_intervals(full_stops_bed, full_stops_fasta, genome_fasta, names=True)

    # build the coding sequences
    cds_sequences_fasta = "{0}/cds_sequences.fasta".format(output_directory)
    if clean_run or not os.path.isfile(cds_sequences_fasta):
        print("Building coding sequences...")
        ops.build_sequences(full_cds_fasta, full_stops_fasta, cds_sequences_fasta)

    # filter the coding sequences
    quality_filtered_cds_sequences_fasta = "{0}/quality_filtered_cds.fasta".format(output_directory)
    if clean_run or not os.path.isfile(quality_filtered_cds_sequences_fasta):
        print("Filtering coding sequences...")
        ops.filter_coding_sequences(cds_sequences_fasta, quality_filtered_cds_sequences_fasta)

    # filter original cds bed file to get filtered bed entries
    quality_filtered_cds_bed = "{0}/quality_filtered_cds.bed".format(output_directory)
    if clean_run or not os.path.isfile(quality_filtered_cds_bed):
        print("Filtering bed file to only contain entries from filtered coding sequences...")
        ops.filter_bed_from_fasta(full_cds_bed, quality_filtered_cds_sequences_fasta, quality_filtered_cds_bed)

    # get one transcript per gene
    unique_filtered_transcripts_fasta = "{0}/unique_filtered_cds.fasta".format(output_directory)
    if clean_run or not os.path.isfile(unique_filtered_transcripts_fasta):
        print("Getting unique transcripts...")
        ops.get_unique_transcripts(quality_filtered_cds_bed, quality_filtered_cds_sequences_fasta, unique_filtered_transcripts_fasta)

    # # # ***
    # # # do I want to filter to one per gene family here ?
    # # # ***

    # from these transcripts, get the exons that make them up
    final_filtered_cds_bed = "{0}/final_filtered_cds.bed".format(output_directory)
    final_filtered_cds_fasta = "{0}/final_filtered_cds.fasta".format(output_directory)
    if not os.path.isfile(final_filtered_cds_bed):
        print("Getting bed entries for unique filtered transcripts...")
        ops.filter_bed_from_fasta(quality_filtered_cds_bed, unique_filtered_transcripts_fasta, final_filtered_cds_bed)
        fo.fasta_from_intervals(final_filtered_cds_bed, final_filtered_cds_fasta, genome_fasta, names=True)

    # get the exon entries for transcripts that remain in the bed entries
    final_filtered_exons_bed = "{0}/final_filtered_exons.bed".format(output_directory)
    final_filtered_exons_fasta = "{0}/final_filtered_exons.fasta".format(output_directory)
    if not os.path.isfile(final_filtered_exons_bed):
        print("Getting the exon entries for corresponding cds bed entries...")
        ops.filter_bed_transcript_id(full_exons_bed, final_filtered_cds_bed, final_filtered_exons_bed)
        fo.fasta_from_intervals(final_filtered_exons_bed, final_filtered_exons_fasta, genome_fasta, names=True)

    # getting the coding of the exons
    print("Getting coding and non coding exons...")
    ops.get_exon_coding(final_filtered_exons_bed, quality_filtered_cds_bed, final_filtered_cds_bed, non_coding_exons_bed, coding_exons_bed)
    # get the sequences for the coding and non coding exons
    fo.fasta_from_intervals(non_coding_exons_bed, non_coding_exons_fasta, genome_fasta, names=True)
    fo.fasta_from_intervals(coding_exons_bed, coding_exons_fasta, genome_fasta, names=True)


def lincRNA_expression(fasta_file, link_file, expression_file, output_file):
    """
    Wrapper to look at expression of lincRNA

    Args:
        fasta_file (str): path to the fasta file containing exon sequences
        link_file (str): path to link file matching transcripts to genes
        expression_file (str): path to file containing expression data
        output_file (str): path to output file
    """

    # links = gen.read_many_fields(link_file, "\t")
    # links = {link[0]: link[1] for link in links}

    names, seqs = gen.read_fasta(fasta_file)
    expression = gen.read_many_fields(expression_file, "\t")

    tissues = expression[2][2:]
    expression_lines = {line[0]: line[2:] for line in expression if "TCONS" in line[0]}

    with open(output_file, "w") as outfile:
        outfile.write("exon,{0}\n".format(",".join(tissues)))
        for exon in names:
            transcript = exon.split(".")[0]
            # if transcript in links:
                # loc = links[transcript]
            if transcript in expression_lines:
                outfile.write("{0},{1}\n".format(exon, ",".join(expression_lines[transcript])))


def non_coding_exons(genome_fasta, gtf, output_directory, sim_stop_count_non_coding_output_file, sim_stop_count_coding_output_file, required_simulations, clean_run=None):
    """
    Wrapper for looking at non coding exons

    Args:
        genome_fasta (str): path to the genome fasta file
        gtf (str): path to the genome gtf file
        output_directory (str): path to the output directory
        sim_stop_count_non_coding_output_file (str): path to the output file for the non coding exons stop count simulation
        sim_stop_count_coding_output_file (str): path to the output file for the coding exons stop count simulation
        required_simulations (int): number of simulations to run
        clean_run (bool): if true, remove the processing sequence files and start fresh
    """

    sequence_output_directory = "{0}/sequence_files".format(output_directory)
    if clean_run:
        gen.remove_directory(sequence_output_directory)

    # create the output directory if it doesnt exist
    gen.create_output_directories(sequence_output_directory)

    # get the coding and non coding exons
    coding_exons_bed = "{0}/coding_exons.bed".format(output_directory)
    non_coding_exons_bed = "{0}/non_coding_exons.bed".format(output_directory)
    coding_exons_fasta = "{0}/coding_exons.fasta".format(output_directory)
    non_coding_exons_fasta = "{0}/non_coding_exons.fasta".format(output_directory)
    # get the coding and non coding exons
    if clean_run or not os.path.isfile(non_coding_exons_fasta) or not os.path.isfile(coding_exons_fasta):
        print("Getting the coding and non coding exons...")
        get_non_coding_exons(genome_fasta, gtf, coding_exons_bed, non_coding_exons_bed, coding_exons_fasta, non_coding_exons_fasta, sequence_output_directory)

    # run the simulations
    print("Simulating the stop counts in non coding exons...")
    simoc.sim_stop_count(non_coding_exons_fasta, required_simulations, sim_stop_count_non_coding_output_file)
    print("Simulating the stop counts in coding exons...")
    simoc.sim_stop_count(coding_exons_fasta, required_simulations, sim_stop_count_coding_output_file)


def stop_density_test(gtf_file, genome_fasta, seqs_fasta, simulations, output_directory, output_file):

    genome_seq_outputs = "{0}/genome_sequence_files".format(output_directory)
    gen.create_output_directories(genome_seq_outputs)

    # get the sequences for non features
    features_bed = "{0}/genome_features.bed".format(genome_seq_outputs)
    non_features_bed = "{0}/non_genome_features.bed".format(genome_seq_outputs)
    non_features_fasta = "{0}/non_genome_features.fasta".format(genome_seq_outputs)
    seqo.get_non_transcribed_regions(gtf_file, genome_fasta, features_bed, non_features_bed, non_features_fasta, genome_seq_outputs)

    threshold = 0.5
    simoc.sim_stop_density(seqs_fasta, non_features_fasta, threshold, simulations, output_directory, output_file)


def cds_motif_test(cds_fasta, output_file):

    nts = ["A", "C", "G", "T"]
    stops = ["TAA", "TAG", "TGA"]

    codon_list = sorted(["".join(codon) for codon in it.product(nts, nts, nts)])
    combinations = [sorted(i) for i in it.combinations(codon_list, 3)]
    # combination_not_all_stops = [i for i in combinations if len(list(set(i) & set(stops))) < 3]

    combinations = combinations[:30]

    temp_dir = "temp_motif_densities"
    gen.create_output_directories(temp_dir)

    args = [cds_fasta, temp_dir]
    outputs = simoc.run_simulation_function(combinations, args, ops.calc_motif_densities, sim_run=False)

    temp_filelist = []
    for output in outputs:
        temp_filelist.append(output)

    densities = collections.defaultdict(lambda: collections.defaultdict())

    for file in temp_filelist:
        motif = file.split("/")[-1].split(".")[0]
        data = gen.read_many_fields(file, ",")[0]
        gc = data[0]
        density = data[1]
        densities[gc][motif] = density

    with open(output_file, "w") as outfile:
        outfile.write("motif,gc,density\n")
        for gc in sorted(densities):
            for motif in sorted(densities[gc]):
                outfile.write("{0},{1},{2}\n".format(motif, gc, densities[gc][motif]))
