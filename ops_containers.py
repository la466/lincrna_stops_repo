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


def cds_motif_test(cds_fasta, output_file):

    nts = ["A", "C", "G", "T"]
    stops = ["TAA", "TAG", "TGA"]

    codon_list = sorted(["".join(codon) for codon in it.product(nts, nts, nts)])
    combinations = [sorted(i) for i in it.combinations(codon_list, 3)]
    # combination_not_all_stops = [i for i in combinations if len(list(set(i) & set(stops))) < 3]

    # combinations = combinations[:30]

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

    iterator = 0
    with open(output_file, "w") as outfile:
        outfile.write("id,motif,gc,density\n")
        for gc in sorted(densities):
            for motif in sorted(densities[gc]):
                iterator += 1
                outfile.write("{0},{1},{2},{3}\n".format(iterator,motif, gc, densities[gc][motif]))



def motif_codon_density(motif_file, output_directory):

    stops = ["TAA", "TAG", "TGA"]
    gc_matchd_motifs_file = "{0}/gc_matched_combinations.bed".format(output_directory)
    if not os.path.isfile(gc_matchd_motifs_file):
        seqo.get_gc_matched_motifs(stops, gc_matchd_motifs_file)

    temp_dir = "temp_motif_density"
    gen.create_output_directories(temp_dir)

    motif_sets = gen.read_many_fields(gc_matchd_motifs_file, "\t")
    motif_sets.append(["TAA", "TAG", "TGA"])

    args = [motif_file, temp_dir]
    outputs = simoc.run_simulation_function(motif_sets, args, ops.calc_codon_density_in_motifs, sim_run=False)

    new_output_dir = "{0}/motif_densities".format(output_directory)
    gen.create_output_directories(new_output_dir)


    output_file = "{0}/{1}.csv".format(new_output_dir, motif_file.split("/")[-1].split(".")[0])
    with open(output_file, "w") as outfile:
        outfile.write("id,motifs,density\n")
        for i,file in enumerate(sorted(outputs)):
            data = gen.read_many_fields(file, ",")[0]
            outfile.write("{0},{1},{2}\n".format(i+1,data[0],data[1]))

    gen.remove_directory(temp_dir)
