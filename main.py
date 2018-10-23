import generic as gen
import sim_ops_containers as simopc
import ops_containers as opsc
import file_ops as fo
import time

def main():

    description = "Main file to ask questions about the prevelance of stop codons in sequences."
    args = gen.parse_arguments(description, ["source_exons_path", "genome_fasta", "gtf", "mapping_file", "expression_file", "motif_file", "output_directory", "required_simulations", "extract_seqs", "sim_orf_length", "sim_stop_count", "lincRNA_expression", "non_coding_exons", "stop_counts_regions", "stop_counts_exon_ends", "sim_motifs", "stop_density_test", "clean_run"], flags = [8,9,10,11,12,13,14,15,16,17], ints=[7])
    source_exons_path, genome_fasta, gtf, mapping_file, expression_file, motif_file, output_directory, required_simulations, extract_seqs, sim_orf_length, sim_stop_count, lincRNA_expression, non_coding_exons, stop_counts_regions, stop_counts_exon_ends, sim_motifs, stop_density_test, clean_run = args.source_exons_path, args.genome_fasta, args.gtf, args.mapping_file, args.expression_file, args.motif_file, args.output_directory, args.required_simulations, args.extract_seqs, args.sim_orf_length, args.sim_stop_count, args.lincRNA_expression, args.non_coding_exons, args.stop_counts_regions, args.stop_counts_exon_ends, args.sim_motifs, args.stop_density_test, args.clean_run

    # set a start time
    start = time.time()

    # create the output_directory if it doenst already exist
    gen.create_output_directories(output_directory)

    # set up the output fasta to contain the exon seqs
    exons_bed = "{0}/lincRNA_exons.bed".format(output_directory)
    exons_fasta = "{0}/lincRNA_exons.fasta".format(output_directory)
    seqs_fasta = "{0}/lincRNA_seqs.fasta".format(output_directory)

    # extract sequences from source file
    if extract_seqs:
        print("Extracting lincRNA seqs...")
        fo.extract_seqs(source_exons_path, genome_fasta, exons_bed, exons_fasta, seqs_fasta, mapping_file, exclude_XY=True)

    # *****
    # might need to do some filtering of sequences in here
    # *****

    sim_orf_length_output_file = "{0}/sim_orf_lengths.csv".format(output_directory)
    if sim_orf_length:
        simopc.sim_orf_length(seqs_fasta, required_simulations, sim_orf_length_output_file)

    sim_stop_count_output_file = "{0}/sim_stop_count.csv".format(output_directory)
    if sim_stop_count:
        simopc.sim_stop_count(seqs_fasta, required_simulations, sim_stop_count_output_file)

    # sim_stop_count_mm_output_file = "{0}/sim_stop_count_mm.csv".format(output_directory)
    # if sim_stop_count_mm:
    #     simopc.sim_stop_count_mm(seqs_fasta, required_simulations, sim_stop_count_mm_output_file)

    lincRNA_expression_output_file = "{0}/lincRNA_expression_links.csv".format(output_directory)
    if lincRNA_expression:
        opsc.lincRNA_expression(seqs_fasta, mapping_file, expression_file, lincRNA_expression_output_file)

    sim_stop_count_non_coding_exons_file = "{0}/sim_stop_count_non_coding_exons.csv".format(output_directory)
    sim_stop_count_coding_exons_file = "{0}/sim_stop_count_coding_exons.csv".format(output_directory)
    if non_coding_exons:
        opsc.non_coding_exons(genome_fasta, gtf, output_directory, sim_stop_count_non_coding_exons_file, sim_stop_count_coding_exons_file, required_simulations, clean_run)

    sto_count_coding_exons_file = "{0}/sim_coding_exons_stop_count_regions.csv".format(output_directory)
    if stop_counts_regions:
        simopc.sim_coding_exon_stop_counts_regions(genome_fasta, gtf, output_directory, sto_count_coding_exons_file, required_simulations, clean_run=None)

    stop_count_coding_exons_ends_file = "{0}/sim_coding_exons_stop_count_exon_ends.csv".format(output_directory)
    if stop_counts_exon_ends:
        simopc.sim_coding_exon_flanks_stop_counts(genome_fasta, gtf, output_directory, stop_count_coding_exons_ends_file, required_simulations, clean_run=None)

    motifs_sim_file = "{0}/sim_rbp_motifs_stop_count.csv".format(output_directory)
    if sim_motifs:
        simopc.sim_motifs(motif_file, motifs_sim_file, required_simulations)

    if stop_density_test:
        opsc.stop_density_test(gtf, genome_fasta, output_directory)



if __name__ == "__main__":
    main()
