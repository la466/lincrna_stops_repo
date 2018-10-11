import generic as gen
import sim_ops_containers as simopc
import ops_containers as opsc
import file_ops as fo
import time

def main():

    description = ""
    args = gen.parse_arguments(description, ["source_exons_path", "genome_fasta", "gtf", "output_directory", "required_simulations", "extract_seqs", "sim_orf_length", "sim_stop_count", "non_coding_exons"], flags = [5,6,7,8], ints=[4])
    source_exons_path, genome_fasta, gtf, output_directory, required_simulations, extract_seqs, sim_orf_length, sim_stop_count, non_coding_exons = args.source_exons_path, args.genome_fasta, args.gtf, args.output_directory, args.required_simulations, args.extract_seqs, args.sim_orf_length, args.sim_stop_count, args.non_coding_exons

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
        fo.extract_seqs(source_exons_path, genome_fasta, exons_bed, exons_fasta, seqs_fasta, exclude_XY=True)

    # *****
    # might need to do some filtering of sequences in here
    # *****

    sim_orf_length_output_file = "{0}/sim_orf_lengths.csv".format(output_directory)
    if sim_orf_length:
        simopc.sim_orf_length(seqs_fasta, required_simulations, sim_orf_length_output_file)

    sim_stop_count_output_file = "{0}/sim_stop_count.csv".format(output_directory)
    if sim_stop_count:
        simopc.sim_stop_count(seqs_fasta, required_simulations, sim_stop_count_output_file)

    if non_coding_exons:
        opsc.non_coding_exons(genome_fasta, gtf, output_directory)



if __name__ == "__main__":
    main()
