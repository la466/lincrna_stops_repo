import generic as gen
import sim_ops as simo
import file_ops as fo
import time

def main():

    description = ""
    args = gen.parse_arguments(description, ["source_exons_path", "genome_fasta", "output_directory", "required_simulations", "extract_seqs", "sim_orf_length"], flags = [4,5], ints=[3])
    source_exons_path, genome_fasta, output_directory, required_simulations, extract_seqs, sim_orf_length = args.source_exons_path, args.genome_fasta, args.output_directory, args.required_simulations, args.extract_seqs, args.sim_orf_length

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
        fo.extract_seqs(source_exons_path, genome_fasta, exons_bed, exons_fasta, seqs_fasta)

    # *****
    # might need to do some filtering of sequences in here
    # *****

    sim_orf_length_output_file = "{0}/sim_orf_lengths.csv".format(output_directory)
    if sim_orf_length:
        simo.sim_orf_length(seqs_fasta, required_simulations, sim_orf_length_output_file)



if __name__ == "__main__":
    main()
