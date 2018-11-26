import generic as gen
import sequence_ops as sequo

# TODO: simulation to get the normalised densities of lincRNA stop codons
# TODO: compare stop codon density in exons of lincRNA

def stop_density_test(input_fasta, input_families):
    exon_regions = sequo.get_exon_region_seqs(input_fasta)
