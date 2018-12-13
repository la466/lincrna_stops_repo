import generic as gen
import containers as cont
import sim_ops_containers as simopc
import main_tests_ops as mto
import ops_containers as opsc
import seq_ops as seqo
import file_ops as fo
import time
from useful_motif_sets import stops

def main():

    arguments = ["input_directory", "output_directory", "generate_gc_matched_stop_sets", "generate_motif_dinucleotide_controls", "compare_stop_density", "compare_codon_density", "coding_exons", "generate_gc_controls_exons", "generate_gc_controls_introns", "generate_dint_exon_controls", "generate_dint_intron_controls", "cds_ds", "cds_codon_ds", "cds_density_nd", "stop_density_nd", "without_ese", "exon_region_density"]

    description = ""
    args = gen.parse_arguments(description, arguments, flags = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16], ints=[])
    input_directory, output_directory, generate_gc_matched_stop_sets, generate_motif_dinucleotide_controls, compare_stop_density, compare_codon_density, coding_exons, generate_gc_controls_exons, generate_gc_controls_introns, generate_dint_exon_controls, generate_dint_intron_controls, cds_ds, cds_codon_ds, cds_density_nd, stop_density_nd, without_ese, exon_region_density = args.input_directory, args.output_directory, args.generate_gc_matched_stop_sets, args.generate_motif_dinucleotide_controls, args.compare_stop_density, args.compare_codon_density, args.coding_exons, args.generate_gc_controls_exons, args.generate_gc_controls_introns, args.generate_dint_exon_controls, args.generate_dint_intron_controls, args.cds_ds, args.cds_codon_ds, args.cds_density_nd, args.stop_density_nd, args.without_ese, args.exon_region_density

    # set a start time
    start = time.time()

    # create the output_directory if it doenst already exist
    gen.create_output_directories(output_directory)

    exons_fasta = "{0}/genome_sequences/human/human.cds.clean_coding_exons.fasta".format(input_directory)
    cds_fasta = "{0}/genome_sequences/human/human.cds.clean.fasta".format(input_directory)
    introns_fasta = "{0}/genome_sequences/human/human.clean_introns.fasta".format(input_directory)
    non_transcribed_fasta = "{0}/genome_sequences/human/human.non_transcribed.fasta".format(input_directory)
    families_file = "{0}/genome_sequences/human/human.cds.families.bed".format(input_directory)
    coding_exons_fasta = "{0}/genome_sequences/human/human.cds.clean_coding_exons.fasta".format(input_directory)

    ese_file = "source_data/motif_sets/int3.txt"

    # get the set of codons with the same gc content as stop codons
    gc_matched_stops_sets_file = "{0}/stops_gc_matched_motifs.bed".format(output_directory)
    if generate_gc_matched_stop_sets:
        seqo.get_gc_matched_motifs(stops, gc_matched_stops_sets_file)

    motif_simulations_directory = "{0}/{1}_dinucleotide_controls".format(output_directory, ese_file.split("/")[-1].split(".")[0])
    if generate_motif_dinucleotide_controls:
        simopc.generate_motif_dinucleotide_controls(ese_file, 1000, motif_simulations_directory)

    gc_control_exon_output_directory = "{0}/clean_exon_gc_controls".format(output_directory)
    if generate_gc_controls_exons:
        simopc.generate_gc_controls(exons_fasta, non_transcribed_fasta, gc_control_exon_output_directory)

    gc_control_intron_output_directory = "{0}/clean_intron_gc_controls".format(output_directory)
    if generate_gc_controls_introns:
        simopc.generate_gc_controls(introns_fasta, non_transcribed_fasta, gc_control_intron_output_directory)

    dint_control_cds_output_directory = "{0}/clean_cds_dint_controls".format(output_directory)
    if generate_dint_exon_controls:
        simopc.generate_dint_controls(cds_fasta, dint_control_cds_output_directory)

    dint_control_intron_output_directory = "{0}/clean_intron_dint_controls".format(output_directory)
    if generate_dint_intron_controls:
        simopc.generate_dint_intron_controls(introns_fasta, dint_control_intron_output_directory)

    # get the stop density
    output_file = "{0}/_test_compare_exon_intron_stop_density.csv".format(output_directory)
    if compare_stop_density:
        mto.compare_stop_density(exons_fasta, introns_fasta, output_file, families_file = families_file)

    # get the stop density
    output_codon_density_directory = "{0}/tests/compare_exon_intron_codon_densities".format(output_directory)
    if compare_codon_density:
        mto.compare_codon_density(exons_fasta, introns_fasta, output_codon_density_directory, families_file = families_file)

    output_file = "{0}/exonic_stop_density_nd.csv".format(output_directory)
    if stop_density_nd:
        mto.stop_density_nd(exons_fasta, cds_fasta, introns_fasta, dint_control_cds_output_directory, output_file, families_file = families_file)


    if without_ese:
        mto.compare_density_no_ese(exons_fasta, cds_fasta, ese_file, families_file = families_file)

    # if coding_exons:
    #     mto.coding_exons(input_file1, input_file2, output_directory)

    cds_fasta = "{0}/genome_sequences/human/human.cds.clean.fasta".format(input_directory)
    ortholog_cds_fasta = "{0}/genome_sequences/macaque/macaque.cds.quality_filtered.step1.fasta".format(input_directory)
    ortholog_transcript_links = "{0}/genome_sequences/human.macaque.conservation_filtering.step3.bed".format(input_directory)
    alignments_file = "{0}/genome_sequences/human.macaque.alignments.fasta".format(input_directory)
    ese_file = "source_data/motif_sets/int3.txt"
    ds_output_directory = "{0}/tests/ese_ds".format(output_directory)
    output_file = "{0}/codon_ds.csv".format(ds_output_directory)
    output_file1 = "{0}/codons_ds.csv".format(ds_output_directory)
    if cds_ds:
        mto.calc_ds(alignments_file, cds_fasta, ortholog_cds_fasta, ortholog_transcript_links, ese_file, ds_output_directory, output_file, motif_simulations_directory = motif_simulations_directory, families_file = families_file)

    if cds_codon_ds:
        mto.calc_codon_ds(alignments_file, cds_fasta, ortholog_cds_fasta, ortholog_transcript_links, ese_file, ds_output_directory, output_file1, families_file = families_file, codon_sets_file = gc_matched_stops_sets_file)

    if exon_region_density:
        mto.exon_region_density(cds_fasta, coding_exons_fasta, gc_matched_stops_sets_file, families_file = families_file)


if __name__ == "__main__":
    main()
