import generic as gen
import containers as cont
import sim_ops_containers as simopc
import main_tests_ops as mto
import ops_containers as opsc
import seq_ops as seqo
import sequence_ops as sequo
import file_ops as fo
import time
from useful_motif_sets import stops
import os


def main():


    arguments = ["input_directory", "output_directory", "simulations", "ese_file", "clean_run",  "generate_gc_matched_stop_sets", "generate_motif_dinucleotide_controls", "compare_stop_density", "compare_codon_density", "coding_exons", "generate_gc_controls_exons", "generate_gc_controls_introns", "generate_dint_exon_controls", "generate_dint_intron_controls", "cds_ds", "cds_ds_all", "cds_codon_ds", "cds_density_nd", "stop_density_nd", "without_ese", "exon_region_density", "intron_density", "calc_purine_content", "exon_region_excess", "non_coding_exons", "cds_ds_mutation", "ese_ds", "ese_ds_randomise_overlap", "ese_ds_mutation", "non_ese_stop_ds", "ese_ds_test", "match_density"]

    description = ""
    args = gen.parse_arguments(description, arguments, flags = [4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31], ints=[], opt_flags=[2,3])
    input_directory, \
    output_directory, \
    simulations, \
    ese_file, \
    clean_run, \
    generate_gc_matched_stop_sets, \
    generate_motif_dinucleotide_controls, \
    compare_stop_density, \
    compare_codon_density, \
    coding_exons, \
    generate_gc_controls_exons, \
    generate_gc_controls_introns, \
    generate_dint_exon_controls, \
    generate_dint_intron_controls, \
    cds_ds, \
    cds_ds_all, \
    cds_codon_ds, \
    cds_density_nd, \
    stop_density_nd, \
    without_ese, \
    exon_region_density, \
    intron_density, \
    calc_purine_content, \
    exon_region_excess, \
    non_coding_exons, \
    cds_ds_mutation, \
    ese_ds, \
    ese_ds_randomise_overlap, \
    ese_ds_mutation, \
    non_ese_stop_ds, \
    ese_ds_test, \
    match_density = \
    args.input_directory, \
    args.output_directory, \
    args.simulations, \
    args.ese_file, \
    args.clean_run, \
    args.generate_gc_matched_stop_sets, \
    args.generate_motif_dinucleotide_controls, \
    args.compare_stop_density, \
    args.compare_codon_density, \
    args.coding_exons, \
    args.generate_gc_controls_exons, \
    args.generate_gc_controls_introns, \
    args.generate_dint_exon_controls, \
    args.generate_dint_intron_controls, \
    args.cds_ds, \
    args.cds_ds_all, \
    args.cds_codon_ds, \
    args.cds_density_nd, \
    args.stop_density_nd, \
    args.without_ese, \
    args.exon_region_density, \
    args.intron_density, \
    args.calc_purine_content, \
    args.exon_region_excess, \
    args.non_coding_exons, \
    args.cds_ds_mutation, \
    args.ese_ds, \
    args.ese_ds_randomise_overlap, \
    args.ese_ds_mutation, \
    args.non_ese_stop_ds, \
    args.ese_ds_test, \
    args.match_density

    if simulations:
        simulations = int(simulations)

    # set a start time
    start = time.time()

    # create the output_directory if it doenst already exist
    gen.create_output_directories(output_directory)

    exons_fasta = "{0}/genome_sequences/human/human.cds.clean_coding_exons.fasta".format(input_directory)
    cds_fasta = "{0}/genome_sequences/human/human.cds.clean.fasta".format(input_directory)
    multi_exon_cds_fasta = "{0}/genome_sequences/human/human.cds.multi_exons.fasta".format(input_directory)
    introns_fasta = "{0}/genome_sequences/human/human.clean_introns.fasta".format(input_directory)
    non_transcribed_fasta = "{0}/genome_sequences/human/human.non_transcribed.fasta".format(input_directory)
    families_file = "{0}/genome_sequences/human/human.cds.families.bed".format(input_directory)
    coding_exons_fasta = "{0}/genome_sequences/human/human.cds.clean_coding_exons.fasta".format(input_directory)
    non_coding_exons_fasta = "{0}/genome_sequences/human/human.cds.clean_non_coding_exons.fasta".format(input_directory)

    # ese_file = "source_data/motif_sets/int3.txt"

    # get the set of codons with the same gc content as stop codons
    gc_matched_stops_sets_file = "{0}/stops_gc_matched_motifs.bed".format(output_directory)
    if generate_gc_matched_stop_sets:
        seqo.get_gc_matched_motifs(stops, gc_matched_stops_sets_file)

    motif_simulations_directory = "{0}/dinucleotide_controls/{1}_dinucleotide_controls".format(output_directory, ese_file.split("/")[-1].split(".")[0])
    if match_density:
        motif_simulations_directory += "_matched_stops"

    if generate_motif_dinucleotide_controls:
        print(motif_simulations_directory)
        simopc.generate_motif_dinucleotide_controls(ese_file, simulations, motif_simulations_directory, match_density = match_density)

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
    output_file = "{0}/tests/compare_exon_intron_stop_density.csv".format(output_directory)
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
    # ese_file = "source_data/motif_sets/int3.txt"
    ds_output_directory = "{0}/tests/ese_ds".format(output_directory)
    output_file = "{0}/codon_ds_matched_densities.csv".format(ds_output_directory)
    output_file1 = "{0}/codons_ds.csv".format(ds_output_directory)
    output_file2 = "{0}/codon_ds_stats.csv".format(ds_output_directory)
    output_file3 = "{0}/codon_ds_stats_mutation.csv".format(ds_output_directory)
    if cds_ds:
        mto.calc_ds(alignments_file, cds_fasta, ortholog_cds_fasta, ortholog_transcript_links, ese_file, ds_output_directory, output_file, motif_controls_directory = motif_simulations_directory, families_file = families_file)

    if cds_ds_all:
        mto.calc_ds_all(simulations, alignments_file, cds_fasta, multi_exon_cds_fasta,  ortholog_cds_fasta, ortholog_transcript_links, ese_file, ds_output_directory, output_file2, motif_controls_directory = motif_simulations_directory, families_file = families_file)

    if cds_codon_ds:
        mto.calc_codon_ds(alignments_file, cds_fasta, ortholog_cds_fasta, ortholog_transcript_links, ese_file, ds_output_directory, output_file1, families_file = families_file, codon_sets_file = gc_matched_stops_sets_file)

    if cds_ds_mutation:
        mto.calc_ds_mutation(simulations, alignments_file, cds_fasta, ortholog_cds_fasta, ortholog_transcript_links, ese_file, ds_output_directory, output_file3, motif_controls_directory = motif_simulations_directory, families_file = families_file)


    ese_ds_test_output = "{0}/ese_ds_test.csv".format(ds_output_directory)
    alignment_file = "{0}/genome_sequences/human.macaque.alignment.fasta".format(output_directory)
    if ese_ds_test:
        # if the alignments file doesn't exist, create it
        if not os.path.isfile(alignment_file):
            sequo.retrieve_alignments(cds_fasta, ortholog_cds_fasta, ortholog_transcript_links, alignment_file)
        # create the output directory
        gen.create_output_directories(ds_output_directory)
        # run the test
        mto.ese_ds(alignment_file, multi_exon_cds_fasta, ese_file, ese_ds_test_output, families_file = families_file)



    ese_ds_output_file = "{0}/ese_ds.csv".format(ds_output_directory)
    if ese_ds:
        # if the sequence alignments file doesnt exist, create it
        if not os.path.isfile(alignments_file):
            sequo.extract_alignments_from_file(cds_fasta, ortholog_cds_fasta, ortholog_transcript_links, alignments_file)
        # create the output directory for the test
        gen.create_output_directories(ds_output_directory)
        # run the test
        mto.calc_ese_ds(alignments_file, cds_fasta, multi_exon_cds_fasta, ese_file, ese_ds_output_file, simulations = simulations, controls_directory = motif_simulations_directory, families_file = families_file)

    ese_ds_randomise_overlap_output_file = "{0}/ese_ds_randomise_overlap.csv".format(ds_output_directory)
    ese_overlap_randomise_dir = "{0}/controls/{1}_overlap_sites_randomisation".format(output_directory, ese_file.split("/")[-1].split(".")[0])
    if ese_ds_randomise_overlap:
        # if the sequence alignments file doesnt exist, create it
        if not os.path.isfile(alignments_file):
            sequo.extract_alignments_from_file(cds_fasta, ortholog_cds_fasta, ortholog_transcript_links, alignments_file)
        if clean_run:
            gen.remove_directory(ese_overlap_randomise_dir)
        gen.create_output_directories(ese_overlap_randomise_dir)
        if len(os.listdir(ese_overlap_randomise_dir)) < simulations:
            required = simulations - len(os.listdir(ese_overlap_randomise_dir))
            sequo.randomise_ese_overlaps_wrapper(alignments_file, multi_exon_cds_fasta, ese_file, ese_overlap_randomise_dir, required = required)
        # create the output directory for the test
        gen.create_output_directories(ds_output_directory)
        # run the test
        mto.calc_ese_ds_random_overlaps(alignments_file, cds_fasta, multi_exon_cds_fasta, ese_file, ese_ds_randomise_overlap_output_file, simulations = simulations, controls_directory = ese_overlap_randomise_dir, families_file = families_file)

    ese_ds_mutation_output_file = "{0}/ese_ds_mutation1.csv".format(ds_output_directory)
    if ese_ds_mutation:
        # if the sequence alignments file doesnt exist, create it
        if not os.path.isfile(alignments_file):
            sequo.extract_alignments_from_file(cds_fasta, ortholog_cds_fasta, ortholog_transcript_links, alignments_file)
        # create the output directory for the test
        gen.create_output_directories(ds_output_directory)
        # run the test
        mto.calc_ese_ds_mutation(alignments_file, cds_fasta, multi_exon_cds_fasta, ese_file, ese_ds_mutation_output_file, simulations = simulations, controls_directory = motif_simulations_directory, families_file = families_file)

    non_ese_output_file = "{0}/non_ese_ds.csv".format(ds_output_directory)
    if non_ese_stop_ds:
        # if the sequence alignments file doesnt exist, create it
        if not os.path.isfile(alignments_file):
            sequo.extract_alignments_from_file(cds_fasta, ortholog_cds_fasta, ortholog_transcript_links, alignments_file)
        # create the output directory for the test
        gen.create_output_directories(ds_output_directory)
        # run the test
        mto.calc_non_ese_ds(alignments_file, cds_fasta, multi_exon_cds_fasta, ese_file, non_ese_output_file, simulations = simulations, controls_directory = motif_simulations_directory, families_file = families_file)


    exon_regions_directory = "{0}/tests/exon_regions".format(output_directory)
    gen.create_output_directories(exon_regions_directory)
    exon_region_density_file = "{0}/exon_regions_density.csv".format(exon_regions_directory)
    if exon_region_density:
        mto.exon_region_density(cds_fasta, coding_exons_fasta, gc_matched_stops_sets_file, exon_region_density_file, families_file = families_file)
    exon_region_excess_file = "{0}/exon_regions_excess.csv".format(exon_regions_directory)
    if exon_region_excess:
        mto.exon_region_excess(cds_fasta, coding_exons_fasta, exon_region_excess_file, families_file = families_file)

    non_coding_exons_dir = "{0}/tests/non_coding_exons".format(output_directory)
    non_coding_exons_file = "{0}/non_coding_exons.csv".format(non_coding_exons_dir)
    if non_coding_exons:
        mto.non_coding_exons(non_coding_exons_fasta, non_coding_exons_file, families_file = families_file)

    gen.create_output_directories("{0}/tests/intron_density".format(output_directory))
    intron_density_sim_file = "{0}/intron_density/motif_sets_intron_density.csv".format(output_directory)
    if intron_density:
        if not ese_file:
            print("Please provide an ESE file...")
            raise Exception
        intron_density_controls_dir = "clean_run/motif_controls/{0}".format(ese_file.split("/")[-1].split(".")[0])
        gen.create_output_directories(intron_density_controls_dir)
        mto.calculate_intron_densities(ese_file, introns_fasta, intron_density_sim_file, intron_density_controls_dir, families_file = families_file, required_simulations = simulations, combined = True, clean_run = clean_run)

    output_file = "{0}/purine_content.csv".format(output_directory)
    if calc_purine_content:
        mto.calc_purine_content(coding_exons_fasta, introns_fasta, output_file, families_file = families_file)


if __name__ == "__main__":
    main()
