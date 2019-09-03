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

    arguments = ["input_directory", "output_directory", "simulations", "motif_file", "input_fasta", "input_fasta2", "families_file", "output_prefix", "controls_dir", "output_prefix", "clean_run",  "generate_gc_matched_stop_sets", "generate_motif_dinucleotide_controls", "compare_stop_density", "compare_codon_density", "coding_exons", "generate_gc_controls_exons", "generate_gc_controls_introns", "generate_dint_exon_controls", "generate_dint_intron_controls", "cds_density_nd", "stop_density_nd", "without_ese", "exon_region_density", "intron_density", "calc_purine_content", "exon_region_excess", "non_coding_exons", "match_density", "match_subs", "intron_length_test", "seq_hits", "seq_hits_linc", "overlap", "motif_overlap", "motif_overlap_density", "overlap_diffs", "intron_hexamers", "ese_purine", "exon_intron_purine", "exon_intron_purine_no_motifs"]
    description = ""
    args = gen.parse_arguments(description, arguments, flags = [9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40], ints=[], opt_flags=[2,3,4,5,6,7,8])
    input_directory, \
    output_directory, \
    simulations, \
    motif_file, \
    input_fasta, \
    input_fasta2, \
    families_file, \
    output_prefix, \
    controls_dir, \
    output_prefix, \
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
    cds_density_nd, \
    stop_density_nd, \
    without_ese, \
    exon_region_density, \
    intron_density, \
    calc_purine_content, \
    exon_region_excess, \
    non_coding_exons, \
    match_density, \
    match_subs, \
    intron_length_test, \
    seq_hits, \
    seq_hits_linc, \
    overlap, \
    motif_overlap, \
    motif_overlap_density, \
    overlap_diffs, \
    intron_hexamers, \
    ese_purine, \
    exon_intron_purine, \
    exon_intron_purine_no_motifs = \
    args.input_directory, \
    args.output_directory, \
    args.simulations, \
    args.motif_file, \
    args.input_fasta, \
    args.input_fasta2, \
    args.families_file, \
    args.output_prefix, \
    args.controls_dir, \
    args.output_prefix, \
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
    args.cds_density_nd, \
    args.stop_density_nd, \
    args.without_ese, \
    args.exon_region_density, \
    args.intron_density, \
    args.calc_purine_content, \
    args.exon_region_excess, \
    args.non_coding_exons, \
    args.match_density, \
    args.match_subs, \
    args.intron_length_test, \
    args.seq_hits, \
    args.seq_hits_linc, \
    args.overlap, \
    args.motif_overlap, \
    args.motif_overlap_density, \
    args.overlap_diffs, \
    args.intron_hexamers, \
    args.ese_purine, \
    args.exon_intron_purine, \
    args.exon_intron_purine_no_motifs

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
    if not families_file:
        families_file = "{0}/genome_sequences/human/human.cds.families.bed".format(input_directory)
    coding_exons_fasta = "{0}/genome_sequences/human/human.cds.clean_coding_exons.fasta".format(input_directory)
    non_coding_exons_fasta = "{0}/genome_sequences/human/human.cds.clean_non_coding_exons.fasta".format(input_directory)

    # get the set of codons with the same gc content as stop codons
    gc_matched_stops_sets_file = "{0}/stops_gc_matched_motifs.bed".format(output_directory)
    if generate_gc_matched_stop_sets:
        seqo.get_gc_matched_motifs(stops, gc_matched_stops_sets_file)

    # if wanting to generate control motifs
    if generate_motif_dinucleotide_controls:
        simopc.generate_motif_dinucleotide_controls(motif_file, simulations, output_directory, match_density = match_density, match_subs = match_subs)

    # get hits to motifs
    if seq_hits:
        local_output_dir = "{0}/tests/ese_hits".format(output_directory)
        if output_prefix:
            sim_output_dir = "{0}/{1}_{2}".format(local_output_dir, output_prefix, motif_file.split("/")[-1].split(".")[0])
            final_output_file = "{0}/{1}_{2}_processed1.csv".format(local_output_dir, output_prefix, motif_file.split("/")[-1].split(".")[0])
        else:
            sim_output_dir = "{0}/{1}".format(local_output_dir, motif_file.split("/")[-1].split(".")[0])
            final_output_file = "{0}/{1}_processed1.csv".format(local_output_dir, motif_file.split("/")[-1].split(".")[0])
        gen.create_output_directories(sim_output_dir)

        runs = 10
        for run in range(runs):
            if output_prefix:
                output_file = "{0}/{1}_{2}_hits_{3}.csv".format(sim_output_dir, output_prefix, motif_file.split("/")[-1].split(".")[0], run+1)
            else:
                output_file = "{0}/{1}_hits_{2}.csv".format(sim_output_dir, motif_file.split("/")[-1].split(".")[0], run+1)
            mto.calc_seq_hits(input_fasta, input_fasta2, output_file, motif_file, controls_dir, required_simulations = simulations, families_file = families_file)
        mto.process_seq_hits(sim_output_dir, final_output_file)


    if overlap:
        local_output_dir = "{0}/tests/motif_overlap_capability".format(output_directory)
        gen.create_output_directories(local_output_dir)
        output_file1 = "{0}/{1}_motif_overlap_capability.csv".format(local_output_dir, motif_file.split("/")[-1].split(".")[0])
        mto.calc_overlap_potential(input_fasta, motif_file, output_file1, required_simulations = simulations, controls_directory = motif_simulations_directory, families_file = families_file)
    if overlap_diffs:
        local_output_dir = "{0}/tests/motif_overlap_capability".format(output_directory)
        gen.create_output_directories(local_output_dir)
        output_file1 = "{0}/{1}_motif_overlap_diff.csv".format(local_output_dir, motif_file.split("/")[-1].split(".")[0])
        mto.calc_overlap_diffs(input_fasta, motif_file, output_file1, required_simulations = simulations, controls_directory = motif_simulations_directory, families_file = families_file)

    # need to move to main
    if motif_overlap:
        local_output_dir = "{0}/tests/motif_overlaps".format(output_directory)
        gen.create_output_directories(local_output_dir)
        gen.check_files_exists([input_fasta, motif_file, families_file])
        output_file = "{0}/{1}_{2}_motif_overlap_chisquare.csv".format(local_output_dir, output_prefix, motif_file.split("/")[-1].split(".")[0])
        runs = 10
        ltests.motif_overlap_test(input_fasta, motif_file, output_file, runs = runs, families_file = families_file)

    if motif_overlap_density:
        local_output_dir = "{0}/tests/motif_overlaps".format(output_directory)
        gen.create_output_directories(local_output_dir)
        gen.check_files_exists([input_fasta, motif_file, families_file])
        output_file = "{0}/{1}_{2}_motif_overlap_density.csv".format(local_output_dir, output_prefix, motif_file.split("/")[-1].split(".")[0])
        runs = 10
        ltests.motif_overlap_density_test(input_fasta, motif_file, output_file, runs = runs, families_file = families_file)


    # can use tihs for lincs too
    if intron_length_test:
        local_output_dir = "{0}/tests/intron_length_ese_densities".format(output_directory)
        gen.create_output_directories(local_output_dir)
        output_file = "{0}/{1}_{2}_ese_densities_all_seq.csv".format(local_output_dir, motif_file.split("/")[-1].split(".")[0], output_prefix)
        output_file_all_genes = "{0}/{1}_{2}_ese_densities_all_seq_all_genes.csv".format(local_output_dir, motif_file.split("/")[-1].split(".")[0], output_prefix)
        output_file_matched_size = "{0}/{1}_{2}_ese_densities.csv".format(local_output_dir, motif_file.split("/")[-1].split(".")[0], output_prefix)
        output_file_flanks = "{0}/{1}_{2}_ese_densities_flanks.csv".format(local_output_dir, motif_file.split("/")[-1].split(".")[0], output_prefix)
        # # run on whole sequence
        mto.intron_length_test(input_fasta, input_fasta2, motif_file, output_file, flanks = None, families_file = families_file)
        # dont group by familiy
        mto.intron_length_test(input_fasta, input_fasta2, motif_file, output_file_all_genes, flanks = None, families_file = None)
        # run on just longer ones
        mto.intron_length_test(input_fasta, input_fasta2, motif_file, output_file_matched_size, flanks = None, families_file = families_file, restrict_size = True)
        # # run on flanks
        mto.intron_length_test(input_fasta, input_fasta2, motif_file, output_file_flanks, flanks = True, families_file = families_file)




    # other stuff

    # generate controls from non translated sequenece
    gc_control_exon_output_directory = "{0}/clean_exon_gc_controls".format(output_directory)
    if generate_gc_controls_exons:
        simopc.generate_gc_controls(exons_fasta, non_transcribed_fasta, gc_control_exon_output_directory)

    # generate intron controls frmo unstranslated sequence
    gc_control_intron_output_directory = "{0}/clean_intron_gc_controls".format(output_directory)
    if generate_gc_controls_introns:
        simopc.generate_gc_controls(introns_fasta, non_transcribed_fasta, gc_control_intron_output_directory)

    # dinucleotide matched controls from exons
    dint_control_cds_output_directory = "{0}/clean_cds_dint_controls".format(output_directory)
    if generate_dint_exon_controls:
        gen.create_output_directories(dint_control_cds_output_directory)
        simopc.generate_dint_controls(cds_fasta, dint_control_cds_output_directory)

    # dinuceotide matched controls from introns
    dint_control_intron_output_directory = "{0}/clean_intron_dint_controls".format(output_directory)
    if generate_dint_intron_controls:
        gen.create_output_directories(dint_control_intron_output_directory)
        simopc.generate_dint_intron_controls(introns_fasta, dint_control_intron_output_directory)

    gen.create_output_directories("{0}/tests/intron_density".format(output_directory))
    intron_density_sim_file = "{0}/intron_density/motif_sets_intron_density.csv".format(output_directory)
    if intron_density:
        if not motif_file:
            print("Please provide an ESE file...")
            raise Exception
        intron_density_controls_dir = "clean_run/motif_controls/{0}".format(motif_file.split("/")[-1].split(".")[0])
        gen.create_output_directories(intron_density_controls_dir)
        mto.calculate_intron_densities(motif_file, introns_fasta, intron_density_sim_file, intron_density_controls_dir, families_file = families_file, required_simulations = simulations, combined = True, clean_run = clean_run)

    if calc_purine_content:
        output_file = "{0}/purine_content.csv".format(output_directory)
        mto.calc_purine_content(coding_exons_fasta, introns_fasta, output_file, families_file = families_file)

    if intron_hexamers:
        local_output_directory = "{0}/tests/introns".format(output_directory)
        output_file = "{0}/intron_hexamers.csv".format(local_output_directory)
        mto.intron_hexamer_test(input_fasta, motif_file, local_output_directory, output_file, required_simulations = simulations, families_file = families_file)

    if ese_purine:
        local_output_directory = "{0}/tests/purine_content".format(output_directory)
        output_file = "{0}/{1}_purine_content.csv".format(local_output_directory, motif_file.split("/")[-1].split(".")[0])
        mto.motif_purine_content(input_fasta, motif_file, output_file, families_file = families_file)

    if exon_intron_purine:
        local_output_directory = "{0}/tests/purine_content".format(output_directory)
        output_file = "{0}/exon_intron_purine.csv".format(local_output_directory)
        mto.exon_intron_purine(input_fasta, input_fasta2, output_file, families_file = families_file)

    if exon_intron_purine_no_motifs:
        local_output_directory = "{0}/tests/purine_content".format(output_directory)
        output_file = "{0}/exon_intron_purine_no_eses.csv".format(local_output_directory)
        mto.exon_intron_purine_no_motifs(input_fasta, input_fasta2, motif_file, output_file, families_file = families_file)


if __name__ == "__main__":
    main()
