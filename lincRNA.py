import generic as gen
import containers as cont
import lincRNA_tests_ops as ltests
import sim_ops_containers as simopc
import sequence_ops as sequo
import main_tests_ops as mto
import seq_ops as seqo
import time
import os

def main():

    arguments = ["input_bed", "input_fasta", "output_directory", "input_fasta2", "input_file", "required_simulations", "motif_file", "families_file", "output_prefix", "controls_dir", "extract_sequences", "calc_gc", "density_sim", "get_exon_dint_controls", "get_intron_dint_controls", "exon_region_density", "compare_stop_density", "sim_orf_lengths", "sim_orf_lengths_masked", "sim_stop_density", "sim_stop_density_introns", "sim_stop_density_within_genes", "sim_stop_density_removed_motifs", "sim_stop_density_removed_motifs_sim_seqs", "sim_stop_density_diff", "exon_intron_density", "motif_nd", "excess_test", "single_exon", "motif_overlap", "motif_overlap_density", "clean_alignments", "seq_hits_linc", "excess_length_thresholds", "upstream_atg"]
    description = "Container for analysis on lincRNAs"
    args = gen.parse_arguments(description, arguments, flags = [10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34], opt_flags = [3,4,5,6,7,8,9])
    input_bed, \
    input_fasta, \
    output_directory, \
    input_fasta2, \
    input_file, \
    required_simulations, \
    motif_file, \
    families_file, \
    output_prefix, \
    controls_dir, \
    extract_sequences, \
    calc_gc, \
    density_sim,  \
    get_exon_dint_controls, \
    get_intron_dint_controls, \
    exon_region_density, \
    compare_stop_density, \
    sim_orf_lengths, \
    sim_orf_lengths_masked, \
    sim_stop_density, \
    sim_stop_density_introns, \
    sim_stop_density_within_genes, \
    sim_stop_density_removed_motifs, \
    sim_stop_density_removed_motifs_sim_seqs, \
    sim_stop_density_diff, \
    exon_intron_density, \
    motif_nd, \
    excess_test, \
    single_exon,\
    motif_overlap, \
    motif_overlap_density, \
    clean_alignments, \
    seq_hits_linc, \
    excess_length_thresholds, \
    upstream_atg = \
    args.input_bed, \
    args.input_fasta, \
    args.output_directory, \
    args.input_fasta2, \
    args.input_file, \
    args.required_simulations, \
    args.motif_file, \
    args.families_file, \
    args.output_prefix, \
    args.controls_dir, \
    args.extract_sequences, \
    args.calc_gc, \
    args.density_sim, \
    args.get_exon_dint_controls, \
    args.get_intron_dint_controls, \
    args.exon_region_density, \
    args.compare_stop_density, \
    args.sim_orf_lengths, \
    args.sim_orf_lengths_masked, \
    args.sim_stop_density, \
    args.sim_stop_density_introns, \
    args.sim_stop_density_within_genes, \
    args.sim_stop_density_removed_motifs, \
    args.sim_stop_density_removed_motifs_sim_seqs, \
    args.sim_stop_density_diff, \
    args.exon_intron_density, \
    args.motif_nd, \
    args.excess_test, \
    args.single_exon, \
    args.motif_overlap, \
    args.motif_overlap_density, \
    args.clean_alignments, \
    args.seq_hits_linc, \
    args.excess_length_thresholds, \
    args.upstream_atg

    # make required simultions an int
    required_simulations = int(required_simulations) if required_simulations else None
    # prcoess output prefix
    output_prefix = output_prefix + "_" if output_prefix else ""

    # create output directories
    global_output_directory = "{0}/tests/lincrna".format(output_directory)
    gen.create_output_directories(global_output_directory)

    # set a start time
    start = time.time()

    # create the output_directory if it doenst already exist
    gen.create_output_directories(global_output_directory)

    # get the sequences
    if extract_sequences:
        lincRNA_single_exon_bed = "{0}/lincrna/lincRNA.single_exon.bed".format(output_directory)
        lincRNA_single_exon_fasta = "{0}/lincrna/lincRNA.single_exon.fasta".format(output_directory)
        lincRNA_single_exon_families = "{0}/lincrna/lincRNA.single_exon_families.bed".format(output_directory)
        lincRNA_multi_exon_bed = "{0}/lincrna/lincRNA.multi_exon.bed".format(output_directory)
        lincRNA_multi_exon_intron_bed = "{0}/lincrna/lincRNA.multi_exon.introns.bed".format(output_directory)
        lincRNA_multi_exon_fasta = "{0}/lincrna/lincRNA.multi_exon.fasta".format(output_directory)
        lincRNA_multi_exon_exons_fasta = "{0}/lincrna/lincRNA.multi_exon.exons.fasta".format(output_directory)
        lincRNA_multi_exon_intron_fasta = "{0}/lincrna/lincRNA.multi_exon.introns.fasta".format(output_directory)
        lincRNA_multi_exon_families = "{0}/lincrna/lincRNA.multi_exon_families.bed".format(output_directory)
        cont.extract_lincRNA_sequences(input_bed, input_fasta, lincRNA_single_exon_bed, lincRNA_multi_exon_bed, lincRNA_single_exon_fasta, lincRNA_multi_exon_fasta, lincRNA_multi_exon_intron_bed, lincRNA_multi_exon_intron_fasta, lincRNA_single_exon_families, lincRNA_multi_exon_families, clean_run = None)

    # clean the alignments to get in usable form
    # might need this
    if clean_alignments:
        output_exon_file = "{0}/clean_exon_alignments.fasta"
        output_intron_file = "{0}/clean_intron_alignments.fasta"
        ltests.clean_alignments(input_bed, input_fasta, output_exon_file, output_intron_file)

    if calc_gc:
        output_file = "{0}/{1}_gc.csv".format(global_output_directory, output_prefix)
        ltests.calc_gc(input_fasta, output_file, families_file = families_file)

    # orf length test
    if sim_orf_lengths:
        sim_orf_length_output_file = "{0}/{1}sim_orf_lengths.csv".format(global_output_directory, output_prefix)
        if families_file:
            sim_orf_length_z_file = "{0}/{1}sim_orf_lengths_zs_grouped.csv".format(global_output_directory, output_prefix)
        else:
            sim_orf_length_z_file = "{0}/{1}sim_orf_lengths_zs.csv".format(global_output_directory, output_prefix)
        # run the test
        simopc.sim_orf_length(input_fasta, required_simulations, sim_orf_length_output_file)
        ltests.process_length_sim(sim_orf_length_output_file, sim_orf_length_z_file, families_file = families_file)

    if sim_orf_lengths_masked:
        masked_output_file = "{0}_{1}_masked.csv".format(input_file.split(".")[0],  motif_file.split("/")[-1].split(".")[0])
        # run the test
        simopc.sim_orf_length_masked(input_fasta, required_simulations, motif_file, input_file, controls_dir, masked_output_file, families_file = families_file)

    # stop density test
    if sim_stop_density:
        local_output_directory = "{0}/stop_density".format(global_output_directory)
        gen.create_output_directories(local_output_directory)
        if families_file:
            sim_stop_density_output_dir = "{0}/{1}_stop_density_simulation_all_genes_grouped_families".format(local_output_directory, output_prefix)
            sim_stop_density_output_file = "{0}/{1}_stop_density_simulation_all_genes_grouped_families.csv".format(local_output_directory, output_prefix)
            runs = 10
        else:
            sim_stop_density_output_dir = "{0}/{1}_stop_density_simulation_all_genes".format(local_output_directory, output_prefix)
            sim_stop_density_output_file = "{0}/{1}_stop_density_simulation_all_genes.csv".format(local_output_directory, output_prefix)
            runs = 1

        gen.create_output_directories(sim_stop_density_output_dir)

        for run in list(range(runs)):
            output_file =  "{0}/stop_density_simulation_{1}.csv".format(sim_stop_density_output_dir, run + 1)
            ltests.sim_stop_density(input_fasta, output_file, simulations = int(required_simulations), families_file = families_file)
        # process the outputs
        ltests.process_sim_stop_density_outputs(sim_stop_density_output_dir, sim_stop_density_output_file)

    # within genes
    if sim_stop_density_within_genes:
        local_output_directory = "{0}/stop_density".format(global_output_directory)
        gen.create_output_directories(local_output_directory)
        if families_file:
            sim_stop_density_within_gene_output_dir = "{0}/{1}_stop_density_simulation_within_genes_grouped_families".format(local_output_directory, output_prefix)
            sim_stop_density_within_gene_output_file = "{0}/{1}_stop_density_simulation_within_genes_grouped_families.csv".format(local_output_directory, output_prefix)
            runs = 10
        else:
            sim_stop_density_within_gene_output_dir = "{0}/{1}_stop_density_simulation_within_genes".format(local_output_directory, output_prefix)
            sim_stop_density_within_gene_output_file = "{0}/{1}_stop_density_simulation_within_genes.csv".format(local_output_directory, output_prefix)
            runs = 1
        gen.create_output_directories(sim_stop_density_within_gene_output_dir)
        for run in list(range(runs)):
            output_file =  "{0}/stop_density_simulation_{1}.csv".format(sim_stop_density_within_gene_output_dir, run + 1)
            ltests.sim_stop_density_within_genes(input_fasta, output_file, simulations = int(required_simulations), families_file = families_file)

        # process the outputs
        ltests.process_sim_stop_density_within_gene_outputs(sim_stop_density_within_gene_output_dir, sim_stop_density_within_gene_output_file)

    # stop density test in the introns
    if sim_stop_density_introns:
        local_output_directory = "{0}/stop_density".format(global_output_directory)
        gen.create_output_directories(local_output_directory)
        if families_file:
            sim_stop_density_output_dir = "{0}/{1}_stop_density_introns_simulation_all_genes_grouped_families".format(local_output_directory, output_prefix)
            sim_stop_density_output_file = "{0}/{1}_stop_density_introns_simulation_all_genes_grouped_families.csv".format(local_output_directory, output_prefix)
            runs = 10
        else:
            sim_stop_density_output_dir = "{0}/{1}_stop_density_introns_simulation_all_genes".format(local_output_directory, output_prefix)
            sim_stop_density_output_file = "{0}/{1}_stop_density_introns_simulation_all_genes.csv".format(local_output_directory, output_prefix)
            runs = 1
        gen.create_output_directories(sim_stop_density_output_dir)

        for run in list(range(runs)):
            output_file =  "{0}/stop_density_simulation_{1}.csv".format(sim_stop_density_output_dir, run + 1)
            ltests.sim_stop_density(input_fasta, output_file, simulations = int(required_simulations), families_file = families_file, introns = True, input_fasta2 = input_fasta2)
        # process the outputs
        ltests.process_sim_stop_density_outputs(sim_stop_density_output_dir, sim_stop_density_output_file)


    # remove motifs and test
    if sim_stop_density_removed_motifs:
        local_output_directory = "{0}/stop_density".format(global_output_directory)
        gen.create_output_directories(local_output_directory)
        if families_file:
            sim_output_dir = "{0}/{1}_{2}_stop_density_simulation_all_genes_grouped_families_removed_motifs".format(local_output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            sim_output_file = "{0}/{1}_{2}_stop_density_simulation_all_genes_grouped_families_removed_motifs.csv".format(local_output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            runs = 10
        else:
            sim_output_dir = "{0}/{1}_{2}_stop_density_simulation_all_genes".format(local_output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            sim_output_file = "{0}/{1}_{2}_stop_density_simulation_all_genes.csv".format(local_output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            runs = 1
        # remove any previous runs
        gen.remove_directory(sim_output_dir)
        gen.create_output_directories(sim_output_dir)

        for run in list(range(runs)):
            run_output_file =  "{0}/stop_density_simulation_{1}.csv".format(sim_output_dir, run + 1)
            ltests.sim_stop_density_removed_motifs(input_fasta, run_output_file, motif_file, simulations = int(required_simulations), families_file = families_file)
        # process the outputs
        ltests.process_sim_stop_density_outputs(sim_output_dir, sim_output_file, reverse = True)

    # remove motifs and test within seqs
    if sim_stop_density_removed_motifs_sim_seqs:
        local_output_directory = "{0}/stop_density".format(global_output_directory)
        gen.create_output_directories(local_output_directory)
        if families_file:
            sim_output_dir = "{0}/{1}_{2}_stop_density_simulation_grouped_families_removed_motifs_seq_sim".format(local_output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            sim_output_file = "{0}/{1}_{2}_stop_density_simulation_grouped_families_removed_motifs_seq_sim.csv".format(local_output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            runs = 10
        else:
            sim_output_dir = "{0}/{1}_{2}_stop_density_simulation_all_genes_seq_sim".format(local_output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            sim_output_file = "{0}/{1}_{2}_stop_density_simulation_all_genes_seq_sim.csv".format(local_output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            runs = 1
        # remove any previous runs
        gen.remove_directory(sim_output_dir)
        gen.create_output_directories(sim_output_dir)

        for run in list(range(runs)):
            run_output_file =  "{0}/stop_density_simulation_{1}.csv".format(sim_output_dir, run + 1)
            ltests.sim_stop_density_removed_motifs_seq_sim(input_fasta, run_output_file, motif_file, controls_dir, simulations = int(required_simulations), families_file = families_file)
        # process the outputs
        ltests.process_sim_stop_density_outputs(sim_output_dir, sim_output_file, reverse = True)


    if sim_stop_density_diff:
        local_output_directory = "{0}/stop_density".format(global_output_directory)
        gen.create_output_directories(local_output_directory)
        if families_file:
            sim_output_dir = "{0}/{1}_{2}_stop_density_diff_grouped_families".format(global_output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            sim_output_file = "{0}/{1}_{2}_stop_density_stop_density_diff_grouped_families.csv".format(global_output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            runs = 10
        else:
            sim_output_dir = "{0}/{1}_{2}_stop_density_stop_density_diff_all_genes".format(global_output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            sim_output_file = "{0}/{1}_{2}_stop_density_stop_density_diff_all_genes.csv".format(global_output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            runs = 1
        # remove any previous runs
        gen.remove_directory(sim_output_dir)
        gen.create_output_directories(sim_output_dir)

        for run in list(range(runs)):
            run_output_file =  "{0}/stop_density_simulation_{1}.csv".format(sim_output_dir, run + 1)
            ltests.sim_stop_density_diff(input_fasta, run_output_file, motif_file, controls_dir, simulations = int(required_simulations), families_file = families_file)
        # process the outputs
        ltests.process_sim_stop_density_diffs(sim_output_dir, sim_output_file, greater_than = False)

    # get density in exons and introns
    if exon_intron_density:
        local_output_directory = "{0}/stop_density".format(global_output_directory)
        gen.create_output_directories(local_output_directory)
        output_file = "{0}/exon_intron_stop_density.csv".format(local_output_directory)
        ltests.exon_intron_stop_density(input_fasta, input_fasta2, output_file, families_file = families_file)


    # test whether there is an excess in flanks
    if excess_test:
        gen.check_files_exists([input_fasta, motif_file])
        # local output directory
        local_output_directory = "{0}/stop_excesses".format(global_output_directory)
        gen.create_output_directories(local_output_directory)
        # if the families file exists, group by family
        if families_file:
            excess_test_output_file = "{0}/{1}_stop_codon_excesses_grouped.csv".format(local_output_directory, motif_file.split("/")[-1].split(".")[0])
        else:
            excess_test_output_file = "{0}/{1}_stop_codon_excesses.csv".format(local_output_directory, motif_file.split("/")[-1].split(".")[0])
        # run the test
        ltests.excess_test(input_fasta, motif_file, excess_test_output_file, simulations = required_simulations, families_file = families_file)

    # upstream from the atg
    if upstream_atg:
        output_file = "{0}/tests/lincrna/stop_density/upstream_atg_stop_density.csv".format(global_output_directory)
        ltests.upstream_atg(input_fasta, output_file, simulations = int(required_simulations), families_file = families_file)




    # test hits to seqs
    if seq_hits_linc:
        local_output_dir = "{0}/ese_hits".format(global_output_directory)
        if output_prefix:
            tests_output_dir = "{0}/{1}_{2}".format(local_output_dir, output_prefix, motif_file.split("/")[-1].split(".")[0])
            final_output_file = "{0}/{1}_{2}_processed.csv".format(local_output_dir, output_prefix, motif_file.split("/")[-1].split(".")[0])
        else:
            tests_output_dir = "{0}/{1}_{2}".format(local_output_dir, output_prefix, motif_file.split("/")[-1].split(".")[0])
            final_output_file = "{0}/{1}_{2}_processed.csv".format(local_output_dir, output_prefix, motif_file.split("/")[-1].split(".")[0])
        gen.create_output_directories(tests_output_dir)

        runs = 10
        for run in range(runs):
            if output_prefix:
                output_file = "{0}/{1}_{2}_hits_{3}.csv".format(tests_output_dir, output_prefix, motif_file.split("/")[-1].split(".")[0], run+1)
            else:
                output_file = "{0}/{1}_hits_{2}.csv".format(tests_output_dir, motif_file.split("/")[-1].split(".")[0], run+1)
            mto.calc_seq_hits_linc(input_fasta, output_file, motif_file, controls_dir, required_simulations = required_simulations, families_file = families_file)
        mto.process_seq_hits_linc(tests_output_dir, final_output_file)

    if excess_length_thresholds:
        local_output_dir = "{0}/orf_length_thresholds".format(global_output_directory)
        gen.create_output_directories(local_output_dir)
        ltests.orf_exceed_length_threshold(input_fasta, local_output_directory, required_simulations = required_simulations, families_file = families_file)



if __name__ == "__main__":
    main()
