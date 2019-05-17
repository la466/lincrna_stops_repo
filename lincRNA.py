import generic as gen
import containers as cont
import lincRNA_tests as ltests
import sim_ops_containers as simopc
import sequence_ops as sequo
import main_tests_ops as mto
import seq_ops as seqo
import time
import os

def main():

    arguments = ["input_bed", "input_fasta", "output_directory", "input_fasta2", "required_simulations", "motif_file", "families_file", "output_prefix", "sim_dir", "extract_sequences", "clean_run", "calc_gc", "density_sim", "get_exon_dint_controls", "get_intron_dint_controls", "exon_region_density", "compare_stop_density", "sim_orf_lengths", "sim_stop_density", "sim_stop_density_introns", "sim_stop_density_within_genes", "sim_stop_density_removed_motifs", "sim_stop_density_removed_motifs_sim_seqs", "sim_stop_density_diff", "motif_nd", "excess_test", "single_exon", "substitution_rate", "substitution_rate_motif", "dinucleotide_substitution_rate", "motif_overlap", "motif_overlap_density", "clean_alignments", "upstream_atg"]
    description = "Container for analysis on lincRNAs"
    args = gen.parse_arguments(description, arguments, flags = [9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33], opt_flags = [3,4,5,6,7,8])
    input_bed, \
    input_fasta, \
    output_directory, \
    input_fasta2, \
    required_simulations, \
    motif_file, \
    families_file, \
    output_prefix, \
    sim_dir, \
    extract_sequences, \
    clean_run, \
    calc_gc, \
    density_sim,  \
    get_exon_dint_controls, \
    get_intron_dint_controls, \
    exon_region_density, \
    compare_stop_density, \
    sim_orf_lengths, \
    sim_stop_density, \
    sim_stop_density_introns, \
    sim_stop_density_within_genes, \
    sim_stop_density_removed_motifs, \
    sim_stop_density_removed_motifs_sim_seqs, \
    sim_stop_density_diff, \
    motif_nd, \
    excess_test, \
    single_exon,\
    substitution_rate,\
    substitution_rate_motif,\
    dinucleotide_substitution_rate,\
    motif_overlap, \
    motif_overlap_density, \
    clean_alignments, \
    upstream_atg = \
    args.input_bed, \
    args.input_fasta, \
    args.output_directory, \
    args.input_fasta2, \
    args.required_simulations, \
    args.motif_file, \
    args.families_file, \
    args.output_prefix, \
    args.sim_dir, \
    args.extract_sequences, \
    args.clean_run, \
    args.calc_gc, \
    args.density_sim, \
    args.get_exon_dint_controls, \
    args.get_intron_dint_controls, \
    args.exon_region_density, \
    args.compare_stop_density, \
    args.sim_orf_lengths, \
    args.sim_stop_density, \
    args.sim_stop_density_introns, \
    args.sim_stop_density_within_genes, \
    args.sim_stop_density_removed_motifs, \
    args.sim_stop_density_removed_motifs_sim_seqs, \
    args.sim_stop_density_diff, \
    args.motif_nd, \
    args.excess_test, \
    args.single_exon, \
    args.substitution_rate, \
    args.substitution_rate_motif, \
    args.dinucleotide_substitution_rate, \
    args.motif_overlap, \
    args.motif_overlap_density, \
    args.clean_alignments, \
    args.upstream_atg

    # make required simultions an int
    required_simulations = int(required_simulations) if required_simulations else None

    lincrna_output_directory = "{0}/tests/lincrna".format(output_directory)
    # gen.create_output_directories(lincrna_output_directory)

    # set a start time
    start = time.time()

    # create the output_directory if it doenst already exist
    gen.create_output_directories(output_directory)

    lincRNA_single_exon_bed = "{0}/lincrna/lincRNA.single_exon.bed".format(output_directory)
    lincRNA_single_exon_fasta = "{0}/lincrna/lincRNA.single_exon.fasta".format(output_directory)
    lincRNA_single_exon_families = "{0}/lincrna/lincRNA.single_exon_families.bed".format(output_directory)
    lincRNA_multi_exon_bed = "{0}/lincrna/lincRNA.multi_exon.bed".format(output_directory)
    lincRNA_multi_exon_intron_bed = "{0}/lincrna/lincRNA.multi_exon.introns.bed".format(output_directory)
    lincRNA_multi_exon_fasta = "{0}/lincrna/lincRNA.multi_exon.fasta".format(output_directory)
    lincRNA_multi_exon_exons_fasta = "{0}/lincrna/lincRNA.multi_exon.exons.fasta".format(output_directory)
    lincRNA_multi_exon_intron_fasta = "{0}/lincrna/lincRNA.multi_exon.introns.fasta".format(output_directory)
    lincRNA_multi_exon_families = "{0}/lincrna/lincRNA.multi_exon_families.bed".format(output_directory)
    # get the sequences

    if extract_sequences:
        cont.extract_lincRNA_sequences(input_bed, input_fasta, lincRNA_single_exon_bed, lincRNA_multi_exon_bed, lincRNA_single_exon_fasta, lincRNA_multi_exon_fasta, lincRNA_multi_exon_intron_bed, lincRNA_multi_exon_intron_fasta, lincRNA_single_exon_families, lincRNA_multi_exon_families, clean_run = clean_run)

    # multi_exon_sequence_density_simulation = "{0}/lincrna/tests/lincRNA.multi_exon.full_sequence.stop_density_simulation.csv".format(output_directory)
    # if stop_density_test:
    #     ltests.stop_density_test(lincRNA_multi_exon_fasta, multi_exon_sequence_density_simulation, required_simulations, families_file = lincRNA_multi_exon_families)
    # single_exon_sequence_density_simulation = "{0}/lincrna/tests/lincRNA.single_exon.full_sequence.stop_density_simulation.csv".format(output_directory)
    # if stop_density_test:
    #     ltests.stop_density_test(lincRNA_single_exon_fasta, single_exon_sequence_density_simulation, required_simulations, families_file = lincRNA_single_exon_families)


    if calc_gc:
        output_file = "{0}/{1}_gc.csv".format(output_directory, output_prefix)
        ltests.calc_gc(input_fasta, output_file, families_file = families_file)


    exon_dint_control_directory = "{0}/lincrna/exon_dinucleotide_controls".format(output_directory)
    if get_exon_dint_controls:
        simopc.generate_dint_controls(lincRNA_multi_exon_fasta, exon_dint_control_directory)

    intron_dint_control_directory = "{0}/lincrna/intron_dinucleotide_controls".format(output_directory)
    if get_intron_dint_controls:
        simopc.generate_dint_intron_controls(lincRNA_multi_exon_intron_fasta, intron_dint_control_directory)

    exon_region_density_file = "{0}/exon_region_densities.csv".format(lincrna_output_directory)
    if exon_region_density:
        mto.exon_region_density(lincRNA_multi_exon_fasta, lincRNA_multi_exon_exons_fasta, None, exon_region_density_file, families_file = lincRNA_multi_exon_families)

    # get the stop density
    exon_intron_density_file = "{0}/compare_exon_intron_stop_density.csv".format(lincrna_output_directory)
    if compare_stop_density:
        mto.compare_stop_density(lincRNA_multi_exon_exons_fasta, lincRNA_multi_exon_intron_fasta, exon_intron_density_file, families_file = lincRNA_multi_exon_families)

    # if density_sim:
    #     ltests.density_simulation(lincRNA_multi_exon_exons_fasta, lincRNA_multi_exon_intron_fasta, required_simulations, families_file = lincRNA_multi_exon_families)



    # sim_orf_length_z_file = "{0}/sim_orf_lengths_zs.csv".format(output_directory)
    if sim_orf_lengths:
        sim_orf_length_output_file = "{0}/{1}sim_orf_lengths.csv".format(output_directory, "{0}_".format(output_prefix) if output_prefix else None)
        if families_file:
            sim_orf_length_z_file = "{0}/{1}sim_orf_lengths_zs_grouped.csv".format(output_directory, "{0}_".format(output_prefix) if output_prefix else None)
        else:
            sim_orf_length_z_file = "{0}/{1}sim_orf_lengths_zs.csv".format(output_directory, "{0}_".format(output_prefix) if output_prefix else None)

        # # only run the simulation if not there or wanted
        # if not os.path.isfile(sim_orf_length_output_file) or clean_run:
        #     simopc.sim_orf_length(input_fasta, required_simulations, sim_orf_length_output_file)
        ltests.process_length_sim(sim_orf_length_output_file, sim_orf_length_z_file, families_file = families_file)




    if sim_stop_density:
        if families_file:
            sim_stop_density_output_dir = "{0}/stop_density/{1}_stop_density_simulation_all_genes_grouped_families".format(output_directory, output_prefix)
            sim_stop_density_output_file = "{0}/stop_density/{1}_stop_density_simulation_all_genes_grouped_families.csv".format(output_directory, output_prefix)
            runs = 10
        else:
            sim_stop_density_output_dir = "{0}/stop_density/{1}_stop_density_simulation_all_genes".format(output_directory, output_prefix)
            sim_stop_density_output_file = "{0}/stop_density/{1}_stop_density_simulation_all_genes.csv".format(output_directory, output_prefix)
            runs = 1
        # remove any previous runs
        if clean_run:
            gen.remove_directory(sim_stop_density_output_dir)
        gen.create_output_directories(sim_stop_density_output_dir)

        # # if we need to run the simulations
        # # if len(os.listdir(sim_stop_density_output_dir)) < runs:
        # # required_runs = list(range(runs - len(os.listdir(sim_stop_density_output_dir))))
        for run in list(range(runs)):
            output_file =  "{0}/stop_density_simulation_{1}.csv".format(sim_stop_density_output_dir, run + 1)
            ltests.sim_stop_density(input_fasta, output_file, simulations = int(required_simulations), families_file = families_file)
        # process the outputs
        ltests.process_sim_stop_density_outputs(sim_stop_density_output_dir, sim_stop_density_output_file)

    if sim_stop_density_introns:
        if families_file:
            sim_stop_density_output_dir = "{0}/stop_density/{1}_stop_density_introns_simulation_all_genes_grouped_families".format(output_directory, output_prefix)
            sim_stop_density_output_file = "{0}/stop_density/{1}_stop_density_introns_simulation_all_genes_grouped_families.csv".format(output_directory, output_prefix)
            runs = 10
        else:
            sim_stop_density_output_dir = "{0}/stop_density/{1}_stop_density_introns_simulation_all_genes".format(output_directory, output_prefix)
            sim_stop_density_output_file = "{0}/stop_density/{1}_stop_density_introns_simulation_all_genes.csv".format(output_directory, output_prefix)
            runs = 1
        # remove any previous runs
        if clean_run:
            gen.remove_directory(sim_stop_density_output_dir)
        gen.create_output_directories(sim_stop_density_output_dir)

        # # if we need to run the simulations
        # # if len(os.listdir(sim_stop_density_output_dir)) < runs:
        # # required_runs = list(range(runs - len(os.listdir(sim_stop_density_output_dir))))
        for run in list(range(runs)):
            output_file =  "{0}/stop_density_simulation_{1}.csv".format(sim_stop_density_output_dir, run + 1)
            ltests.sim_stop_density(input_fasta, output_file, simulations = int(required_simulations), families_file = families_file, introns = True, input_fasta2 = input_fasta2)
        # process the outputs
        ltests.process_sim_stop_density_outputs(sim_stop_density_output_dir, sim_stop_density_output_file)


    if sim_stop_density_within_genes:
        if families_file:
            sim_stop_density_within_gene_output_dir = "{0}/stop_density/{1}_stop_density_simulation_within_genes_grouped_families".format(output_directory, output_prefix)
            sim_stop_density_within_gene_output_file = "{0}/stop_density/{1}_stop_density_simulation_within_genes_grouped_families.csv".format(output_directory, output_prefix)
            runs = 10
        else:
            sim_stop_density_within_gene_output_dir = "{0}/stop_density/{1}_stop_density_simulation_within_genes".format(output_directory, output_prefix)
            sim_stop_density_within_gene_output_file = "{0}/stop_density/{1}_stop_density_simulation_within_genes.csv".format(output_directory, output_prefix)
            runs = 1

        # if clean_run:
        #     gen.remove_directory(sim_stop_density_within_gene_output_dir)
        # gen.create_output_directories(sim_stop_density_within_gene_output_dir)
        # # if we need to run the simulations
        # # if len(os.listdir(sim_stop_density_within_gene_output_dir)) < runs:
        # #     required_runs = list(range(runs - len(os.listdir(sim_stop_density_within_gene_output_dir))))
        # for run in list(range(runs)):
        #     output_file =  "{0}/stop_density_simulation_{1}.csv".format(sim_stop_density_within_gene_output_dir, run + 1)
        #     ltests.sim_stop_density_within_genes(input_fasta, output_file, simulations = int(required_simulations), families_file = families_file)
        # process the outputs
        ltests.process_sim_stop_density_within_gene_outputs(sim_stop_density_within_gene_output_dir, sim_stop_density_within_gene_output_file)


    if sim_stop_density_removed_motifs:
        if families_file:
            sim_output_dir = "{0}/stop_density/{1}_{2}_stop_density_simulation_all_genes_grouped_families_removed_motifs".format(output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            sim_output_file = "{0}/stop_density/{1}_{2}_stop_density_simulation_all_genes_grouped_families_removed_motifs.csv".format(output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            runs = 10
        else:
            sim_output_dir = "{0}/stop_density/{1}_{2}_stop_density_simulation_all_genes".format(output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            sim_output_file = "{0}/stop_density/{1}_{2}_stop_density_simulation_all_genes.csv".format(output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            runs = 1
        # remove any previous runs
        gen.remove_directory(sim_output_dir)
        gen.create_output_directories(sim_output_dir)

        # # # if we need to run the simulations
        # # if len(os.listdir(sim_output_dir)) < runs:
        # # required_runs = list(range(runs - len(os.listdir(sim_output_dir))))
        # for run in list(range(runs)):
        #     run_output_file =  "{0}/stop_density_simulation_{1}.csv".format(sim_output_dir, run + 1)
        #     ltests.sim_stop_density_removed_motifs(input_fasta, run_output_file, motif_file, simulations = int(required_simulations), families_file = families_file)
        # process the outputs
        ltests.process_sim_stop_density_outputs(sim_output_dir, sim_output_file, reverse = True)

    if sim_stop_density_removed_motifs_sim_seqs:
        if families_file:
            sim_output_dir = "{0}/stop_density/{1}_{2}_stop_density_simulation_grouped_families_removed_motifs_seq_sim".format(output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            sim_output_file = "{0}/stop_density/{1}_{2}_stop_density_simulation_grouped_families_removed_motifs_seq_sim.csv".format(output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            runs = 10
        else:
            sim_output_dir = "{0}/stop_density/{1}_{2}_stop_density_simulation_all_genes_seq_sim".format(output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            sim_output_file = "{0}/stop_density/{1}_{2}_stop_density_simulation_all_genes_seq_sim.csv".format(output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            runs = 1
        # remove any previous runs
        gen.remove_directory(sim_output_dir)
        gen.create_output_directories(sim_output_dir)

        # if we need to run the simulations
        # if len(os.listdir(sim_output_dir)) < runs:
        #     required_runs = list(range(runs - len(os.listdir(sim_output_dir))))
        #     for run in required_runs:
        for run in list(range(runs)):
            run_output_file =  "{0}/stop_density_simulation_{1}.csv".format(sim_output_dir, run + 1)
            ltests.sim_stop_density_removed_motifs_seq_sim(input_fasta, run_output_file, motif_file, sim_dir, simulations = int(required_simulations), families_file = families_file)
        # process the outputs
        ltests.process_sim_stop_density_outputs(sim_output_dir, sim_output_file, reverse = True)


    if sim_stop_density_diff:
        if families_file:
            sim_output_dir = "{0}/stop_density/{1}_{2}_stop_density_diff_grouped_families".format(output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            sim_output_file = "{0}/stop_density/{1}_{2}_stop_density_stop_density_diff_grouped_families.csv".format(output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            runs = 10
        else:
            sim_output_dir = "{0}/stop_density/{1}_{2}_stop_density_stop_density_diff_all_genes".format(output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            sim_output_file = "{0}/stop_density/{1}_{2}_stop_density_stop_density_diff_all_genes.csv".format(output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
            runs = 1
        # remove any previous runs
        gen.remove_directory(sim_output_dir)
        gen.create_output_directories(sim_output_dir)

        # if we need to run the simulations
        # if len(os.listdir(sim_output_dir)) < runs:
        #     required_runs = list(range(runs - len(os.listdir(sim_output_dir))))
        #     for run in required_runs:
        for run in list(range(runs)):
            run_output_file =  "{0}/stop_density_simulation_{1}.csv".format(sim_output_dir, run + 1)
            ltests.sim_stop_density_diff(input_fasta, run_output_file, motif_file, sim_dir, simulations = int(required_simulations), families_file = families_file)
        # process the outputs
        ltests.process_sim_stop_density_diffs(sim_output_dir, sim_output_file, greater_than = False)


    motif_nd_output_file = "{0}/motif_nd.csv".format(lincrna_output_directory)
    if motif_nd:
        ltests.calculcate_motif_nd(lincRNA_multi_exon_fasta, motif_file, motif_nd_output_file, simulations = int(required_simulations), families_file = lincRNA_multi_exon_families)

    if excess_test:
        gen.check_files_exists([input_fasta, motif_file])
        # if the families file exists, group by family
        if families_file:
            excess_test_output_file = "{0}/{1}_stop_codon_excesses_grouped.csv".format(output_directory, motif_file.split("/")[-1].split(".")[0])
        else:
            excess_test_output_file = "{0}/{1}_stop_codon_excesses.csv".format(output_directory, motif_file.split("/")[-1].split(".")[0])
        # run the test
        ltests.excess_test(input_fasta, motif_file, excess_test_output_file, simulations = required_simulations, families_file = families_file)

    if motif_overlap:
        gen.check_files_exists([input_fasta, motif_file, families_file])
        output_file = "{0}/motif_overlaps/{1}_{2}_motif_overlap_chisquare.csv".format(output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
        runs = 10
        ltests.motif_overlap_test(input_fasta, motif_file, output_file, runs = runs, families_file = families_file)

    if motif_overlap_density:
        gen.check_files_exists([input_fasta, motif_file, families_file])
        output_file = "{0}/motif_overlaps/{1}_{2}_motif_overlap_density.csv".format(output_directory, output_prefix, motif_file.split("/")[-1].split(".")[0])
        runs = 10
        ltests.motif_overlap_density_test(input_fasta, motif_file, output_file, runs = runs, families_file = families_file)

    if clean_alignments:
        output_exon_file = "{0}/clean_exon_alignments.fasta"
        output_intron_file = "{0}/clean_intron_alignments.fasta"
        ltests.clean_alignments(input_bed, input_fasta, output_exon_file, output_intron_file)

    if upstream_atg:
        output_file = "{0}/tests/lincrna/stop_density/upstream_atg_stop_density.csv".format(output_directory)
        ltests.upstream_atg(input_fasta, output_file, simulations = int(required_simulations), families_file = families_file)


if __name__ == "__main__":
    main()
