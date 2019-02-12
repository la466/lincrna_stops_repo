import generic as gen
import containers as cont
import lincRNA_tests as ltests
import sim_ops_containers as simopc
import main_tests_ops as mto
import seq_ops as seqo
import time
import os

def main():

    arguments = ["input_bed", "input_fasta", "output_directory", "required_simulations", "extract_sequences", "clean_run", "density_sim", "get_exon_dint_controls", "get_intron_dint_controls", "exon_region_density", "compare_stop_density", "sim_orf_lengths", "sim_stop_density", "sim_stop_density_within_genes"]
    description = "Container for analysis on lincRNAs"
    args = gen.parse_arguments(description, arguments, flags = [4,5,6,7,8,9,10,11,12,13], ints=[], opt_flags = [3])
    input_bed, \
    input_fasta, \
    output_directory, \
    required_simulations, \
    extract_sequences, \
    clean_run, density_sim,  \
    get_exon_dint_controls, \
    get_intron_dint_controls, \
    exon_region_density, \
    compare_stop_density, \
    sim_orf_lengths, \
    sim_stop_density, \
    sim_stop_density_within_genes = \
    args.input_bed, \
    args.input_fasta, \
    args.output_directory, \
    args.required_simulations, \
    args.extract_sequences, \
    args.clean_run, \
    args.density_sim, \
    args.get_exon_dint_controls, \
    args.get_intron_dint_controls, \
    args.exon_region_density, \
    args.compare_stop_density, \
    args.sim_orf_lengths, \
    args.sim_stop_density, \
    args.sim_stop_density_within_genes

    lincrna_output_directory = "{0}/tests/lincrna".format(output_directory)
    gen.create_output_directories(lincrna_output_directory)

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

    if extract_sequences or not os.path.isfile(lincRNA_single_exon_fasta) or not os.path.isfile(lincRNA_multi_exon_fasta):
        cont.extract_lincRNA_sequences(input_bed, input_fasta, lincRNA_single_exon_bed, lincRNA_multi_exon_bed, lincRNA_single_exon_fasta, lincRNA_multi_exon_fasta, lincRNA_multi_exon_intron_bed, lincRNA_multi_exon_intron_fasta, lincRNA_single_exon_families, lincRNA_multi_exon_families, clean_run = clean_run)

    # multi_exon_sequence_density_simulation = "{0}/lincrna/tests/lincRNA.multi_exon.full_sequence.stop_density_simulation.csv".format(output_directory)
    # if stop_density_test:
    #     ltests.stop_density_test(lincRNA_multi_exon_fasta, multi_exon_sequence_density_simulation, required_simulations, families_file = lincRNA_multi_exon_families)
    # single_exon_sequence_density_simulation = "{0}/lincrna/tests/lincRNA.single_exon.full_sequence.stop_density_simulation.csv".format(output_directory)
    # if stop_density_test:
    #     ltests.stop_density_test(lincRNA_single_exon_fasta, single_exon_sequence_density_simulation, required_simulations, families_file = lincRNA_single_exon_families)

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

    if density_sim:
        ltests.density_simulation(lincRNA_multi_exon_exons_fasta, lincRNA_multi_exon_intron_fasta, required_simulations, families_file = lincRNA_multi_exon_families)

    sim_orf_length_output_file = "{0}/sim_orf_lengths.csv".format(lincrna_output_directory)
    lincRNA_length_output_file = "{0}/lincRNA_lengths.csv".format(lincrna_output_directory)
    sim_orf_length_z_file = "{0}/sim_orf_lengths_zs.csv".format(lincrna_output_directory)
    if sim_orf_lengths:
        if not os.path.isfile(sim_orf_length_output_file) or clean_run:
            simopc.sim_orf_length(lincRNA_multi_exon_fasta, required_simulations, sim_orf_length_output_file)
        ltests.process_length_sim(sim_orf_length_output_file, sim_orf_length_z_file)
        ltests.calculate_lengths(lincRNA_multi_exon_fasta, lincRNA_length_output_file, families_file = lincRNA_multi_exon_families)

    sim_stop_density_output_dir = "{0}/stop_density_simulations_all_genes".format(lincrna_output_directory)
    sim_stop_density_output_file = "{0}/stop_density_simulation_all_genes_outputs.csv".format(lincrna_output_directory)
    if sim_stop_density:
        runs = 10
        if clean_run:
            gen.remove_directory(sim_stop_density_output_dir)
        gen.create_output_directories(sim_stop_density_output_dir)
        # if we need to run the simulations
        if len(os.listdir(sim_stop_density_output_dir)) < runs:
            required_runs = list(range(runs - len(os.listdir(sim_stop_density_output_dir))))
            for run in required_runs:
                output_file =  "{0}/stop_density_simulation_{1}.csv".format(sim_stop_density_output_dir, run + 1)
                ltests.sim_stop_density(lincRNA_multi_exon_fasta, output_file, simulations = int(required_simulations), families_file = lincRNA_multi_exon_families)
        # process the outputs
        ltests.process_sim_stop_density_outputs(sim_stop_density_output_dir, sim_stop_density_output_file)

    sim_stop_density_within_gene_output_dir = "{0}/stop_density_simulations_within_genes".format(lincrna_output_directory)
    sim_stop_density_within_gene_output_file = "{0}/stop_density_simulation_within_genes_outputs.csv".format(lincrna_output_directory)
    if sim_stop_density_within_genes:
        runs = 10
        if clean_run:
            gen.remove_directory(sim_stop_density_within_gene_output_dir)
        gen.create_output_directories(sim_stop_density_within_gene_output_dir)
        # if we need to run the simulations
        if len(os.listdir(sim_stop_density_within_gene_output_dir)) < runs:
            required_runs = list(range(runs - len(os.listdir(sim_stop_density_within_gene_output_dir))))
            for run in required_runs:
                output_file =  "{0}/stop_density_simulation_{1}.csv".format(sim_stop_density_within_gene_output_dir, run + 1)
                ltests.sim_stop_density_within_genes(lincRNA_multi_exon_fasta, output_file, simulations = int(required_simulations), families_file = lincRNA_multi_exon_families)
        # process the outputs
        ltests.process_sim_stop_density_within_gene_outputs(sim_stop_density_within_gene_output_dir, sim_stop_density_within_gene_output_file)


if __name__ == "__main__":
    main()
