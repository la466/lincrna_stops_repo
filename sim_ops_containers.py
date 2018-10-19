import generic as gen
import sim_ops as simo
import seq_ops as seqo
import file_ops as fo
import ops_containers as opsc
import ops
import os


def run_simulation_function(required_simulations, sim_args, function_to_run, parallel = True):
    """
    Wrapper to run simulation function

    Args:
        required_simulations (int): the number of times to run simulation
        sim_args (list): list of arguments for the simulation function
        function_to_run (func): simulation_function
        parallel (bool): if true, run in parallel

    Returns:
        outputs (list): list containing the result of the simulation
    """

    # get a list of simulations to iterate over
    simulations = list(range(required_simulations))
    # run the simulations
    if parallel:
        # add foo to argument list for parallelisation
        sim_args.insert(0, "foo")
        workers = os.cpu_count() - 2
        processes = gen.run_in_parallel(simulations, sim_args, function_to_run, workers = workers)
        outputs = []
        for process in processes:
            outputs.extend(process.get())
    else:
        outputs = simo.sim_orf_lengths(simulations, *sim_args)

    return outputs


def sim_coding_exon_stop_counts_regions(genome_fasta, gtf, output_directory, output_file, simulations, clean_run=None):
    """
    """

    sequence_output_directory = "{0}/sequence_files".format(output_directory)
    if clean_run:
        gen.remove_directory(sequence_output_directory)

    # create the output directory if it doesnt exist
    gen.create_output_directories(sequence_output_directory)

    # get the coding and non coding exons
    coding_exons_bed = "{0}/coding_exons.bed".format(output_directory)
    non_coding_exons_bed = "{0}/non_coding_exons.bed".format(output_directory)
    coding_exons_fasta = "{0}/coding_exons.fasta".format(output_directory)
    non_coding_exons_fasta = "{0}/non_coding_exons.fasta".format(output_directory)
    # get the coding and non coding exons
    if clean_run or not os.path.isfile(non_coding_exons_fasta) or not os.path.isfile(coding_exons_fasta):
        print("Getting the coding and non coding exons...")
        opsc.get_non_coding_exons(genome_fasta, gtf, coding_exons_bed, non_coding_exons_bed, coding_exons_fasta, non_coding_exons_fasta, sequence_output_directory)

    # # get which frame each of the exons resides
    final_filtered_cds_file = "{0}/final_filtered_cds.fasta".format(sequence_output_directory)
    coding_exons_frames_file = "{0}/coding_exons_reading_frames.fasta".format(output_directory)
    if clean_run or not os.path.isfile(coding_exons_frames_file):
        print("Getting the reading frame each exon starts in...")
        ops.get_exon_reading_frame(coding_exons_fasta, final_filtered_cds_file, coding_exons_frames_file)

    # simulate the coding exons
    simulate_coding_exon_regions(coding_exons_fasta, coding_exons_frames_file, output_file, simulations)


def simulate_coding_exon_regions(exons_fasta, exons_frames_file, output_file, required_simulations, parallel=True, seeds=None, seq_seeds=None):

    window_start, window_end = 3, 69

    # get the dinucleotide content
    names, seqs_list = gen.read_fasta(exons_fasta)
    exon_frame_names, exon_frames = gen.read_fasta(exons_frames_file)

    # get the frames of the exons
    exon_frames_list = {name: int(exon_frames[i]) for i, name in enumerate(exon_frame_names)}

    # return a dictionary and list of unique sequences
    seqs = {name: seqs_list[i] for i, name in enumerate(names) if len(seqs_list[i]) > window_end*2 }
    unique_seqs = [seqs[i] for i in seqs]

    # get the dicnucleotide and nucleotide content of the sequences
    dinucleotide_content = seqo.get_dinucleotide_content(unique_seqs)
    nucleotide_content = seqo.get_nucleotide_content(unique_seqs)

    real_regions_stop_counts = ops.get_region_stop_counts(seqs, window_start, window_end)

    # create a temporary output directory
    temp_dir = "temp_data_regions"
    gen.create_output_directories(temp_dir)

    sim_args = [seqs, exon_frames_list, dinucleotide_content, nucleotide_content, window_start, window_end, temp_dir, seeds, seq_seeds]
    # run the simulation
    outputs = run_simulation_function(required_simulations, sim_args, simo.sim_exon_region_counts, parallel=parallel)

    # get temp filelist so we can order simulants for tests
    temp_filelist = fo.order_temp_files(outputs)

    # write to file
    with open(output_file, "w") as outfile:
        # write the header
        outfile.write("sim_no,flank_count,core_count\n")
        # write the real line
        outfile.write("real,{0},{1}\n".format(real_regions_stop_counts[0], real_regions_stop_counts[1]))
        # write the simulants
        for sim in sorted(temp_filelist):
            # get the simulation number
            data = gen.read_many_fields(temp_filelist[sim], "\t")[0]
            outfile.write("{0},{1},{2}\n".format("sim_{0}".format(sim), data[0], data[1]))

    # remove the temp files
    gen.remove_directory(temp_dir)



def sim_orf_length(seq_fasta, required_simulations, output_file, parallel=True, seeds=None, seq_seeds=None):
    """
    Simulation to look at the lengths of ORFs in sequences
    compared with simulated null sequences

    Args:
        seq_fasta (str): path to fasta file containing sequences
        required_simulations (int): number of simulations to do
        output_file (str): path to output file
        parallel (bool): opt whether to use multiprocessing
        seeds (list): seed numbers for the simulations
        seq_seeds (list): a list of lists containing seeds for testing
    """

    names, seqs_list = gen.read_fasta(seq_fasta)

    # return a dictionary and list of unique sequences
    seqs = seqo.get_unique_seqs(names, seqs_list)
    unique_seqs = [seqs[i] for i in seqs]

    # get the dicnucleotide and nucleotide content of the sequences
    dinucleotide_content = seqo.get_dinucleotide_content(unique_seqs)
    nucleotide_content = seqo.get_nucleotide_content(unique_seqs)

    # get the longest orfs and gc content for the real sequences
    longest_orfs = seqo.get_longest_orfs(seqs)
    gc_contents = {}
    for seq in seqs:
        gc_contents[seq] = seqo.calc_seq_gc(seqs[seq])

    # create a temporary output directory
    temp_dir = "temp_data_length_sim"
    gen.create_output_directories(temp_dir)

    # run the simulations
    simulations = list(range(required_simulations))
    sim_args = [seqs, dinucleotide_content, nucleotide_content, temp_dir, seeds, seq_seeds]

    # run the simulations
    if parallel:
        # add foo to argument list for parallelisation
        sim_args.insert(0, "foo")
        workers = os.cpu_count() - 2
        processes = gen.run_in_parallel(simulations, sim_args, simo.sim_orf_lengths, workers = workers)
        outputs = []
        for process in processes:
            outputs.extend(process.get())
    else:
        outputs = simo.sim_orf_lengths(simulations, *sim_args)

    # get temp filelist so we can order simulants for tests
    temp_filelist = fo.order_temp_files(outputs)

    # write to file
    with open(output_file, "w") as outfile:
        # write the header
        outfile.write("sim_no,{0}\n".format(",".join(name for name in sorted(seqs))))
        # write the gc content
        outfile.write("gc,{0}\n".format(",".join(gen.stringify([gc_contents[id] for id in sorted(gc_contents)]))))
        # write the real line
        outfile.write("real,{0}\n".format(",".join(gen.stringify([longest_orfs[id] for id in sorted(longest_orfs)]))))
        # write the simulants
        for sim in sorted(temp_filelist):
            # get the simulation number
            data = seqo.fasta_to_dict(temp_filelist[sim])
            outfile.write("{0},{1}\n".format("sim_{0}".format(sim), ",".join([data[id] for id in sorted(data)])))

    # remove the temp files
    gen.remove_directory(temp_dir)


def sim_stop_count(seq_fasta, required_simulations, output_file, parallel=True, seeds=None, seq_seeds=None):
    """
    Simulation asking whether sequences have fewer stop codons than
    nucleotide matched controls.

    Args:
        seq_fasta (str): path to fasta file containing sequences
        required_simulations (int): number of simulations to do
        output_file (str): path to output file
        parallel (bool): opt whether to use multiprocessing
        seeds (list): seed numbers for the simulations
        seq_seeds (list): a list of lists containing seeds for testing
    """

    names, seqs_list = gen.read_fasta(seq_fasta)

    # return a dictionary and list of unique sequences
    seqs = seqo.get_unique_seqs(names, seqs_list)
    unique_seqs = [seqs[i] for i in seqs]

    # get the gc contents of the sequences
    gc_contents = {}
    for seq in seqs:
        gc_contents[seq] = seqo.calc_seq_gc(seqs[seq])

    # get the dicnucleotide and nucleotide content of the sequences
    dinucleotide_content = seqo.get_dinucleotide_content(unique_seqs)
    nucleotide_content = seqo.get_nucleotide_content(unique_seqs)

    # get the stop codon counts for the real sequences
    stop_counts = seqo.get_stop_counts(seqs)

    # create a temporary output directory
    temp_dir = "temp_data_stop_count_sim"
    gen.create_output_directories(temp_dir)

    # run the simulations
    simulations = list(range(required_simulations))
    sim_args = [seqs, dinucleotide_content, nucleotide_content, temp_dir, seeds, seq_seeds]

    # run the simulations
    if parallel:
        # add foo to argument list for parallelisation
        sim_args.insert(0, "foo")
        workers = os.cpu_count() - 2
        processes = gen.run_in_parallel(simulations, sim_args, simo.sim_stop_counts, workers = workers)
        outputs = []
        for process in processes:
            outputs.extend(process.get())
    else:
        outputs = simo.sim_stop_counts(simulations, *sim_args)

    # get temp filelist so we can order simulants for tests
    temp_filelist = fo.order_temp_files(outputs)

    # write to file
    with open(output_file, "w") as outfile:
        # write the header
        outfile.write("sim_no,{0}\n".format(",".join(name for name in sorted(seqs))))
        # write the gc content
        outfile.write("gc,{0}\n".format(",".join(gen.stringify([gc_contents[id] for id in sorted(gc_contents)]))))
        # write the real line
        outfile.write("real,{0}\n".format(",".join(gen.stringify([stop_counts[id] for id in sorted(stop_counts)]))))
        # write the simulants
        for sim in sorted(temp_filelist):
            # get the simulation number
            data = seqo.fasta_to_dict(temp_filelist[sim])
            outfile.write("{0},{1}\n".format("sim_{0}".format(sim), ",".join([data[id] for id in sorted(data)])))

    # remove the temp files
    gen.remove_directory(temp_dir)
