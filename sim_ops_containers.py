import generic as gen
import sim_ops as simo
import seq_ops as seqo
import file_ops as fo
import ops_containers as opsc
import ops
import os
import collections
import numpy as np

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
        outputs = function_to_run(simulations, *sim_args)

    return outputs


def sim_coding_exon_stop_counts_regions(genome_fasta, gtf, output_directory, output_file, simulations, clean_run=None):
    """
    Simulate the number of stop codons found in regions of coding exons.

    Args:
        genome_fasta (str): path to genome fasta file
        gtf (str): path to gtf file
        output_directory (str): path to output directory
        output_file (str): path to output file
        simulations (int): number of simulations to do
        clean_run (bool): if true, delete the sequence files compiled and start again
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


def sim_coding_exon_flanks_stop_counts(genome_fasta, gtf, output_directory, output_file, simulations, clean_run=None):
    """
    Simulate the number of stop codons found in coding exons flanks.

    Args:
        genome_fasta (str): path to genome fasta file
        gtf (str): path to gtf file
        output_directory (str): path to output directory
        output_file (str): path to output file
        simulations (int): number of simulations to do
        clean_run (bool): if true, delete the sequence files compiled and start again
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
    coding_exons_flanks_frames_file = "{0}/coding_exons_flanks_reading_frames.fasta".format(output_directory)
    if clean_run or not os.path.isfile(coding_exons_flanks_frames_file):
        print("Getting the reading frame each exon starts in...")
        ops.get_exon_flank_reading_frame(coding_exons_fasta, final_filtered_cds_file, coding_exons_flanks_frames_file)

    # simulate the coding exons
    simulate_coding_exon_flanks(coding_exons_fasta, coding_exons_flanks_frames_file, output_file, simulations)


def simulate_coding_exon_flanks(exons_fasta, exons_frames_file, output_file, required_simulations, parallel=True, seeds=None, seq_seeds=None):

    window_start, window_end = 3, 69

    # get the dinucleotide content
    names, seqs_list = gen.read_fasta(exons_fasta)
    exon_frame_names, exon_frame_starts = gen.read_fasta(exons_frames_file)

    # get the frames of the exons
    exon_frames_list = {}
    for i, name in enumerate(exon_frame_names):
        exon_frames_list[name] = [int(i) for i in exon_frame_starts[i].split(",")]

    # return a dictionary and list of unique sequences
    full_seqs = {name: seqs_list[i] for i, name in enumerate(names) if name in exon_frames_list}
    seqs = {name: [full_seqs[name][2:69], full_seqs[name][-69:-2]] for name in full_seqs}
    unique_seqs = []
    for name in seqs:
        unique_seqs.extend(seqs[name])

    # get the dicnucleotide and nucleotide content of the sequences
    dinucleotide_content = seqo.get_dinucleotide_content(unique_seqs)
    nucleotide_content = seqo.get_nucleotide_content(unique_seqs)


    # create a temporary output directory
    temp_dir = "temp_data_flanks"
    gen.create_output_directories(temp_dir)

    sim_args = [seqs, exon_frames_list, dinucleotide_content, nucleotide_content, window_start, window_end, temp_dir, seeds, seq_seeds]
    # run the simulation
    outputs = run_simulation_function(required_simulations, sim_args, simo.sim_exon_flank_counts, parallel=parallel)

    # get temp filelist so we can order simulants for tests
    filelist = []
    exons_to_exclude = []
    for i, output in enumerate(outputs):
        if len(output) and i % 2 == 0:
            filelist.extend(output)
        elif len(output):
            exons_to_exclude.extend(output)


    temp_filelist = fo.order_temp_files(filelist)
    exons_to_exclude = list(set(exons_to_exclude))


    # get a list of sequences that doesnt include the exons to avoid
    restricted_list_counts = collections.defaultdict(lambda: [])
    for name in seqs:
        for index, seq in enumerate(seqs[name]):
            query = "{0}.{1}".format(name, index)
            if query not in exons_to_exclude:
                restricted_list_counts[name].append(seqo.get_stop_counts([seqs[name][index]])[0])
            else:
                restricted_list_counts[name].append(np.nan)

    real_counts = [np.nansum([restricted_list_counts[name][0] for name in restricted_list_counts]),np.nansum([restricted_list_counts[name][1] for name in restricted_list_counts])]
    real_na_counts = np.nansum([len([restricted_list_counts[name][0] for name in restricted_list_counts if np.isnan(restricted_list_counts[name][0])]),len([restricted_list_counts[name][1] for name in restricted_list_counts if np.isnan(restricted_list_counts[name][1])])])

    # write to file
    with open(output_file, "w") as outfile:
        # write the header
        outfile.write("sim_no,five_prime_count,three_prime_count,flanks_not_counted\n")
        # write the real line
        outfile.write("real,{0},{1},{2}\n".format(real_counts[0], real_counts[1], real_na_counts))
        # write the simulants
        for sim in sorted(temp_filelist):
            # get the simulation number
            names, counts = gen.read_fasta(temp_filelist[sim])
            sim_counts = [[],[]]
            sim_na = 0
            for i, name in enumerate(names):
                exon_counts = counts[i].split(',')
                for index, count in enumerate(exon_counts):
                    query = "{0}.{1}".format(name, index)
                    if query not in exons_to_exclude:
                        if count != "nan":
                            sim_counts[index].append(int(count))
                        else:
                            sim_counts[index].append(np.nan)
                    else:
                        sim_na += 1
            outfile.write("{0},{1},{2},{3}\n".format("sim_{0}".format(sim), np.nansum(sim_counts[0]), np.nansum(sim_counts[1]), sim_na))

    # remove the temp files
    gen.remove_directory(temp_dir)


def simulate_coding_exon_regions(exons_fasta, exons_frames_file, output_file, required_simulations, parallel=True, seeds=None, seq_seeds=None):

    window_start, window_end = 3, 69

    # get the dinucleotide content
    names, seqs_list = gen.read_fasta(exons_fasta)
    exon_frame_names, exon_frame_starts = gen.read_fasta(exons_frames_file)

    # get the frames of the exons
    exon_frames_list = {}
    for i, name in enumerate(exon_frame_names):
        # print(name, exon_frame_starts[i])
        exon_frames_list[name] = int(exon_frame_starts[i])
    # exon_frames = {name: int(exon_frame_starts[i]) for i, name in enumerate(exon_frame_names)}

    # return a dictionary and list of unique sequences
    seqs = {name: seqs_list[i] for i, name in enumerate(names) if len(seqs_list[i]) > window_end*2}
    unique_seqs = [seqs[i] for i in seqs]

    # get the dicnucleotide and nucleotide content of the sequences
    dinucleotide_content = seqo.get_dinucleotide_content(unique_seqs)
    nucleotide_content = seqo.get_nucleotide_content(unique_seqs)


    # create a temporary output directory
    temp_dir = "temp_data_regions"
    gen.create_output_directories(temp_dir)

    sim_args = [seqs, exon_frames_list, dinucleotide_content, nucleotide_content, window_start, window_end, temp_dir, seeds, seq_seeds]
    # run the simulation
    outputs = run_simulation_function(required_simulations, sim_args, simo.sim_exon_region_counts, parallel=parallel)

    # get temp filelist so we can order simulants for tests
    filelist = []
    exons_to_exclude = []
    for i, output in enumerate(outputs):
        if len(output) and i % 2 == 0:
            filelist.extend(output)
        elif len(output):
            exons_to_exclude.extend(output)


    temp_filelist = fo.order_temp_files(filelist)
    exons_to_exclude = list(set(exons_to_exclude))

    # get a list of sequences that doesnt include the exons to avoid
    restricted_real_seqs = {name: seqs[name] for name in seqs if name not in exons_to_exclude}
    real_regions_stop_counts = ops.get_region_stop_counts(restricted_real_seqs, window_start, window_end)

    # write to file
    with open(output_file, "w") as outfile:
        # write the header
        outfile.write("sim_no,flank_count,core_count,exons_counted\n")
        # write the real line
        outfile.write("real,{0},{1},{2}\n".format(real_regions_stop_counts[0], real_regions_stop_counts[1], len(restricted_real_seqs)))
        # write the simulants
        for sim in sorted(temp_filelist):
            # get the simulation number
            names, counts = gen.read_fasta(temp_filelist[sim])
            required = {name: counts[i] for i, name in enumerate(names) if name not in exons_to_exclude}
            flank_count = sum([int(required[name].split(",")[0]) for name in required])
            core_count = sum([int(required[name].split(",")[1]) for name in required])
            outfile.write("{0},{1},{2},{3}\n".format("sim_{0}".format(sim), flank_count, core_count, len(required)))

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


def sim_stop_count_mm(seq_fasta, required_simulations, output_file, parallel=True, seeds=None, seq_seeds=None):
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

    # names, seqs_list = names[:100], seqs_list[:100]

    # return a dictionary and list of unique sequences
    seqs = seqo.get_unique_seqs(names, seqs_list)



    # get the gc contents of the sequences
    gc_contents = {}
    for seq in seqs:
        gc_contents[seq] = seqo.calc_seq_gc(seqs[seq])

    # get the stop codon counts for the real sequences
    stop_counts = seqo.get_stop_counts(seqs)

    # create a temporary output directory
    temp_dir = "temp_data_stop_count_sim"
    gen.create_output_directories(temp_dir)

    # run the simulations
    sim_args = [seqs, temp_dir, seeds, seq_seeds]
    # run the simulation
    outputs = run_simulation_function(required_simulations, sim_args, simo.sim_stop_counts_mm, parallel=parallel)


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
def sim_motifs(motifs_file, output_file, required_simulations, seeds=None, seq_seeds=None, parallel=True):

    # get a list of motifs
    motifs = [i[0] for i in gen.read_many_fields(motifs_file, "\t")]
    # get the dicnucleotide and nucleotide content of the sequences
    dinucleotide_content = seqo.get_dinucleotide_content(motifs)
    nucleotide_content = seqo.get_nucleotide_content(motifs)

    # get the stop codon counts for the real sequences
    real_stop_counts = sum(seqo.get_stop_counts(motifs))

    # create a temporary output directory
    temp_dir = "temp_data_motif_sim"
    gen.create_output_directories(temp_dir)

    # run the simulations
    simulations = list(range(required_simulations))
    sim_args = [motifs, dinucleotide_content, nucleotide_content, temp_dir, seeds, seq_seeds]

    outputs = run_simulation_function(required_simulations, sim_args, simo.sim_motif_stop_counts, parallel=parallel)

    # get temp filelist so we can order simulants for tests
    filelist = []
    for output in outputs:
        filelist.extend(output)
    temp_filelist = fo.order_temp_files(outputs)

    with open(output_file, "w") as outfile:
        outfile.write("sim_no,stop_count\n")
        outfile.write("real,{0}\n".format(real_stop_counts))
        for file in temp_filelist:
            data = gen.read_fasta(temp_filelist[file])[1]
            outfile.write("{0},{1}\n".format(file, data[0]))
    gen.remove_directory(temp_dir)



def sim_stop_density(seqs_fasta, non_features_fasta, threshold, required_simulations, output_directory, output_file, parallel=True):

    matched_seqs_directory = "temp_matched_gc_dir"
    gen.create_output_directories(matched_seqs_directory)

    simulations = list(range(required_simulations))
    sim_args = [seqs_fasta, non_features_fasta, 0.5, matched_seqs_directory]
    outputs = run_simulation_function(required_simulations, sim_args, simo.generate_matched_gc_seqs, parallel=parallel)

    # get temp filelist so we can order simulants for tests
    filelist = []
    for output in outputs:
        filelist.extend(output)
    temp_filelist = fo.order_temp_files(outputs)

    seqs = {name: gen.read_fasta(seqs_fasta)[1][i] for i, name in enumerate(gen.read_fasta(seqs_fasta)[0])}

    with open(output_file, "w") as outfile:
        outfile.write("id,gc,real,{0}\n".format(",".join(["sim_{0}".format(i) for i in sorted(temp_filelist)])))
        real_gc = seqo.calc_seqs_gc(seqs)
        real_stop_densities = {name: seqo.calc_stop_density(seqs[name]) for name in seqs}

        sim_densities_list = collections.defaultdict(lambda: [])
        for file in sorted(temp_filelist):
            sim_densities = {name: seqo.calc_stop_density(gen.read_fasta(temp_filelist[file])[1][i]) for i, name in enumerate(gen.read_fasta(temp_filelist[file])[0])}
            for name in sim_densities:
                sim_densities_list[name].append(sim_densities[name])

        for name in sorted(real_gc):
            outfile.write("{0},{1},{2},{3}\n".format(name, real_gc[name], real_stop_densities[name], ",".join(gen.stringify(sim_densities_list[name]))))
