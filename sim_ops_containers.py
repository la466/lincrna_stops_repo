import generic as gen
import sim_ops as simo
import seq_ops as seqo
import sequence_ops as sequo
import file_ops as fo
import ops_containers as opsc
import ops
import numpy as np
import pandas as pd
import os
import collections
import time
import csv
from useful_motif_sets import dinucleotides, nucleotides, stops
from itertools import zip_longest


def get_temp_filelist(outputs):
    filelist = []
    for output in outputs:
        filelist.extend(output)
    temp_filelist = fo.order_temp_files(outputs)
    return temp_filelist

def run_simulation_function(required_simulations, sim_args, function_to_run, kwargs_dict = None, parallel = True, workers=None, sim_run=True):
    """
    Wrapper to run simulation function

    Args:
        required_simulations (int): the number of times to run simulation
        sim_args (list): list of arguments for the simulation function
        function_to_run (func): simulation_function
        parallel (bool): if true, run in parallel
        sim_run (bool): if true, running a simulation, else running over another list

    Returns:
        outputs (list): list containing the result of the simulation
    """

    if sim_run:
        # get a list of simulations to iterate over
        simulations = list(range(required_simulations))
    else:
        simulations = required_simulations
    # run the simulations
    if parallel:
        # add foo to argument list for parallelisation
        sim_args.insert(0, "foo")
        if not workers:
            workers = os.cpu_count() - 2
        processes = gen.run_in_parallel(simulations, sim_args, function_to_run, kwargs_dict = kwargs_dict, workers = workers)
        results = []
        for process in processes:
            results.append(process.get())
        # now merge into one
        if isinstance(results[0],list):
            flattened_outputs = []
            [flattened_outputs.extend(i) for i in results]
            outputs = flattened_outputs
        elif isinstance(results[0], dict):
            keys = list(results[0].keys())
            if len(keys):
                if isinstance(results[0][keys[0]], list):
                    flattened_outputs = collections.defaultdict(lambda: [])
                    for result in results:
                        for key in result:
                            flattened_outputs[key].extend(result[key])
                    #unpickle
                    outputs = {i: flattened_outputs[i] for i in flattened_outputs}
                else:
                    flattened_outputs = {}
                    [flattened_outputs.update(i) for i in results]
                    outputs = flattened_outputs
            else:
                outputs = {}

    else:
        if kwargs_dict:
            outputs = function_to_run(simulations, *sim_args, **kwargs_dict)
        else:
            outputs = function_to_run(simulations, *sim_args)

    return outputs


def generate_gc_controls(input_file, untranslated_sequence_file, output_directory, required_simulations = None):

    if not required_simulations:
        required_simulations = 1000

    gen.create_output_directories(output_directory)
    # get a list of the sequences
    sequence_names, sequences = gen.read_fasta(input_file)
    sequence_names = [i.split("(")[0] for i in sequence_names]
    sequence_list = {name: sequences[i] for i, name in enumerate(sequence_names)}
    # get the untransribed regions and join
    untransribed_sequences = [i for i in gen.read_fasta(untranslated_sequence_file)[1] if "N" not in i]
    untranscribed_sequence = "".join(untransribed_sequences)
    # set up to generate sequences
    args = [sequence_list, untranscribed_sequence, output_directory, required_simulations]
    run_simulation_function(sequence_names, args, simo.get_gc_matched_seqs, sim_run = False)




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
    sim_args = [seqs, temp_dir, seeds, seq_seeds]

    outputs = run_simulation_function(simulations, sim_args, simo.sim_orf_lengths, sim_run = False)
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
    motifs = [i[0] for i in gen.read_many_fields(motifs_file, "\t") if "#" not in i[0]]
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


def simulate_sequence_stop_density(input_fasta, output_file, required_simulations, families_file = None):
    """
    Given a fasta file, simulate the sequences through dinucleotide matched controls
    and calculate the stop codon densities.

    Args:
        input_fasta (str): path to fasta file
        output_file (str): path to output file
        required_simulations (int): number of simulations to run
        families_file (str): if set, path to families file. Take the mean of sequences grouped into family
    """

    # create the output directory if not already exists
    gen.create_output_directories("/".join(output_file.split("/")[:-1]))
    # create temp simulation directory
    temp_dir = "temp_sims"
    gen.create_output_directories(temp_dir)


    names, seqs = gen.read_fasta(input_fasta)
    # get the dicnucleotide and nucleotide content of the sequences
    dinucleotide_content = seqo.get_dinucleotide_content(seqs)
    dinucleotide_probabilities = [dinucleotide_content[i] for i in sorted(dinucleotide_content)]
    nucleotide_content = seqo.get_nucleotide_content(seqs)
    nucleotide_probabilities = [nucleotide_content[i] for i in sorted(nucleotide_content)]

    query_sequences = {name: seqs[i] for i, name in enumerate(names)}

    # run simulation
    args = [query_sequences, nucleotide_probabilities, dinucleotide_probabilities, temp_dir]
    outputs = run_simulation_function(required_simulations, args, simo.simulate_sequence_stop_density)
    # calculate the real density
    real_densities = {name: seqo.calc_stop_density(seqs[i]) for i, name in enumerate(sorted(names))}
    gcs = {name: seqo.calc_seq_gc(seqs[i]) for i, name in enumerate(sorted(names))}

    if families_file:
        family_list = gen.read_many_fields(families_file, "\t")

    with open(output_file, "w") as outfile:
        # get families and group if required
        if family_list:
            gcs = sequo.group_family_results(gcs, family_list)
            gcs = {i: np.mean(gcs[i]) for i in gcs}
            real_densities = sequo.group_family_results(real_densities, family_list)
            real_densities = {i: np.mean(real_densities[i]) for i in real_densities}

        outfile.write("sim_id,{0}\n".format(",".join([i for i in sorted(real_densities)])))
        outfile.write("gc,{0}\n".format(",".join(gen.stringify([gcs[i] for i in sorted(gcs)]))))
        outfile.write("real,{0}\n".format(",".join(gen.stringify([real_densities[i] for i in sorted(real_densities)]))))


        for i, file in enumerate(outputs):
            sim_densities = {i[0]: float(i[1]) for i in gen.read_many_fields(file, "\t")}
            if family_list:
                sim_densities = sequo.group_family_results(sim_densities, family_list)
                sim_densities = {i: np.mean(sim_densities[i]) for i in real_densities}
            outfile.write("sim_{0},{1}\n".format(i+1, ",".join(gen.stringify([sim_densities[i] for i in sorted(sim_densities)]))))

    gen.remove_directory(temp_dir)

    output_results_file = "{0}.nds.csv".format(".".join(output_file.split(".")[:-1]))
    # get the results
    results = pd.read_csv(output_file)
    cols = list(results)[1:]
    with open(output_results_file, "w") as outfile:
        outfile.write("id,gc,real_density,nd\n")
        for col in cols:
            gc = results[col][0]
            real = results[col][1]
            sims = results[col][2:]
            nd = np.divide(real - np.mean(sims), np.mean(sims))
            outfile.write("{0},{1},{2},{3}\n".format(col, gc, real, nd))



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


def sim_cds_stop_density(source_file, required_simulations, output_file, parallel=True):
    """
    Given a file containing sequences, simulate them by shuffling codons and calculate
    stop codon densitiesself.

    Args:
        source_file (str): path to fasta file containing sequences
        required_simulations (int): number of simulations to run
        output_file (str): path to output file
        parallel (bool): if true, run in parallel
    """

    start_time = time.time()

    temp_output_directory = "temp_cds_dir"
    gen.create_output_directories(temp_output_directory)

    cds_names, cds_seqs = gen.read_fasta(source_file)


    simulations = list(range(required_simulations))
    sim_args = [cds_seqs, temp_output_directory]
    outputs = run_simulation_function(required_simulations, sim_args, simo.sim_cds_seqs_stop_counts, parallel=parallel)

    real_stop_counts = seqo.get_stop_counts([seq[:-3] for seq in cds_seqs])
    real_densities = [np.divide(count, len(cds_seqs[i])) for i, count in enumerate(real_stop_counts)]
    # get temp filelist so we can order simulants for tests
    temp_filelist = get_temp_filelist(outputs)

    real_gc = [seqo.calc_seq_gc(seq[:-3]) for seq in cds_seqs]

    with open(output_file, "w") as outfile:
        outfile.write("sim_no,{0}\n".format(",".join([name for name in cds_names[:len(cds_seqs)]])))
        outfile.write("gc,{0}\n".format(",".join([gc for gc in gen.stringify(real_gc[:len(cds_seqs)])])))
        outfile.write("real,{0}\n".format(",".join(gen.stringify(real_densities))))
        for file in temp_filelist:
            data = gen.read_many_fields(temp_filelist[file], ",")
            outfile.write("{0},{1}\n".format(file, ",".join(data[0])))

    gen.remove_directory(temp_output_directory)
    gen.get_time(start_time)


def sim_exon_stop_density(source_cds_fasta, source_cds_bed, source_coding_exons_bed, required_simulations, output_directory, output_file, parallel=True):
    """
    Given a file containing sequences, simulate them by shuffling codons and calculate
    stop codon densitiesself.

    Args:
        source_file (str): path to fasta file containing sequences
        source_file (str): path to fasta file containing sequences
        source_file (str): path to fasta file containing sequences
        required_simulations (int): number of simulations to run
        output_dir (str): path to output directory
        output_file (str): path to output file
        parallel (bool): if true, run in parallel
    """

    start_time = time.time()

    temp_output_directory = "temp_exons_dir"
    gen.create_output_directories(temp_output_directory)

    exon_positions_bed = "{0}/exon_positions.bed".format(output_directory)
    if not os.path.isfile(exon_positions_bed):
        seqo.get_exon_positions_bed(source_cds_bed, source_coding_exons_bed, exon_positions_bed)

    cds_names, cds_seqs = gen.read_fasta(source_cds_fasta)
    cds_list = {name: cds_seqs[i] for i, name in enumerate(cds_names[:1500])}

    exon_info = simo.get_exon_info(exon_positions_bed)
    real_densities, real_core_densities = simo.get_exon_stop_densities(cds_list, exon_info)


    simulations = list(range(required_simulations))
    sim_args = [cds_list, exon_positions_bed, temp_output_directory]
    outputs = run_simulation_function(required_simulations, sim_args, simo.sim_exon_stop_density, parallel=parallel)

    # get temp filelist so we can order simulants for tests
    temp_filelist = get_temp_filelist(outputs)

    sim_densities = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: [])))
    sim_core_densities = collections.defaultdict(lambda: [])

    for file in temp_filelist:
        data = gen.read_many_fields(temp_filelist[file], ",")
        for line in data[:-1]:
            end = int(line[0])
            region = int(line[1])
            sim_density_list = [float(i) for i in line[2:]]
            for i, sim_density in enumerate(sim_density_list):
                sim_densities[end][region][i].append(sim_density)
        sim_cores = [float(i) for i in data[-1][1:]]
        for i, sim_core_density in enumerate(sim_cores):
            sim_core_densities[i].append(sim_core_density)

    with open(output_file, "w") as outfile:

        outlist = []
        for end in real_densities:
            for region in real_densities[end]:
                outdata = ["{0}_{1}".format(end, region)]
                for i, real_density in enumerate(real_densities[end][region]):
                    nd = np.divide(real_density - np.mean(sim_densities[end][region][i]), np.mean(sim_densities[end][region][i]))
                    outdata.append(nd)
                outlist.append(outdata)
        core = ["core"]
        for i, core_density in enumerate(real_core_densities):
            nd = np.divide(core_density - np.mean(sim_core_densities[i]), np.mean(sim_core_densities[i]))
            core.append(nd)
        outlist.append(core)
        outrows = zip_longest(*outlist, fillvalue = '')
        for item in outrows:
            item = [i for i in item]
            outfile.write("{0}\n".format(",".join(gen.stringify(item))))

    gen.remove_directory(temp_output_directory)
    gen.get_time(start_time)


def sim_intron_density(coding_exons_bed, genome_fasta, required_simulations, output_directory, output_file, parallel=True):

    start_time = time.time()

    intron_bed = "{0}/internal_introns.bed".format(output_directory)
    intron_fasta = "{0}/internal_introns.fasta".format(output_directory)
    if not os.path.isfile(intron_bed):
        print("Getting intron coordinates...")
        ops.get_introns_from_bed(coding_exons_bed, intron_bed)
        fo.fasta_from_intervals(intron_bed, intron_fasta, genome_fasta, names = True)


    temp_output_directory = "temp_intron_dir"
    gen.create_output_directories(temp_output_directory)

    intron_seq_names, intron_seqs = gen.read_fasta(intron_fasta)
    intron_seqs = intron_seqs[:10000]

    simulations = list(range(required_simulations))
    sim_args = [intron_seqs, temp_output_directory]
    outputs = run_simulation_function(required_simulations, sim_args, simo.sim_intron_seqs_stop_density, parallel=parallel)


    temp_filelist = get_temp_filelist(outputs)

    real_stop_counts = seqo.get_stop_counts(intron_seqs)
    real_densities = [np.divide(count, len(intron_seqs[i])) for i, count in enumerate(real_stop_counts)]

    real_gc = [seqo.calc_seq_gc(seq) for seq in intron_seqs]

    with open(output_file, "w") as outfile:
        outfile.write("sim_no,{0}\n".format(",".join([name for name in intron_seq_names[:len(intron_seqs)]])))
        outfile.write("gc,{0}\n".format(",".join([gc for gc in gen.stringify(real_gc[:len(intron_seqs)])])))
        outfile.write("real,{0}\n".format(",".join(gen.stringify(real_densities))))
        for file in temp_filelist:
            data = gen.read_many_fields(temp_filelist[file], ",")
            outfile.write("{0},{1}\n".format(file, ",".join(data[0])))

    gen.remove_directory(temp_output_directory)
    gen.get_time(start_time)


def generate_exon_dinucleotide_controls(motifs_file, exons_fasta, required_simulations, output_directory, match_density = None):

    exon_seqs = gen.read_fasta(exons_fasta)[1]
    exon_seqs = [i for i in exon_seqs if len(i) >= 138]

    exon_flanks = []
    for seq in exon_seqs:
        exon_flanks.append(seq[:69])
        exon_flanks.append(seq[-69:])

    dinucleotide_content = seqo.get_dinucleotide_content(exon_flanks)
    nucleotide_content = seqo.get_nucleotide_content(exon_flanks)

    if os.path.isdir(output_directory):
        gen.remove_directory(output_directory)
    gen.create_output_directories(output_directory)

    start = time.time()

    motifs = sequo.read_motifs(motifs_file)

    # create sets of motifs
    print("Simluating sequences...")
    kwargs = {"match_density": match_density}
    args = [motifs, dinucleotide_content, nucleotide_content, output_directory]
    simulated_seqs = run_simulation_function(required_simulations, args, simo.generate_dinucleotide_matched_seqs, kwargs_dict = kwargs)

def generate_motif_dinucleotide_controls(motifs_file, required_simulations, output_directory, match_density = None, match_subs = None):

    motifs = [i[0] for i in gen.read_many_fields(motifs_file, "\t") if "#" not in i[0]]
    dinucleotide_content = seqo.get_dinucleotide_content(motifs)
    nucleotide_content = seqo.get_nucleotide_content(motifs)

    if os.path.isdir(output_directory):
        gen.remove_directory(output_directory)
    gen.create_output_directories(output_directory)

    start = time.time()

    # create sets of motifs
    print("Simluating sequences...")
    kwargs = {"match_density": match_density, "match_subs": match_subs}
    args = [motifs, dinucleotide_content, nucleotide_content, output_directory]
    simulated_seqs = run_simulation_function(required_simulations, args, simo.generate_dinucleotide_matched_seqs, kwargs_dict = kwargs)



def sim_motif_codon_densities(seqs_file, required_simulations, output_directory):

    # get gc matched motifs
    stops = ["TAA", "TAG", "TGA"]
    gc_matchd_motifs_file = "{0}/gc_matched_combinations.bed".format(output_directory)
    if not os.path.isfile(gc_matchd_motifs_file):
        seqo.get_gc_matched_motifs(stops, gc_matchd_motifs_file)

    seqs = [i[0] for i in gen.read_many_fields(seqs_file, "\t") if "#" not in i[0]]
    dinucleotide_content = seqo.get_dinucleotide_content(seqs)
    nucleotide_content = seqo.get_nucleotide_content(seqs)

    query_sets = gen.read_many_fields(gc_matchd_motifs_file, "\t")
    query_sets = [stops] + query_sets

    temp_dir = "temp_sim_motif_codon_densities"
    temp_seq_dir = "temp_seq_dir"
    gen.create_output_directories(temp_dir)
    gen.create_output_directories(temp_seq_dir)

    start = time.time()

    # create sets of motifs
    print("Simluating sequences...")
    args = [seqs, dinucleotide_content, nucleotide_content, temp_seq_dir]
    simulated_seqs = run_simulation_function(required_simulations, args, simo.generate_dinucleotide_matched_seqs)

    # now get densities for each set
    density_args = [seqs, simulated_seqs, temp_dir]
    density_files = run_simulation_function(query_sets, density_args, opsc.get_sequence_densities, parallel=True, sim_run=False)

    output_file = "{0}/{1}_sim_motif_densities.csv".format(output_directory, seqs_file.split("/")[-1].split(".")[0])

    with open(output_file, "w") as outfile:
        outfile.write("id,codons,density\n")
        for i,file in enumerate(sorted(density_files)):
            ids, densities = gen.read_fasta(file)
            real_density = float(densities[0])
            sim_densities = densities[1]
            sim_densities = [float(i) for i in densities[1].split(",")]
            nd = np.divide(real_density - np.mean(sim_densities), np.mean(sim_densities))
            outfile.write("{0},{1},{2}\n".format(i+1, ids[0], nd))
    #
    gen.get_time(start)

    gen.remove_directory(temp_dir)
    gen.remove_directory(temp_seq_dir)


def sim_motif_codon_densities_genome(genome_fasta, gtf, seqs_file, required_simulations, output_directory):

    # get gc matched motifs
    stops = ["TAA", "TAG", "TGA"]
    gc_matchd_motifs_file = "{0}/gc_matched_combinations.bed".format(output_directory)
    if not os.path.isfile(gc_matchd_motifs_file):
        seqo.get_gc_matched_motifs(stops, gc_matchd_motifs_file)

    seqs = [i[0] for i in gen.read_many_fields(seqs_file, "\t") if "#" not in i[0]]
    # now write to fasta for getting gc matched sequences
    gen.create_output_directories("temp_data")
    temp_seqs_fasta = "temp_data/temp_seqs.fasta"
    with open(temp_seqs_fasta, "w") as temp_seqs_file:
        for i, seq in enumerate(seqs):
            temp_seqs_file.write(">{0}\n{1}\n".format(i+1, seq))
    # dinucleotide_content = seqo.get_dinucleotide_content(seqs)
    # nucleotide_content = seqo.get_nucleotide_content(seqs)

    query_sets = gen.read_many_fields(gc_matchd_motifs_file, "\t")
    query_sets = [stops] + query_sets

    temp_dir = "temp_sim_motif_codon_densities"
    temp_seq_dir = "temp_gc_matched_seq_dir"
    gen.create_output_directories(temp_dir)
    gen.create_output_directories(temp_seq_dir)

    start = time.time()

    # create sets of motifs
    print("Simluating sequences...")
    genome_seq_outputs = "{0}/genome_sequence_files".format(output_directory)
    gen.create_output_directories(genome_seq_outputs)

    # get the sequences for non features
    features_bed = "{0}/genome_features.bed".format(genome_seq_outputs)
    non_features_bed = "{0}/non_genome_features.bed".format(genome_seq_outputs)
    non_features_fasta = "{0}/non_genome_features.fasta".format(genome_seq_outputs)
    if not os.path.isfile(non_features_fasta):
        seqo.get_non_transcribed_regions(gtf_file, genome_fasta, features_bed, non_features_bed, non_features_fasta, genome_seq_outputs)

    threshold = 0.05
    names, non_features_seqs = gen.read_fasta(non_features_fasta)
    non_features_seqs= non_features_seqs[:100]
    args = [temp_seqs_fasta, non_features_seqs, threshold, temp_seq_dir]
    simulated_seqs = run_simulation_function(required_simulations, args, simo.generate_matched_gc_seqs)

    # now get densities for each set
    density_args = [seqs, simulated_seqs, temp_dir]
    density_files = run_simulation_function(query_sets, density_args, opsc.get_sequence_densities, parallel=True, sim_run=False)

    output_file = "{0}/{1}_sim_motif_genome_gc_matched_densities.csv".format(output_directory, seqs_file.split("/")[-1].split(".")[0])

    with open(output_file, "w") as outfile:
        outfile.write("id,codons,density\n")
        for i,file in enumerate(sorted(density_files)):
            ids, densities = gen.read_fasta(file)
            real_density = float(densities[0])
            sim_densities = densities[1]
            sim_densities = [float(i) for i in densities[1].split(",")]
            nd = np.divide(real_density - np.mean(sim_densities), np.mean(sim_densities))
            outfile.write("{0},{1},{2}\n".format(i+1, ids[0], nd))
    # #
    gen.get_time(start)
    #
    # gen.remove_file(temp_seqs_fasta)
    gen.remove_directory(temp_dir)
    # gen.remove_directory(temp_seq_dir)


def generate_dint_controls(input_fasta, output_directory, simulations = None):

    gen.create_output_directories(output_directory)
    names, seqs = gen.read_fasta(input_fasta)
    # names = names[:10]

    seq_list = {name: seqs[i] for i, name in enumerate(names)}

    dinucleotide_content = seqo.get_dinucleotide_content([seq_list[i] for i in seq_list])
    nucleotide_content = seqo.get_nucleotide_content([seq_list[i] for i in seq_list])

    args = [seq_list, dinucleotide_content, nucleotide_content, output_directory]
    kwargs_dict = {"simulations": simulations}
    run_simulation_function(names, args, simo.generate_dint_controls,  kwargs_dict = kwargs_dict, sim_run = False)


def generate_dint_intron_controls(input_fasta, output_directory):

    # need to concatenate introns to make it easier to simulate
    intron_names, intron_seqs = gen.read_fasta(input_fasta)
    intron_list = collections.defaultdict(lambda: [])
    for i, name in enumerate(intron_names):
        intron_list[name.split(".")[0]].append(intron_seqs[i])

    temp_dir = "temp_files"
    gen.create_output_directories(temp_dir)
    temp_fasta = "{0}/temp_introns.fasta".format(temp_dir)
    with open(temp_fasta, "w") as temp:
        [temp.write(">{0}\n{1}\n".format(i, "".join(intron_list[i]))) for i in intron_list]

    generate_dint_controls(temp_fasta, output_directory)
    gen.remove_file(temp_fasta)

def calculate_substitution_rates(iteration_list, motifs, seq_list, sim_randomisations, output_directory):
    """
    For a set of motifs, predict hits to the sequence and calcualte the rate of
    substitution between stop nucleotides and non stop nucleotides both in the
    hits and outside of the hits.

    Args:
        iteration_list (list): list of items to iterate over
        motifs (list): list of motifs to predict hits to
        seq_list (dict): set of alignment sequences
        sim_randomisations (dict): set defining simulation or not
        output_directory (str): path to output directory

    Returns:
        outputs (list): list of outputs
    """

    outputs = []
    # check if the iteration exists
    if iteration_list:
        np.random.seed()
        for iter, iteration in enumerate(iteration_list):
            gen.print_parallel_status(iter, iteration_list)

            all_hit_human_motifs = []
            all_hit_mac_motifs = []
            all_non_hit_human_motifs = []
            all_non_hit_mac_motifs = []
            all_hit_stop_human_motifs = []
            all_hit_stop_mac_motifs = []
            all_hit_non_stop_human_motifs = []
            all_hit_non_stop_mac_motifs = []
            all_non_hit_stop_human_motifs = []
            all_non_hit_stop_mac_motifs = []
            all_non_hit_non_stop_human_motifs = []
            all_non_hit_non_stop_mac_motifs = []

            # for each sequence id
            for id in seq_list:
                sequences = seq_list[id][0]
                human = sequences[0]
                mac = sequences[1]
                # get the overlaps and chunk into motif overlaps
                hits = sequo.sequence_overlap_indicies(human, motifs)
                chunked_hits = sequo.chunk_overlaps(hits)
                # if a simulation, pick random overlaps
                if sim_randomisations[iteration]:
                    hits = simo.get_random_overlaps(human, chunked_hits)
                    chunked_hits = sequo.chunk_overlaps(hits)

                # now get all the non hit indicies and chunk together
                non_hits = [i for i in range(len(human)) if i not in hits]
                chunked_non_hits = sequo.chunk_overlaps(non_hits)

                # get all the motif parts of the sequences
                hit_human_motifs = ["".join([human[i] for i in chunk]) for chunk in chunked_hits]
                hit_mac_motifs = ["".join([mac[i] for i in chunk]) for chunk in chunked_hits]
                non_hit_human_motifs = ["".join([human[i] for i in chunk]) for chunk in chunked_non_hits]
                non_hit_mac_motifs = ["".join([mac[i] for i in chunk]) for chunk in chunked_non_hits]

                # now get all the parts of the motif hits that also hit stops
                hit_stop_hits = [sequo.sequence_overlap_indicies(motif, stops) for motif in hit_human_motifs]
                all_hit_stop_human_motifs.append("".join(["".join([hit_human_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(hit_stop_hits)]))
                all_hit_stop_mac_motifs.append("".join(["".join([hit_mac_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(hit_stop_hits)]))
                # now get all the parts of the motif hits that dont hit stops
                hit_non_stop_hits = [[i for i in list(range(len(motif))) if i not in hit_stop_hits[no]] for no, motif in enumerate(hit_human_motifs)]
                all_hit_non_stop_human_motifs.append("".join(["".join([hit_human_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(hit_non_stop_hits)]))
                all_hit_non_stop_mac_motifs.append("".join(["".join([hit_mac_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(hit_non_stop_hits)]))
                # get all the parts of non motif hits that hit stops
                non_hit_stop_hits = [sequo.sequence_overlap_indicies(motif, stops) for motif in non_hit_human_motifs]
                all_non_hit_stop_human_motifs.append("".join(["".join([non_hit_human_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(non_hit_stop_hits)]))
                all_non_hit_stop_mac_motifs.append("".join(["".join([non_hit_mac_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(non_hit_stop_hits)]))
                # get all the parts of the non motif hits that dont hit stops
                non_hit_non_stop_hits = [[i for i in list(range(len(motif))) if i not in non_hit_stop_hits[no]] for no, motif in enumerate(non_hit_human_motifs)]
                all_non_hit_non_stop_human_motifs.append("".join(["".join([non_hit_human_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(non_hit_non_stop_hits)]))
                all_non_hit_non_stop_mac_motifs.append("".join(["".join([non_hit_mac_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(non_hit_non_stop_hits)]))
                # concatenate the motif hits and append
                all_hit_human_motifs.append("".join(hit_human_motifs))
                all_hit_mac_motifs.append("".join(hit_mac_motifs))
                all_non_hit_human_motifs.append("".join(non_hit_human_motifs))
                all_non_hit_mac_motifs.append("".join(non_hit_mac_motifs))

            # now concatenate all of the various parts
            all_hit_human_motifs = "".join(all_hit_human_motifs)
            all_hit_mac_motifs = "".join(all_hit_mac_motifs)
            all_non_hit_human_motifs = "".join(all_non_hit_human_motifs)
            all_non_hit_mac_motifs = "".join(all_non_hit_mac_motifs)
            all_hit_stop_human_motifs = "".join(all_hit_stop_human_motifs)
            all_hit_stop_mac_motifs = "".join(all_hit_stop_mac_motifs)
            all_hit_non_stop_human_motifs = "".join(all_hit_non_stop_human_motifs)
            all_hit_non_stop_mac_motifs = "".join(all_hit_non_stop_mac_motifs)
            all_non_hit_stop_human_motifs = "".join(all_non_hit_stop_human_motifs)
            all_non_hit_stop_mac_motifs = "".join(all_non_hit_stop_mac_motifs)
            all_non_hit_non_stop_human_motifs = "".join(all_non_hit_non_stop_human_motifs)
            all_non_hit_non_stop_mac_motifs = "".join(all_non_hit_non_stop_mac_motifs)

            # get the substitution rate for all of the parts
            hit_rate = sequo.get_sub_rate(all_hit_human_motifs, all_hit_mac_motifs)
            non_hit_rate = sequo.get_sub_rate(all_non_hit_human_motifs, all_non_hit_mac_motifs)
            hit_stop_rate = sequo.get_sub_rate(all_hit_stop_human_motifs, all_hit_stop_mac_motifs)
            hit_non_stop_rate = sequo.get_sub_rate(all_hit_non_stop_human_motifs, all_hit_non_stop_mac_motifs)
            non_hit_stop_rate = sequo.get_sub_rate(all_non_hit_stop_human_motifs, all_non_hit_stop_mac_motifs)
            non_hit_non_stop_rate = sequo.get_sub_rate(all_non_hit_non_stop_human_motifs, all_non_hit_non_stop_mac_motifs)

            relative_hit_rate = np.divide(hit_stop_rate, hit_non_stop_rate)
            relative_non_hit_rate = np.divide(non_hit_stop_rate, non_hit_non_stop_rate)
            realtive_diff = relative_hit_rate - relative_non_hit_rate

            local_output = [iteration if iteration == "real" else "sim_{0}".format(iteration), hit_rate, non_hit_rate, "", hit_stop_rate, hit_non_stop_rate, "", non_hit_stop_rate, non_hit_non_stop_rate, "", relative_hit_rate, relative_non_hit_rate, realtive_diff]
            output_file = "{0}/{1}.txt".format(output_directory, iteration)
            with open(output_file, "w") as outfile:
                outfile.write("{0}\n".format(",".join(gen.stringify(local_output))))
            outputs.append(output_file)

    return outputs

def calculate_substitution_rates_motifs(iteration_list, motif_file, motif_files, seq_list, output_directory):
    """
    For a set of motifs, predict hits to the sequence and calcualte the rate of
    substitution between stop nucleotides and non stop nucleotides both in the
    hits and outside of the hits.

    Args:
        iteration_list (list): list of items to iterate over
        motifs (list): list of motifs to predict hits to
        seq_list (dict): set of alignment sequences
        sim_randomisations (dict): set defining simulation or not
        output_directory (str): path to output directory

    Returns:
        outputs (list): list of outputs
    """

    outputs = []
    # check if the iteration exists
    if iteration_list:
        np.random.seed()
        for iter, iteration in enumerate(iteration_list):
            gen.print_parallel_status(iter, iteration_list)

            motifs = sequo.read_motifs(motif_files[iteration])
            motif_lengths = list(set([len(i) for i in motifs]))
            # motifs = sequo.read_motifs(motif_file)

            all_hit_human_motifs = []
            all_hit_mac_motifs = []
            all_non_hit_human_motifs = []
            all_non_hit_mac_motifs = []
            all_hit_stop_human_motifs = []
            all_hit_stop_mac_motifs = []
            all_hit_non_stop_human_motifs = []
            all_hit_non_stop_mac_motifs = []
            all_non_hit_stop_human_motifs = []
            all_non_hit_stop_mac_motifs = []
            all_non_hit_non_stop_human_motifs = []
            all_non_hit_non_stop_mac_motifs = []

            all_sequence = []
            total_hits = []

            # for each sequence id
            for i, id in enumerate(seq_list):
                if i<50:
                    sequences = seq_list[id][0]
                    human = sequences[0]
                    mac = sequences[1]
                    # get the overlaps and chunk into motif overlaps
                    hits = sequo.sequence_overlap_indicies(human, motifs)
                    chunked_hits = sequo.chunk_overlaps(hits)

                    #all hits
                    all_hits = hits
                    # all_hits = sequo.sequence_overlap_indicies(human, motifs)
                    chunked_hits = sequo.chunk_overlaps(all_hits)

                    hit_human_motifs = ["".join([human[i] for i in chunk]) for chunk in chunked_hits]
                    hit_mac_motifs = ["".join([mac[i] for i in chunk]) for chunk in chunked_hits]

                    retained_hits = all_hits

                    # retained_hits = []
                    # for i, human_motif in enumerate(hit_human_motifs):
                    #     mac_motif = hit_mac_motifs[i]
                    #     chunk = chunked_hits[i]
                    #
                    #     # the motifs are both in the set, or the motifs are the same
                    #     # in the case wherer ESEs overlap
                    #     if human_motif in motifs and mac_motif in motifs or human_motif == mac_motif:
                    #         retained_hits.extend(chunk)
                    #     # otherwise elimimate sites that arent
                    #     else:
                    #         motif_retained = []
                    #         for length in motif_lengths:
                    #             for pos in range(0, len(human_motif)-length+1):
                    #                 test_human = human_motif[pos:pos+length]
                    #                 test_mac = mac_motif[pos:pos+length]
                    #                 # print(test_human)
                    #                 # print(test_mac)
                    #                 if test_human in motifs and test_mac in motifs:
                    #                     motif_retained.extend(chunk[pos:pos+length])
                    #         motif_retained = sorted(list(set(motif_retained)))
                    #         retained_hits.extend(motif_retained)
                            # print(motif_lengths)
                            # print(human_motif, mac_motif, chunk)

                    #get the overlaps and chunk into motif overlaps for motifs that are an ese in both cases
                    retained_hits = []
                    for motif in motifs:
                        # predict hits to each motif
                        hits = sequo.sequence_motif_overlap(human, motif)
                        chunked_hits = sequo.chunk_overlaps(hits)
                        # get the motifs
                        hit_human_motifs = ["".join([human[i] for i in chunk]) for chunk in chunked_hits]
                        hit_mac_motifs = ["".join([mac[i] for i in chunk]) for chunk in chunked_hits]
                        # ask whether they are both in the motifs
                        for j, motif in enumerate(hit_human_motifs):
                            if motif in motifs and hit_mac_motifs[j] in motifs:
                                retained_hits.extend(chunked_hits[j])


                    # get a set of retained hits and chunk them
                    retained_hits = sorted(list(set(retained_hits)))

                    # print(all_hits, "\n")
                    # print(retained_hits)
                    # # print(test_retained)
                    # print("\n\n")

                    chunked_retained_hits = sequo.chunk_overlaps(retained_hits)
                    non_hits = [i for i in range(len(human)) if i not in all_hits]
                    chunked_non_hits = sequo.chunk_overlaps(non_hits)


                    all_sequence.append(human)
                    total_hits.append(len(chunked_retained_hits))

                    ## all motifs
                    # get all the parts of the sequences that are/are not a motif in both cases
                    hit_human_motifs = ["".join([human[i] for i in chunk]) for chunk in chunked_retained_hits]
                    hit_mac_motifs = ["".join([mac[i] for i in chunk]) for chunk in chunked_retained_hits]
                    non_hit_human_motifs = ["".join([human[i] for i in chunk]) for chunk in chunked_non_hits]
                    non_hit_mac_motifs = ["".join([mac[i] for i in chunk]) for chunk in chunked_non_hits]


                    # add to global list
                    all_hit_human_motifs.append("".join(hit_human_motifs))
                    all_hit_mac_motifs.append("".join(hit_mac_motifs))
                    all_non_hit_human_motifs.append("".join(non_hit_human_motifs))
                    all_non_hit_mac_motifs.append("".join(non_hit_mac_motifs))

                    ### stops specifically


                    ## within motifs
                    # now get all the parts of the motif hits that also hit stops
                    hit_stop_hits = [sequo.sequence_overlap_indicies(motif, stops) for motif in hit_human_motifs]
                    # get the parts of the motifs that hit stops
                    hit_stop_human_motifs = ["".join([hit_human_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(hit_stop_hits)]
                    hit_stop_mac_motifs = ["".join([hit_mac_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(hit_stop_hits)]
                    # and append to global list
                    all_hit_stop_human_motifs.append("".join(hit_stop_human_motifs))
                    all_hit_stop_mac_motifs.append("".join(hit_stop_mac_motifs))
                    # get the parts of the motifs that dont hit stops
                    # get the positions that arent part of a stop
                    hit_non_stop_hits = [[i for i in range(len(motif)) if i not in hit_stop_hits[no]] for no, motif in enumerate(hit_human_motifs)]
                    # get the parts of the motifs that dont hit stops
                    hit_non_stop_human_motifs = ["".join([hit_human_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(hit_non_stop_hits)]
                    hit_non_stop_mac_motifs = ["".join([hit_mac_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(hit_non_stop_hits)]
                    # and append to the global list
                    all_hit_non_stop_human_motifs.append("".join(hit_non_stop_human_motifs))
                    all_hit_non_stop_mac_motifs.append("".join(hit_non_stop_mac_motifs))


                    ## for non motif hits
                    non_hit_stop_hits = [sequo.sequence_overlap_indicies(motif, stops) for motif in non_hit_human_motifs]
                    # get the parts of the motifs that hit stops
                    non_hit_stop_human_motifs = ["".join([non_hit_human_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(non_hit_stop_hits)]
                    non_hit_stop_mac_motifs = ["".join([non_hit_mac_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(non_hit_stop_hits)]
                    # add to global list
                    all_non_hit_stop_human_motifs.append("".join(non_hit_stop_human_motifs))
                    all_non_hit_stop_mac_motifs.append("".join(non_hit_stop_mac_motifs))
                    # get the parts of the motifs that dont hits stops
                    # get the positions
                    non_hit_non_stop_hits = [[i for i in list(range(len(motif))) if i not in non_hit_stop_hits[no]] for no, motif in enumerate(non_hit_human_motifs)]
                    # get the motif parts
                    non_hit_non_stop_human_motifs = ["".join([non_hit_human_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(non_hit_non_stop_hits)]
                    non_hit_non_stop_mac_motifs = ["".join([non_hit_mac_motifs[i][pos] for pos in chunk]) for i, chunk in enumerate(non_hit_non_stop_hits)]
                    # add to global list
                    all_non_hit_non_stop_human_motifs.append("".join(non_hit_non_stop_human_motifs))
                    all_non_hit_non_stop_mac_motifs.append("".join(non_hit_non_stop_mac_motifs))



            # now concatenate all of the various parts
            all_hit_human_motifs = "".join(all_hit_human_motifs)
            all_hit_mac_motifs = "".join(all_hit_mac_motifs)
            all_non_hit_human_motifs = "".join(all_non_hit_human_motifs)
            all_non_hit_mac_motifs = "".join(all_non_hit_mac_motifs)
            all_hit_stop_human_motifs = "".join(all_hit_stop_human_motifs)
            all_hit_stop_mac_motifs = "".join(all_hit_stop_mac_motifs)
            all_hit_non_stop_human_motifs = "".join(all_hit_non_stop_human_motifs)
            all_hit_non_stop_mac_motifs = "".join(all_hit_non_stop_mac_motifs)
            all_non_hit_stop_human_motifs = "".join(all_non_hit_stop_human_motifs)
            all_non_hit_stop_mac_motifs = "".join(all_non_hit_stop_mac_motifs)
            all_non_hit_non_stop_human_motifs = "".join(all_non_hit_non_stop_human_motifs)
            all_non_hit_non_stop_mac_motifs = "".join(all_non_hit_non_stop_mac_motifs)


            # get the substitution rate for all of the parts
            hit_rate = sequo.get_sub_rate(all_hit_human_motifs, all_hit_mac_motifs, p = True)
            non_hit_rate = sequo.get_sub_rate(all_non_hit_human_motifs, all_non_hit_mac_motifs)
            hit_stop_rate = sequo.get_sub_rate(all_hit_stop_human_motifs, all_hit_stop_mac_motifs)
            hit_non_stop_rate = sequo.get_sub_rate(all_hit_non_stop_human_motifs, all_hit_non_stop_mac_motifs)
            non_hit_stop_rate = sequo.get_sub_rate(all_non_hit_stop_human_motifs, all_non_hit_stop_mac_motifs)
            non_hit_non_stop_rate = sequo.get_sub_rate(all_non_hit_non_stop_human_motifs, all_non_hit_non_stop_mac_motifs)
            #
            relative_hit_rate = np.divide(hit_stop_rate, hit_non_stop_rate)
            relative_non_hit_rate = np.divide(non_hit_stop_rate, non_hit_non_stop_rate)
            relative_diff = np.log(np.divide(relative_hit_rate, relative_non_hit_rate))

            all_sequence = "".join(all_sequence)
            density = np.divide(len(all_hit_human_motifs), len(all_sequence))

            local_output = [iteration if iteration == "real" else "sim_{0}".format(iteration), hit_rate, non_hit_rate,  hit_stop_rate, hit_non_stop_rate,  non_hit_stop_rate, non_hit_non_stop_rate,  relative_hit_rate, relative_non_hit_rate, relative_diff]
            output_file = "{0}/{1}.txt".format(output_directory, iteration)
            with open(output_file, "w") as outfile:
                outfile.write("{0}\n".format(",".join(gen.stringify(local_output))))
            outputs.append(output_file)

    return outputs
