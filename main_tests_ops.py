import generic as gen
import containers as cont
import file_ops as fo
import seq_ops as seqo
import sequence_ops as sequo
import sim_ops_containers as simoc
import conservation as cons
import time
import random
import numpy as np
import collections
import zipfile
import os
import re
import multiprocessing as mp
import pandas as pd
from progressbar import ProgressBar
from useful_motif_sets import stops
import pickle
import scipy.stats

def calc_codon_set_density(exon_list, intron_list, codon_set = None):
    """
    Given a set of codons, calculate their density in exons and introns

    Args:
        exon_list (dict): dict containing exon sequences for each transcript
        intron_list (dict): dict containing intron sequences for each transcript
        codon_set (list): list of codons to test

    Returns:
        exon_densities (dict): density of codon set in exons
        intron_densities (dict): density of codon set in introns
        intron_densities_scaled (dict): density of codon set in introns but scaled to remove reading frame
    """

    if not codon_set:
        print("Please define a set of codons to test...")
        raise Exception

    exon_densities = {name: seqo.calc_seqs_codon_set_density(exon_list[name], codon_set = codon_set) for name in exon_list}
    intron_densities = {name: seqo.calc_seqs_codon_set_density(intron_list[name], codon_set = codon_set) for name in intron_list}
    intron_densities_scaled = {name: seqo.calc_intron_seqs_stop_density(intron_list[name], codon_set = codon_set) for name in intron_list}

    return exon_densities, intron_densities, intron_densities_scaled


def compare_codon_density(exons_fasta, introns_fasta, output_directory, families_file = None):
    """
    Calculate the density of sets of codons with similar GC content to stop codons
    in exons and introns

    Args:
        exons_fasta (str): path to fasta containing exon sequences
        introns_fasta (str): path to fasta containing intron sequences
        output_directory (str): path to output directory
        families_file (str): if set, group results into paralagous_families
    """

    # create the output directory
    gen.create_output_directories(output_directory)


    # get the exons and introns for each transcript
    exon_list, intron_list = get_transcript_exons_and_introns(exons_fasta, introns_fasta)
    # calculate the gc contents of the sets
    exons_gc = {name: seqo.calc_gc_seqs_combined(exon_list[name]) for name in exon_list}
    introns_gc = {name: seqo.calc_gc_seqs_combined(intron_list[name]) for name in intron_list}

    # get the set of codons with the same gc content as stop codons
    gc_matched_motifs_file = "{0}/gc_matched_motifs.bed".format("/".join(output_directory.split("/")[:-1]))
    if not os.path.isfile(gc_matched_motifs_file):
        seqo.get_gc_matched_motifs(stops, gc_matched_motifs_file)
    motif_sets = gen.read_many_fields(gc_matched_motifs_file, "\t")

    args = [len(motif_sets), exons_gc, introns_gc, exon_list, intron_list, output_directory, families_file]
    simoc.run_simulation_function(motif_sets, args, calculate_densities, sim_run = False, workers = int(os.cpu_count()/2))



def calculate_densities(codon_sets, codon_set_count, exons_gc, introns_gc, exons_list, introns_list, output, families_file):
    """
    For each set of codons provided, calculate the density of that set in exons,
    introns and introns when scaled to account for only 2 reading frames.

    Args:
        codon_sets (list): list containing lists of codons
        codon_set_count (int): get global number of codon sets
        exons_gc (dict): dict containig the gc content for each exon set
        introns_gc (dict): dict containig the gc content for each intron set
        exons_list (dict): dict containing the sequences for each exon in a transcript
        introns_list (dict): dict containing the sequences for each intron in a transcript
        output (str): either 1) path to output file or 2) path to output directory
        families_file (str): if set, group results by paralagous families
    """

    if families_file:
        families = gen.read_many_fields(families_file, "\t")
        exons_gc = sequo.group_family_results(exons_gc, families)
        introns_gc = sequo.group_family_results(introns_gc, families)

    for i, codon_set in enumerate(codon_sets):
        print("(W{0}) {1}/{2}: {3}".format(mp.current_process().name.split("-")[-1], i+1, len(codon_sets), "_".join(sorted(codon_set))))

        if codon_set_count > 1:
            output_file = "{0}/{1}.csv".format(output, "_".join(sorted(codon_set)))
        else:
            output_file = output

        # get the densities
        exon_densities, intron_densities, intron_densities_scaled = calc_codon_set_density(exons_list, introns_list, codon_set = codon_set)

        # group by family
        if families_file:
            exon_densities = sequo.group_family_results(exon_densities, families)
            intron_densities = sequo.group_family_results(intron_densities, families)
            intron_densities_scaled = sequo.group_family_results(intron_densities_scaled, families)

        # write to file
        with open(output_file, "w") as outfile:
            outfile.write("id,exon_gc,intron_gc,exon_density,intron_density,intron_density_scaled\n")
            for id in exon_densities:
                args = [id, np.median(exons_gc[id]), np.median(introns_gc[id]), np.median(exon_densities[id]), np.median(intron_densities[id]), np.median(intron_densities_scaled[id])]
                outfile.write("{0}\n".format(",".join(gen.stringify(args))))

    return []

def get_transcript_exons_and_introns(exons_fasta, introns_fasta):
    """
    Given a file containing exon and intron sequences, return two lists
    grouped by transcript of those that for transcripts that only occur in both files

    Args:
        exons_fasta (str): path to exons fasta
        introns_fasta (str): path to introns fasta

    Returns:
        transcript_exon_list (dict): dict containing exon sequences for each transcript
        transcript_intron_list (dict): dict containing intron sequences for each transcript
    """


    print("Getting exons and introns...")
    # get a list of introns
    intron_names, intron_seqs = gen.read_fasta(introns_fasta)
    intron_list = {name.split("(")[0]: intron_seqs[i] for i, name in enumerate(intron_names)}
    # get a list of introns grouped by transcript
    transcript_intron_list = collections.defaultdict(lambda: [])
    [transcript_intron_list[name.split(".")[0]].append(intron_list[name]) for name in intron_list]
    # get a list of exons that have given introns
    exon_names, exon_seqs = gen.read_fasta(exons_fasta)
    exon_list = {name.split("(")[0]: exon_seqs[i] for i, name in enumerate(exon_names) if name.split(".")[0] in transcript_intron_list}
    # get a list of exons grouped by transcript
    transcript_exon_list = collections.defaultdict(lambda: [])
    [transcript_exon_list[name.split(".")[0]].append(exon_list[name]) for name in exon_list]

    # unpickle
    transcript_exon_list = {i: transcript_exon_list[i] for i in transcript_exon_list}
    transcript_intron_list = {i: transcript_intron_list[i] for i in transcript_intron_list}

    return transcript_exon_list, transcript_intron_list

def compare_stop_density(exons_fasta, introns_fasta, output_file, families_file = None):
    """
    Calculate the stop codon density in exons and introns

    Args:
        exons_fasta (str): path to fasta containing exon sequences
        introns_fasta (str): path to fasta containing intron sequences
        output_file (str): path to output file
        families_file (str): if set, group results into paralagous_families
    """

    # get the exons and introns for each transcript
    exon_list, intron_list = get_transcript_exons_and_introns(exons_fasta, introns_fasta)
    # calculate the gc contents of the sets
    exons_gc = {name: seqo.calc_gc_seqs_combined(exon_list[name]) for name in exon_list}
    introns_gc = {name: seqo.calc_gc_seqs_combined(intron_list[name]) for name in intron_list}

    codon_sets = [stops]
    args = [len(codon_sets), exons_gc, introns_gc, exon_list, intron_list, output_file, families_file]
    simoc.run_simulation_function(codon_sets, args, calculate_densities, parallel = False, sim_run = False)


def calc_nds(transcript_ids, transcript_list, gc_controls_zip):

    outputs = []

    with zipfile.ZipFile(gc_controls_zip) as z:
        filelist = z.namelist()
        filelist = {i.split("/")[-1][:-4]: i for i in filelist if i.split("/")[-1][:-4] in transcript_ids}

        for id in filelist:
            name = id.split(".")[0]
            exon_id = int(id.split(".")[1])
            focal_seq = transcript_list[id]
            # get the reading frame of the sequence
            reading_frame = seqo.get_exon_reading_frame(focal_seq)
            # get the density of the exon
            density = seqo.calc_seq_stop_density(focal_seq, exclude_frame = reading_frame)
            # now get the density of the gc matched sequences
            sims_ds = []
            with z.open(filelist[id]) as sim_file:
                sims = str(sim_file.readlines()[0].decode("utf-8")).split(",")
                [sims_ds.append(seqo.calc_seq_stop_density(i, exclude_frame = reading_frame)) for i in sims]
            nd = np.divide(density - np.mean(sims_ds), np.mean(sims_ds))
            outputs.append([name, exon_id, nd])

    return outputs

def cds_density_nd(exons_fasta, families_file, gc_controls_zip, output_directory):

    exon_names, exon_seqs = gen.read_fasta(exons_fasta)
    transcript_list = {name.split("(")[0]: exon_seqs[i] for i, name in enumerate(exon_names[:5000])}

    transcript_ids = [i for i in transcript_list]
    args = [transcript_list, gc_controls_zip]
    outputs = simoc.run_simulation_function(transcript_ids, args, calc_nds, sim_run = False)

    nd_list = collections.defaultdict(lambda: [])
    for output in outputs:
        nd_list[output[0]].append(output[2])

    nd_scores = {i: np.mean(nd_list[i]) for i in nd_list}

    families = gen.read_many_fields(families_file, "\t")

    grouped_nd_scores = sequo.group_family_results(nd_scores, families)

    output_file = "{0}/nd_gc_controls.csv".format(output_directory)
    with open(output_file, "w") as outfile:
        outfile.write("id,nd\n")
        [outfile.write("{0},{1}\n".format(i, np.median(grouped_nd_scores[i]))) for i in grouped_nd_scores]


def calc_stop_densities(id_list, exon_list, cds_list, filelist):

    outputs = []



    for i, id in enumerate(id_list):
        # exons = [i for i in exon_list if id in i]
        exons = []
        for i in exon_list:
            if id in i:
                exons.append(i)

        if len(exons):
            print(id)
            cds_seq = cds_list[id]

            # get the start index and length of each exon in the cds
            exon_info = [[cds_seq.index(exon_list[i]), len(exon_list[i])] for i in exons]
            # read in the simulated cds seqs
            if id in filelist:
                sim_cds_seqs = gen.read_many_fields(filelist[id], ",")
                if sim_cds_seqs[0] and len(sim_cds_seqs[0]) > 0:
                    sim_cds_seqs = sim_cds_seqs[0]
                    # create an empty list to hold the simulated exons
                    sim_list = []
                    # for each of the simulated sequences, get the simulated exon sequences

                    for sim_cds in sim_cds_seqs:
                        sim_exon_sequences = [sim_cds[i[0]:i[0] + i[1]] for i in exon_info]
                        sim_list.append(sim_exon_sequences)

                    real_exons = [exon_list[i] for i in exons]
                    gc = seqo.calc_gc_seqs_combined(real_exons)

                    real_density = seqo.calc_seqs_stop_density(real_exons)
                    sim_densities = [seqo.calc_seqs_stop_density(i) for i in sim_list]

                    nd = np.divide(real_density - np.mean(sim_densities), np.mean(sim_densities))
                    outputs.append([id, gc, nd])

    return outputs

def stop_density_nd(exons_fasta, cds_fasta, introns_fasta, dint_control_cds_output_directory, output_file, families_file = None):

    exon_names, exon_seqs = gen.read_fasta(exons_fasta)
    cds_names, cds_seqs = gen.read_fasta(cds_fasta)
    intron_names, intron_seqs = gen.read_fasta(introns_fasta)

    filelist = {i.split(".")[0]: "{0}/{1}".format(dint_control_cds_output_directory, i) for i in os.listdir(dint_control_cds_output_directory) if len(i.split(".")[0])}

    intron_list = {name.split("(")[0]: intron_seqs[i] for i, name in enumerate(intron_names) if name.split(".")[0] in filelist}
    intron_ids = [i.split(".")[0] for i in intron_list]
    exon_list = {name.split("(")[0]: exon_seqs[i] for i, name in enumerate(exon_names) if name.split(".")[0] in filelist and name.split(".")[0] in intron_ids}
    cds_list = {name: cds_seqs[i] for i, name in enumerate(cds_names) if name.split(".")[0] in filelist and name.split(".")[0] in intron_ids}

    # calculate the nd scores for each data point
    id_list = [i for i in filelist]
    args = [exon_list, cds_list, filelist]
    outputs = simoc.run_simulation_function(id_list, args, calc_stop_densities, sim_run = False)

    gc_list = {}
    nd_list = {}
    for output in outputs:
        gc_list[output[0]] = output[1]
        nd_list[output[0]] = output[2]

    grouped_exons = collections.defaultdict(lambda: [])
    [grouped_exons[i.split(".")[0]].append(exon_list[i]) for i in exon_list]
    grouped_introns = collections.defaultdict(lambda: [])
    [grouped_introns[i.split(".")[0]].append(intron_list[i]) for i in intron_list]


    intron_numbers = {}
    intron_sizes = {}
    intron_density = {}
    intron_ratio = {}
    for id in nd_list:
        intron_number = len(grouped_introns[id])
        exon_bp = sum([len(i) for i in grouped_exons[id]])
        intron_bp = sum([len(i) for i in grouped_introns[id]])
        intron_numbers[id] = intron_number
        intron_sizes[id] = np.mean([len(i) for i in grouped_introns[id]])
        intron_density[id] = np.divide(intron_number, exon_bp)
        intron_ratio[id] = np.divide(exon_bp, intron_bp)


    # group by family if exists
    if families_file:
        families = gen.read_many_fields(families_file, "\t")
        gc_list = sequo.group_family_results(gc_list, families)
        nd_list = sequo.group_family_results(nd_list, families)
        intron_numbers = sequo.group_family_results(intron_numbers, families)
        intron_sizes = sequo.group_family_results(intron_sizes, families)
        intron_density = sequo.group_family_results(intron_density, families)
        intron_ratio = sequo.group_family_results(intron_ratio, families)

    # write to file
    with open(output_file, "w") as outfile:
        outfile.write("id,gc,nd,intron_count,mean_intron_size,intron_density,intron_ratio\n")
        for id in nd_list:
            args = [id, np.median(gc_list[id]), np.median(nd_list[id]), np.median(intron_numbers[id]), np.median(intron_sizes[id]), np.median(intron_density[id]), np.median(intron_ratio[id])]
            outfile.write("{0}\n".format(",".join(gen.stringify(args))))



def compare_density_no_ese(exons_fasta, cds_fasta, ese_file, families_file = None):

    # get eses
    eses = [i[0] for i in gen.read_many_fields(ese_file, "\t") if "#" not in i[0]]
    # get cds sequences
    cds_seqs = gen.fasta_to_list(cds_fasta)

    cds_ids = list(cds_seqs.keys())[:30]
    cds_seqs = {i: cds_seqs[i] for i in cds_ids}

    exon_list = gen.fasta_to_list(exons_fasta, split = "(")
    transcript_exon_list = collections.defaultdict(lambda: [])
    [transcript_exon_list[i.split(".")[0]].append(exon_list[i]) for i in exon_list if i.split(".")[0] in cds_seqs]

    # remove all nucleotides that overlap with an intron_densities_scaled
    ese_removed_seqs = sequo.replace_motifs_in_seqs(cds_seqs, motif_set = eses)

    for id in transcript_exon_list:
        cds_seq = cds_seqs[id]
        exon_coordinates = [[cds_seq.index(exon), cds_seq.index(exon) + len(exon)] for exon in transcript_exon_list[id]]
        removed_ese_seq = ese_removed_seqs[id]
        removed_ese_exons = [removed_ese_seq[i[0]:i[1]] for i in exon_coordinates]
        densities = seqo.calc_seqs_codon_set_density(transcript_exon_list[id], codon_set = stops)
        repl_densities = seqo.calc_seq_replaced_codon_set_density(removed_ese_exons, codon_set = stops)
        print(densities, repl_densities)




def calc_motif_sets_ds(motif_sets, codon_sets, sequence_alignments):

    outputs = []

    for i, motif_set in enumerate(motif_sets):
        print("(W{0}) {1}/{2}".format(mp.current_process().name.split("-")[-1], i+1, len(motif_sets)))

        stops_ds, no_stops_ds = sequo.calc_ese_stop_ds(motif_set, codon_sets, sequence_alignments)
        temp_file = "temp_files/{0}.txt".format(random.random())
        outputs.append(temp_file)
        with open(temp_file, "w") as temp:
            temp.write("{0}\t{1}\n".format(stops_ds, no_stops_ds))

    return outputs


def calc_ds(alignment_file, cds_fasta, ortholog_cds_fasta, ortholog_transcript_links, motif_file, output_directory, output_file, motif_controls_directory = None, families_file = None, run_number = None, codon_sets_file = None):

    start_time = time.time()

    # if the sequence alignment file doesnt exist, create it
    if not os.path.isfile(alignment_file):
        sequo.extract_alignments_from_file(cds_fasta, ortholog_cds_fasta, ortholog_transcript_links, alignment_file)

    # set up the output directory
    gen.create_output_directories(output_directory)
    # get all the sequence alignments
    sequence_alignment_names, sequence_alignment_seqs = gen.read_fasta(alignment_file)
    sequence_alignments = {name: sequence_alignment_seqs[i].split(",") for i, name in enumerate(sequence_alignment_names)}

    # if families file, pick a random member of the family
    if families_file:
        if not run_number:
            family_output_choices_file = "{0}/family_choices.txt".format(output_directory)
        else:
            family_output_choices_file = "{0}/family_choices_{1}.txt".format(output_directory, run_number)
        sequence_alignments = sequo.pick_random_family_member(families_file, sequence_alignments, output_file = family_output_choices_file)


    codon_sets = [stops]
    if codon_sets_file:
        gc_matched_sets = gen.read_many_fields(codon_sets_file, "\t")
        purine_matched = sequo.get_purine_matched_sets(stops, gc_matched_sets)
        extra_sets = [i for i in purine_matched if len(list(set(i) & set(stops))) == 0]
        codon_sets = codon_sets + extra_sets

    # create a file for the real outputs
    temp_dir = "temp_ds"
    gen.create_output_directories(temp_dir)
    #
    # # create a list containing sets to test
    real_name = "real"
    motif_sets = {real_name: motif_file}
    if motif_controls_directory:
        for i, file in enumerate(os.listdir(motif_controls_directory)):
            motif_sets[i] = "{0}/{1}".format(motif_controls_directory, file)
    motif_set_list = [i for i in motif_sets]

    # set up the control runs
    kwargs_dict = {"output_directory": temp_dir}
    args = [motif_sets, codon_sets, sequence_alignments]
    # run on the controls
    outputs = simoc.run_simulation_function(motif_set_list, args, sequo.calc_motif_sets_codons_ds_wrapper, kwargs_dict = kwargs_dict, sim_run = False, parallel = True)


    # now write all the results to the output file
    sorted_codon_sets = ["_".join(j) for j in sorted([sorted(i) for i in codon_sets])]
    with open(output_file, "w") as outfile:
        sorted_names_hits_ds = [",".join(["{0}_hits_ds".format(i), "{0}_no_hits_ds".format(i), "{0}_hits_query_count".format(i), "{0}_no_hits_query_count".format(i)]) for i in sorted_codon_sets]
        outfile.write("sim_id,{0}\n".format(",".join(sorted_names_hits_ds)))
        for file in outputs:
            results = gen.read_many_fields(file, ",")[1:]
            results = {i[0]: i[1:] for i in results}
            # get the id of the results
            sim_no = file.split("/")[-1].split(".")[0]
            if sim_no == real_name:
                id = "real"
            else:
                id = "sim_{0}".format(int(sim_no)+1)
            outfile.write("{0}".format(id))
            for codon_set in sorted_codon_sets:
                outfile.write(",{0}".format(",".join(gen.stringify(results[codon_set]))))
            outfile.write("\n")

    # # remove all the temp files
    [gen.remove_file(i) for i in outputs]

    gen.get_time(start_time)


def calc_codon_ds(alignment_file, cds_fasta, ortholog_cds_fasta, ortholog_transcript_links, motif_file, output_directory, output_file, motif_controls_directory = None, families_file = None, run_number = None, codon_sets_file = None):

    start_time = time.time()

    # if the sequence alignment file doesnt exist, create it
    if not os.path.isfile(alignment_file):
        sequo.extract_alignments_from_file(cds_fasta, ortholog_cds_fasta, ortholog_transcript_links, alignment_file)

    # set up the output directory
    gen.create_output_directories(output_directory)
    # get all the sequence alignments
    sequence_alignment_names, sequence_alignment_seqs = gen.read_fasta(alignment_file)
    sequence_alignments = {name: sequence_alignment_seqs[i].split(",") for i, name in enumerate(sequence_alignment_names)}

    # if families file, pick a random member of the family
    if families_file:
        if not run_number:
            family_output_choices_file = "{0}/family_choices_codons.txt".format(output_directory)
        else:
            family_output_choices_file = "{0}/family_choices_codons_{1}.txt".format(output_directory, run_number)
        sequence_alignments = sequo.pick_random_family_member(families_file, sequence_alignments, output_file = family_output_choices_file)


    codon_sets = [stops]
    if codon_sets_file:
        gc_matched_sets = gen.read_many_fields(codon_sets_file, "\t")
        purine_matched = sequo.get_purine_matched_sets(stops, gc_matched_sets)
        # extra_sets = [i for i in purine_matched if len(list(set(i) & set(stops))) == 0]
        extra_sets = purine_matched
        codon_sets = codon_sets + extra_sets

    # create a file for the real outputs
    temp_dir = "temp_ds"
    gen.create_output_directories(temp_dir)

    motif_set = [i[0] for i in gen.read_many_fields(motif_file, "\t") if "#" not in i[0] and ">" not in i[0]]
    # now get all the parts of the sequences where only the motif sets occur
    restricted_sequences = sequo.extract_motif_sequences_from_alignment(sequence_alignments, motif_set)

    # set up the control runs
    kwargs_dict = {"output_directory": temp_dir}
    args = [restricted_sequences]
    # run on the controls
    outputs = simoc.run_simulation_function(codon_sets, args, sequo.calc_codon_sets_ds_wrapper, kwargs_dict = kwargs_dict, sim_run = False, parallel = True)


    # now write all the results to the output file
    with open(output_file, "w") as outfile:
        outfile.write("codon_set,hits_ds,no_hits_ds,hits_query_count,no_hits_query_count\n")
        for file in sorted(outputs):
            results = gen.read_many_fields(file, ",")[1:]
            results = {i[0]: i[1:] for i in results}
            [outfile.write("{0},{1}\n".format(i, ",".join(gen.stringify(results[i])))) for i in sorted(results)]

    # # remove all the temp files
    [gen.remove_file(i) for i in outputs]

    gen.get_time(start_time)


def exon_region_excess(cds_fasta, exons_fasta, output_file, families_file = None):

    cds_seqs = gen.fasta_to_list(cds_fasta)
    exon_names, exon_seqs = gen.read_fasta(exons_fasta)
    exons_list = collections.defaultdict(lambda: [])
    [exons_list[name.split(".")[0]].append(exon_seqs[i]) for i, name in enumerate(exon_names)]

    exon_seqs = {i: [seq for seq in exons_list[i] if len(seq) > 211] for i in exons_list}

    # if families file, pick a random member of the family
    if families_file:
        exon_seqs = sequo.pick_random_family_member(families_file, exon_seqs)


    # set up the dictionaries to hold info
    upstream_flanks = collections.defaultdict(lambda: [])
    downstream_flanks = collections.defaultdict(lambda: [])
    cores = collections.defaultdict(lambda: [])


    for o, id in enumerate(exon_seqs):
        for i, seq in enumerate(exon_seqs[id]):
            # upstream_flanks[id].append(seq)
            cds_index = cds_seqs[id].index(seq)
            exon_start_frame = cds_index % 3

            if exon_start_frame == 0:
                upstream_start_index = 0
            elif exon_start_frame == 1:
                upstream_start_index = 2
            elif exon_start_frame == 2:
                upstream_start_index = 1
            upstream_end_index = upstream_start_index + 69

            exon_end_frame = (cds_index + len(seq)) % 3
            if exon_end_frame == 2:
                # +1 here becuase you want the full end codon
                downstream_end_index = len(seq) + 1
            elif exon_end_frame == 1:
                downstream_end_index = len(seq) - 1
            elif exon_end_frame == 0:
                downstream_end_index = len(seq)
            downstream_start_index = downstream_end_index - 69

            core_region = seq[upstream_end_index:downstream_start_index]

            remainder = len(core_region)-69
            if remainder % 2 == 0:
                core_start_index = upstream_end_index + int(np.divide(remainder,2))
            else:
                # to keep in frame, take distance of n nucleotides from end of upstream exon
                # and n + 3 nucleotides from start of downstream exon
                # this is the equivalent of (remainder - 3)/2
                core_start_index = upstream_end_index + int(np.divide(remainder-3, 2))

            core_end_index = core_start_index + 69


            upstream_flanks[id].append(seq[upstream_start_index:upstream_end_index])
            downstream_flanks[id].append(seq[downstream_start_index:downstream_end_index])
            cores[id].append(seq[core_start_index:core_end_index])
    #
    # for id in upstream_flanks:
    #     for i in upstream_flanks[id]:
    #         print(len(i))

    core_densities = {id: seqo.calc_seqs_stop_density(cores[id]) for id in cores}
    upstream_counts = {id: seqo.calc_seqs_stop_count(upstream_flanks[id]) for id in upstream_flanks}
    downstream_counts = {id: seqo.calc_seqs_stop_count(downstream_flanks[id]) for id in downstream_flanks}

    expected_upstream_counts = {}
    expected_downstream_counts = {}
    for id in core_densities:
        expected_upstream_counts[id] = np.divide(core_densities[id] * sum([len(i) for i in upstream_flanks[id]]), 3)
        expected_downstream_counts[id] = np.divide(core_densities[id] * sum([len(i) for i in downstream_flanks[id]]), 3)


    with open(output_file, "w") as outfile:
        outfile.write("id,core_density,expected_upstream_counts,upstream_counts,expected_downstream_counts,downstream_counts\n")
        for id in core_densities:
            outfile.write("{0},{1},{2},{3},{4},{5}\n".format(id, core_densities[id], expected_upstream_counts[id], upstream_counts[id], expected_downstream_counts[id], downstream_counts[id]))
    #
    # expected_upstream_densities = {id: [len(i) for i in upstream_flanks[id]] for id in upstream_flanks}
    # print(expected_upstream_densities)




def exon_region_density(cds_fasta, exons_fasta, gc_matched_stops_file, output_file, families_file = None):


    cds_seqs = gen.fasta_to_list(cds_fasta)

    exon_names, exon_seqs = gen.read_fasta(exons_fasta)
    exons_list = collections.defaultdict(lambda: [])
    [exons_list[name.split(".")[0]].append(exon_seqs[i]) for i, name in enumerate(exon_names)]

    exon_seqs = {i: [seq for seq in exons_list[i] if len(seq) > 211] for i in exons_list}


    # if families file, pick a random member of the family
    if families_file:
        exon_seqs = sequo.pick_random_family_member(families_file, exon_seqs)

    # set up the dictionaries to hold info
    upstream_flanks = collections.defaultdict(lambda: [])
    downstream_flanks = collections.defaultdict(lambda: [])
    cores = collections.defaultdict(lambda: [])


    for o, id in enumerate(exon_seqs):
        for i, seq in enumerate(exon_seqs[id]):
            # upstream_flanks[id].append(seq)
            cds_index = cds_seqs[id].index(seq)
            exon_start_frame = cds_index % 3

            if exon_start_frame == 0:
                upstream_start_index = 0
            elif exon_start_frame == 1:
                upstream_start_index = 2
            elif exon_start_frame == 2:
                upstream_start_index = 1
            upstream_end_index = upstream_start_index + 69

            exon_end_frame = (cds_index + len(seq)) % 3
            if exon_end_frame == 2:
                # +1 here becuase you want the full end codon
                downstream_end_index = len(seq) + 1
            elif exon_end_frame == 1:
                downstream_end_index = len(seq) - 1
            elif exon_end_frame == 0:
                downstream_end_index = len(seq)
            downstream_start_index = downstream_end_index - 69

            core_region = seq[upstream_end_index:downstream_start_index]

            remainder = len(core_region)-69
            if remainder % 2 == 0:
                core_start_index = upstream_end_index + int(np.divide(remainder,2))
            else:
                # to keep in frame, take distance of n nucleotides from end of upstream exon
                # and n + 3 nucleotides from start of downstream exon
                # this is the equivalent of (remainder - 3)/2
                core_start_index = upstream_end_index + int(np.divide(remainder-3, 2))

            core_end_index = core_start_index + 69

            # print(upstream_start_index, upstream_end_index)

            upstream_flanks[id].append(seq[upstream_start_index:upstream_end_index])
            downstream_flanks[id].append(seq[downstream_start_index:downstream_end_index])
            cores[id].append(seq[core_start_index:core_end_index])


    # gc_matched_stops = gen.read_many_fields(gc_matched_stops_file, "\t")
    # matched_stops = sequo.get_purine_matched_sets(stops, gc_matched_stops)
    #
    # matched_stops = [i for i in matched_stops if len(list(set(stops) & set(i))) == 0]
    # codon_sets = [stops] + matched_stops
    codon_sets = [stops]

    upstream_flanks = {i: upstream_flanks[i] for i in upstream_flanks}
    cores = {i: cores[i] for i in cores}
    downstream_flanks = {i: downstream_flanks[i] for i in downstream_flanks}

    args = [upstream_flanks, cores, downstream_flanks]
    outputs = simoc.run_simulation_function(codon_sets, args, calc_exon_region_densities, sim_run = False)

    # output_list = {i.split("/")[-1].split(".")[0]: i for i in outputs}
    #
    # density_list = collections.defaultdict(lambda: [])
    # for codon_set in output_list:
    densities = gen.read_many_fields(outputs[0], ",")
    densities = {i[0]: i[1:] for i in densities}
    # print(densities)
    # families = gen.read_many_fields(families_file, "\t")
    # densities = sequo.group_family_results(densities, families)
    # density_outputs = collections.defaultdict(lambda: collections.defaultdict(lambda: []))
    # for id in densities:
    #     cases = densities[id]
    #     for i in cases:
    #         for j in range(len(i)):
    #             density_outputs[id][j].append(float(i[j]))
        # for case in densities[id]:
        #     print(id, case, densities[id][case])
            # for i in range(len(densities[id][case])):
            #     density_outputs[id][i].append(densities[id][case][i])


        # [density_list[i[0]].extend(i[1:]) for i in densities]

    with open(output_file, "w") as outfile:
        outfile.write("id,upstream_density,core_density,downstream_density\n")
        for id in densities:
            outfile.write("{0},{1}\n".format(id, ",".join(densities[id])))

    [gen.remove_file(i) for i in outputs]

def calc_exon_region_densities(codon_sets, upstream_flanks, cores, downstream_flanks):

    temp_dir = "temp_files"
    gen.create_output_directories(temp_dir)

    outputs = []

    for i, codon_set in enumerate(codon_sets):
        print("(W{0}) {1}/{2}".format(mp.current_process().name.split("-")[-1], i+1, len(codon_sets)))

        temp_file = "{0}/{1}.txt".format(temp_dir, "_".join(sorted(codon_set)), random.random())
        outputs.append(temp_file)
        if not os.path.isfile(temp_file):
            upstream_densities = {id: seqo.calc_seqs_codon_set_density(upstream_flanks[id], codon_set = codon_set, exclude_frames = 0) for id in upstream_flanks}
            downstream_densities = {id: seqo.calc_seqs_codon_set_density(downstream_flanks[id], codon_set = codon_set, exclude_frames = 0) for id in downstream_flanks}
            core_densities = {id: seqo.calc_seqs_codon_set_density(cores[id], codon_set = codon_set, exclude_frames = 0) for id in cores}
            with open(temp_file, "w") as outfile:
                for id in sorted(upstream_densities):
                    outfile.write("{0},{1},{2},{3}\n".format(id, np.median(upstream_densities[id]), np.median(core_densities[id]), np.median(downstream_densities[id])))

    return outputs


def calc_ese_ds(alignment_file, cds_fasta, exons_fasta, motif_file, output_file, simulations = None, controls_directory = None, families_file = None):
    """
    Calculate the ds scores in ESEs.

    Args:
        alignment_file (str): path to file containing sequence alignments
        exons_fasta (str): path to file contaning exon sequences
        motif_file (str): path to file containing motifs
        output_file (str): path to output file
        simulations (int): if set, the number of simulants to run
        controls_directory (path): if set, the path to simulants
        families_file (path): if set, use to choose one gene member per family
    """

    start_time = time.time()

    # get all the names of the multi exon sequences
    exon_names = gen.read_fasta(exons_fasta)[0]
    # now get all the sequence alignments which match the names of the multi exon genes
    alignment_names, alignment_seqs = gen.read_fasta(alignment_file)
    alignments = {name: alignment_seqs[i].split(",") for i, name in enumerate(alignment_names) if name in exon_names}

    # if we want to group into paralagous families
    if families_file:
        alignments = sequo.pick_random_family_member(families_file, alignments)

    # now create a list filepaths to motif sets test
    # we first want the real motifs, then the number of simulated sets
    motif_set_filepaths = {"real": motif_file}
    if simulations and controls_directory:
        for i, file in enumerate(os.listdir(controls_directory)[:simulations]):
            motif_set_filepaths[i+1] = "{0}/{1}".format(controls_directory, file)
    elif simulations and not controls_directory:
        print("\nPlease provide both the number of simulations and the controls directory...\n")
        raise Exception

    # set up the temporary output_directory to hold results
    temp_output_directory = "temp_ds_outputs"
    gen.create_output_directories(temp_output_directory)

    # set up the test
    args = [motif_set_filepaths, alignments, temp_output_directory]
    # run the function
    outputs = simoc.run_simulation_function(list(motif_set_filepaths), args, sequo.ese_ds_wrapper, sim_run = False)

    with open(output_file, "w") as outfile:
        outfile.write("id,ese_ds,non_ese_ds,ese_stop_ds,ese_non_stop_ds\n")
        for file in outputs:
            outfile.write("{0}\n".format(",".join(gen.read_many_fields(file, "\t")[0])))

    # remove temp output directory
    gen.remove_directory(temp_output_directory)
    # get the finishing time
    gen.get_time(start_time)

def calc_ese_ds_random_overlaps(alignment_file, cds_fasta, exons_fasta, motif_file, output_file, simulations = None, controls_directory = None, families_file = None):
    """
    Calculate the ds scores in ESEs.

    Args:
        alignment_file (str): path to file containing sequence alignments
        exons_fasta (str): path to file contaning exon sequences
        motif_file (str): path to file containing motifs
        output_file (str): path to output file
        simulations (int): if set, the number of simulants to run
        controls_directory (path): if set, the path to simulants
        families_file (path): if set, use to choose one gene member per family
    """

    start_time = time.time()

    # get all the names of the multi exon sequences
    exon_names = gen.read_fasta(exons_fasta)[0]
    # now get all the sequence alignments which match the names of the multi exon genes
    alignment_names, alignment_seqs = gen.read_fasta(alignment_file)
    alignments = {name: alignment_seqs[i].split(",") for i, name in enumerate(alignment_names) if name in exon_names}

    # if we want to group into paralagous families
    if families_file:
        alignments = sequo.pick_random_family_member(families_file, alignments)

    # now create a list filepaths to motif sets test
    # we first want the real motifs, then the number of simulated sets
    motif_set_filepaths = {"real": motif_file}
    if simulations and controls_directory:
        for i, file in enumerate(os.listdir(controls_directory)[:simulations]):
            motif_set_filepaths[i+1] = "{0}/{1}".format(controls_directory, file)
    elif simulations and not controls_directory:
        print("\nPlease provide both the number of simulations and the controls directory...\n")
        raise Exception

    # set up the temporary output_directory to hold results
    temp_output_directory = "temp_ds_random_overlaps_outputs"
    gen.create_output_directories(temp_output_directory)

    # set up the test
    args = [motif_set_filepaths, alignments, temp_output_directory]
    # run the function
    outputs = simoc.run_simulation_function(list(motif_set_filepaths), args, sequo.ese_ds_random_overlaps_wrapper, sim_run = False)

    with open(output_file, "w") as outfile:
        outfile.write("id,ese_ds,non_ese_ds,ese_stop_ds,ese_non_stop_ds\n")
        for file in outputs:
            outfile.write("{0}\n".format(",".join(gen.read_many_fields(file, "\t")[0])))

    # remove temp output directory
    gen.remove_directory(temp_output_directory)
    # get the finishing time
    gen.get_time(start_time)

def calc_non_ese_ds(alignment_file, cds_fasta, exons_fasta, motif_file, output_file, simulations = None, controls_directory = None, families_file = None):
    """
    Calculate the ds scores for stop codons not in ESEs.

    Args:
        alignment_file (str): path to file containing sequence alignments
        exons_fasta (str): path to file contaning exon sequences
        motif_file (str): path to file containing motifs
        output_file (str): path to output file
        simulations (int): if set, the number of simulants to run
        controls_directory (path): if set, the path to simulants
        families_file (path): if set, use to choose one gene member per family
    """

    start_time = time.time()

    # get all the names of the multi exon sequences
    exon_names = gen.read_fasta(exons_fasta)[0]
    # now get all the sequence alignments which match the names of the multi exon genes
    alignment_names, alignment_seqs = gen.read_fasta(alignment_file)
    alignments = {name: alignment_seqs[i].split(",") for i, name in enumerate(alignment_names) if name in exon_names}

    # if we want to group into paralagous families
    if families_file:
        alignments = sequo.pick_random_family_member(families_file, alignments)

    # now create a list filepaths to motif sets test
    # we first want the real motifs, then the number of simulated sets
    motif_set_filepaths = {"real": motif_file}
    if simulations and controls_directory:
        for i, file in enumerate(os.listdir(controls_directory)[:simulations]):
            motif_set_filepaths[i+1] = "{0}/{1}".format(controls_directory, file)
    elif simulations and not controls_directory:
        print("\nPlease provide both the number of simulations and the controls directory...\n")
        raise Exception

    # set up the temporary output_directory to hold results
    temp_output_directory = "temp_ds_outputs"
    gen.create_output_directories(temp_output_directory)

    # set up the test
    args = [motif_set_filepaths, alignments, temp_output_directory]
    # run the function
    outputs = simoc.run_simulation_function(list(motif_set_filepaths), args, sequo.non_ese_ds_wrapper, sim_run = False)

    with open(output_file, "w") as outfile:
        outfile.write("id,non_ese_stop_ds,non_ese_non_stop_ds\n")
        for file in outputs:
            outfile.write("{0}\n".format(",".join(gen.read_many_fields(file, "\t")[0])))

    # remove temp output directory
    gen.remove_directory(temp_output_directory)
    # get the finishing time
    gen.get_time(start_time)

def calc_ese_ds_mutation(alignment_file, cds_fasta, exons_fasta, motif_file, output_file, simulations = None, controls_directory = None, families_file = None):
    """
    Calculate the ds scores in ESEs.

    Args:
        alignment_file (str): path to file containing sequence alignments
        exons_fasta (str): path to file contaning exon sequences
        motif_file (str): path to file containing motifs
        output_file (str): path to output file
        simulations (int): if set, the number of simulants to run
        controls_directory (path): if set, the path to simulants
        families_file (path): if set, use to choose one gene member per family
    """

    start_time = time.time()

    # get all the names of the multi exon sequences
    exon_names = gen.read_fasta(exons_fasta)[0]
    # now get all the sequence alignments which match the names of the multi exon genes
    alignment_names, alignment_seqs = gen.read_fasta(alignment_file)
    alignments = {name: alignment_seqs[i].split(",") for i, name in enumerate(alignment_names) if name in exon_names}

    # if we want to group into paralagous families
    if families_file:
        alignments = sequo.pick_random_family_member(families_file, alignments)

    # now create a list filepaths to motif sets test
    # we first want the real motifs, then the number of simulated sets
    motif_set_filepaths = {"real": motif_file}
    if simulations and controls_directory:
        for i, file in enumerate(os.listdir(controls_directory)[:simulations]):
            motif_set_filepaths[i+1] = "{0}/{1}".format(controls_directory, file)
    elif simulations and not controls_directory:
        print("\nPlease provide both the number of simulations and the controls directory...\n")
        raise Exception

    # set up the temporary output_directory to hold results
    temp_output_directory = "temp_d__mutation_outputs"
    gen.create_output_directories(temp_output_directory)

    # set up the test
    args = [motif_set_filepaths, alignments, temp_output_directory]
    # run the function
    outputs = simoc.run_simulation_function(list(motif_set_filepaths), args, sequo.ese_ds_mutation_wrapper, sim_run = False)

    with open(output_file, "w") as outfile:
        outfile.write("id,one_away_ds,others_ds\n")
        for file in outputs:
            outfile.write("{0}\n".format(",".join(gen.read_many_fields(file, "\t")[0])))

    # remove temp output directory
    gen.remove_directory(temp_output_directory)
    # get the finishing time
    gen.get_time(start_time)





def calc_ds_all(simulations, alignment_file, cds_fasta, mutli_exon_fasta, ortholog_cds_fasta, ortholog_transcript_links, motif_file, output_directory, output_file, motif_controls_directory = None, families_file = None, run_number = None, codon_sets_file = None):

    start_time = time.time()

    # if the sequence alignment file doesnt exist, create it
    if not os.path.isfile(alignment_file):
        sequo.extract_alignments_from_file(cds_fasta, ortholog_cds_fasta, ortholog_transcript_links, alignment_file)

    # set up the output directory
    gen.create_output_directories(output_directory)
    # get all the sequence alignments
    # first get all the names of the cases that are multi exon genes
    multi_exons_cds_names = gen.read_fasta(mutli_exon_fasta)[0]
    # then get the alignments that match these genes
    sequence_alignment_names, sequence_alignment_seqs = gen.read_fasta(alignment_file)
    sequence_alignments = {name: sequence_alignment_seqs[i].split(",") for i, name in enumerate(sequence_alignment_names) if name in multi_exons_cds_names}

    # if families file, pick a random member of the family
    if families_file:
        family_output_choices_file = "{0}/family_choices.txt".format(output_directory)
        sequence_alignments = sequo.pick_random_family_member(families_file, sequence_alignments, output_file = family_output_choices_file)

    codon_sets = [stops]

    # create a file for the real outputs
    temp_dir = "temp_ds"
    gen.create_output_directories(temp_dir)

    # create a list containing sets to test
    real_name = "real"
    motif_sets = {real_name: motif_file}
    if motif_controls_directory:
        for i, file in enumerate(os.listdir(motif_controls_directory)[:simulations]):
            motif_sets[i] = "{0}/{1}".format(motif_controls_directory, file)
    motif_set_list = [i for i in motif_sets]
    #
    # set up the control runs
    kwargs_dict = {"output_directory": temp_dir}
    args = [motif_sets, codon_sets, sequence_alignments]
    # run on the controls
    outputs = simoc.run_simulation_function(motif_set_list, args, sequo.calc_motif_sets_all_ds_wrapper, kwargs_dict = kwargs_dict, sim_run = False, parallel = True)

    # now write all the results to the output file
    # sorted_codon_sets = ["_".join(j) for j in sorted([sorted(i) for i in codon_sets])]
    with open(output_file, "w") as outfile:
        outfile.write("id,all_ese_ds,non_ese_ds,stops_ese_ds,non_stops_ese_ds\n")
        for i, file in enumerate(outputs):
            results = gen.read_many_fields(file, ",")[0]
            if results[0] != "real":
                results[0] = int(results[0])+1
            outfile.write("{0},{1}\n".format(results[0], ",".join(results[1:])))

    [gen.remove_file(i) for i in outputs]

    gen.get_time(start_time)

def non_coding_exons(fasta_file, output_file, families_file = None):

    exons = gen.fasta_to_list(fasta_file, split = "(")
    exon_list = collections.defaultdict(lambda: [])
    for i, id in enumerate(exons):
        exon_list[id.split(".")[0]].append(exons[id])

    density = {id: seqo.calc_seqs_stop_density(exon_list[id]) for id in exon_list}
    scaled_density = {id: seqo.calc_intron_seqs_stop_density(exon_list[id], codon_set = stops) for id in exon_list}

    if families_file:
        families = gen.read_many_fields(families_file, "\t")
        density = sequo.group_family_results(density, families)
        scaled_density = sequo.group_family_results(scaled_density, families)

    with open(output_file, "w") as outfile:
        outfile.write("id,density,scaled_density\n")
        [outfile.write("{0},{1},{2}\n".format(id, np.median(density[id]), np.median(scaled_density[id]))) for id in density]

def calculate_motif_densities(iteration_list, filelist, sequence_list, clean_run, output_directory = None):

    if not output_directory:
        output_directory = "temp_dir"
    gen.create_output_directories(output_directory)

    output_files = []


    for i, iteration in enumerate(iteration_list):
        gen.print_parallel_status(i, iteration_list)
        output_file = "{0}/sim{1}.txt".format(output_directory, iteration)
        output_files.append(output_file)
        if not os.path.isfile(output_file) or clean_run:
            motif_set = [i[0] for i in gen.read_many_fields(filelist[iteration], "\t") if "#" not in i[0] and ">" not in i[0]]
            density = {id: seqo.calc_motif_density(sequence_list[id], motif_set) for id in sequence_list}
            pickle.dump(density, open(output_file, "wb" ) )
        # densities[iteration] = density
    return output_files

def calculate_intron_densities(motif_file, introns_fasta, output_file, controls_dir, families_file = None, required_simulations = None, combined = None, clean_run = None):

    if not required_simulations:
        required_simulations = 1000

    # generate random motifs if they dont already exist
    if len(os.listdir(controls_dir)) < required_simulations:
        print("Generating control motif sets...")
        required_sets = required_simulations - len(os.listdir(controls_dir))
        sequo.generate_random_motifs(motif_file, controls_dir, required_sets, stop_restricted = True, exclude_motif_set = True)

    # get the list of introns
    intron_names, intron_seqs = gen.read_fasta(introns_fasta)
    introns = collections.defaultdict(lambda: [])
    [introns[name.split(".")[0]].append(intron_seqs[i]) for i, name in enumerate(intron_names)]

    introns = sequo.pick_random_family_member(families_file, introns)

    combined = False
    if combined:
        intron_list = {"all": []}
        [intron_list["all"].extend(introns[i]) for i in introns]
        introns = intron_list
    else:
        introns = {i: introns[i] for i in introns}

    # get a list of random motif files, and a list of ids
    motif_sets = {i+1: "{0}/{1}".format(controls_dir, file) for i, file in enumerate(os.listdir(controls_dir)[:required_simulations])}
    motif_set_list = [i for i in motif_sets]

    real_motifs = sequo.read_motifs(motif_file)
    real_densities = {i: seqo.calc_motif_density(introns[i], real_motifs) for i in introns}

    sim_output_directory = "temp_motif_dir"
    args = [motif_sets, introns, clean_run]
    kwargs = {"output_directory": sim_output_directory}
    outputs = simoc.run_simulation_function(motif_set_list, args, calculate_motif_densities, kwargs_dict = kwargs, sim_run = False, parallel = True)

    results = collections.defaultdict(lambda: [])
    for id in real_densities:
        results[id].append(real_densities[id])

    # print(len(real_densities))

    for output in outputs:
        sim_id = int(output.split("/")[-1].split(".")[0][3:])
        densities = pickle.load( open(output, "rb" ) )

        for id in results:
            results[id].append(densities[id])

    with open(output_file, "w") as outfile:
        outfile.write("id,real,{0}\n".format(",".join(["sim_{0}".format(i+1) for i in range(len(outputs))])))
        for id in sorted(results):
            outfile.write("{0},{1},{2}\n".format(id, results[id][0], ",".join(gen.stringify(results[id][1:]))))



def calc_purine_content(exons_fasta, introns_fasta, output_file, families_file = None):

    intron_names, intron_seqs = gen.read_fasta(introns_fasta)
    introns = collections.defaultdict(lambda: [])
    [introns[name.split('.')[0]].append(intron_seqs[i]) for i, name in enumerate(intron_names)]

    exon_names, exon_seqs = gen.read_fasta(exons_fasta)
    exons = collections.defaultdict(lambda: [])
    [exons[name.split('.')[0]].append(exon_seqs[i]) for i, name in enumerate(exon_names) if name.split(".")[0] in introns]


    exon_purine_content = {i: sequo.calc_purine_content(exons[i]) for i in exons}
    intron_purine_content = {i: sequo.calc_purine_content(introns[i]) for i in introns}

    if families_file:
        families = gen.read_many_fields(families_file, "\t")
        exon_purine_content = sequo.group_family_results(exon_purine_content, families)
        intron_purine_content = sequo.group_family_results(intron_purine_content, families)

    with open(output_file, "w") as outfile:
        outfile.write("id,exon_purine_content,intron_purine_content\n")
        [outfile.write("{0},{1},{2}\n".format(i, np.median(exon_purine_content[i]), np.median(intron_purine_content[i]))) for i in sorted(intron_purine_content) if i in exon_purine_content]



def calc_ds_mutation(simulations, alignment_file, cds_fasta, ortholog_cds_fasta, ortholog_transcript_links, motif_file, output_directory, output_file, motif_controls_directory = None, families_file = None, run_number = None, codon_sets_file = None):

    start_time = time.time()

    # if the sequence alignment file doesnt exist, create it
    if not os.path.isfile(alignment_file):
        sequo.extract_alignments_from_file(cds_fasta, ortholog_cds_fasta, ortholog_transcript_links, alignment_file)

    # set up the output directory
    gen.create_output_directories(output_directory)
    # get all the sequence alignments
    sequence_alignment_names, sequence_alignment_seqs = gen.read_fasta(alignment_file)
    sequence_alignments = {name: sequence_alignment_seqs[i].split(",") for i, name in enumerate(sequence_alignment_names)}

    # if families file, pick a random member of the family
    if families_file:
        if not run_number:
            family_output_choices_file = "{0}/family_choices.txt".format(output_directory)
        else:
            family_output_choices_file = "{0}/family_choices_stats_{1}.txt".format(output_directory, run_number)
        sequence_alignments = sequo.pick_random_family_member(families_file, sequence_alignments, output_file = family_output_choices_file)

    codon_sets = [stops]

    # create a file for the real outputs
    temp_dir = "temp_ds"
    gen.create_output_directories(temp_dir)

    # create a list containing sets to test
    real_name = "real"
    motif_sets = {real_name: motif_file}
    if motif_controls_directory:
        for i, file in enumerate(os.listdir(motif_controls_directory)[:simulations]):
            motif_sets[i] = "{0}/{1}".format(motif_controls_directory, file)
    motif_set_list = [i for i in motif_sets]
    #
    # set up the control runs
    kwargs_dict = {"output_directory": temp_dir}
    args = [motif_sets, codon_sets, sequence_alignments]
    # run on the controls
    outputs = simoc.run_simulation_function(motif_set_list, args, sequo.calc_ds_mutation_wrapper, kwargs_dict = kwargs_dict, sim_run = False, parallel = True)

    # now write all the results to the output file
    # sorted_codon_sets = ["_".join(j) for j in sorted([sorted(i) for i in codon_sets])]
    with open(output_file, "w") as outfile:
        outfile.write("id,stop_mutation_ds,non_stops_mutation_ds\n")
        for i, file in enumerate(outputs):
            results = gen.read_many_fields(file, ",")[0]
            if results[0] != "real":
                results[0] = int(results[0])+1
            outfile.write("{0},{1}\n".format(results[0], ",".join(results[1:])))


def ese_ds(alignments_fasta, cds_fasta, motif_file, output_file, families_file = None):
    """
    Wrapper for calculating the ds scores in motifs and also ds scores for stop
    codons in motifs.

    Args:
        alignments_fasta (str): path to alignments file
        cds_fasta (str): path to file containing cds sequences
        motif_set (str): path to file containing motifs
        output_file (str): path to output file
        families_file (str): if set, the path
    """
    # first get a list of all multi-exon cds
    cds_names = gen.read_fasta(cds_fasta)[0]
    # now get the alignments that correspond to these sequences
    alignment_ids, alignment_sequences = gen.read_fasta(alignments_fasta)
    alignment_sequences = {id: alignment_sequences[i].split(",") for i, id in enumerate(alignment_ids) if id in cds_names}
    # if the families file exists, pick only one member at random to keep
    if families_file:
        sequo.pick_random_family_member(families_file, alignment_sequences)
    # now generate the output for the ds scores
    # first set up a dictionary containing the paths to the motif files
    motif_sets = {"real": motif_file}
    # and now get the ids of these motif sets
    motif_set_ids = list(motif_sets)
    # set up a temporary directory to hold the outputs
    temp_dir = "temp_ds_outputs_dir"
    gen.create_output_directories(temp_dir)
    # run the test
    args = [alignment_sequences, motif_sets, temp_dir]
    outputs = simoc.run_simulation_function(motif_set_ids, args, sequo.calculate_motifs_stop_ds, sim_run = False)

    with open(output_file, "w") as outfile:
        outfile.write("id,hit_ds,hit_stop_ds,hit_non_stop_ds,non_hit_ds,non_hit_stop_ds,non_hit_non_stop_ds\n")
        for file in outputs:
            data = gen.read_many_fields(file, "\t")[0][0]
            outfile.write("{0}\n".format(data))
    # remove the temp directory
    gen.remove_directory(temp_dir)


def intron_length_test(exons_fasta, introns_fasta, motif_file, output_file, flanks = None, families_file = None, restrict_size = None):

    # get the exon fasta
    exon_names, exon_seqs = gen.read_fasta(exons_fasta)
    # get the introns fasta
    intron_names, intron_seqs = gen.read_fasta(introns_fasta)

    # if limited to flanks, extract the flanks
    if flanks:
        # get all the exon flanks
        exon_list = collections.defaultdict(lambda: collections.defaultdict())
        for i, name in enumerate(exon_names):
            if len(exon_seqs[i]) > 211 and "N" not in exon_seqs[i]:
                exon_list[name.split(".")[0]][int(name.split(".")[-1].split("(")[0])] = exon_seqs[i]
        exon_list = {i: exon_list[i] for i in exon_list}


        intron_list = collections.defaultdict(lambda: [])

        # get only those that flank an included exon
        for i, name in enumerate(intron_names):
            intron_id = name.split(".")[0]
            if intron_id in exon_list:
                intron_flanking_exons = [int(i) for i in name.split(".")[-1].split("(")[0].split("-")]
                if intron_flanking_exons[0] in exon_list[intron_id] or intron_flanking_exons[1] in exon_list[intron_id]:
                    intron_list[name.split(".")[0]].append(intron_seqs[i])
        intron_list = {i: intron_list[i] for i in intron_list}

        temp_exon_list = collections.defaultdict(lambda: [])
        for id in exon_list:
            flanks = []
            for exon_id in exon_list[id]:
                flanks.append(exon_list[id][exon_id][2:69])
                flanks.append(exon_list[id][exon_id][-69:-2])
            temp_exon_list[id] = flanks
        exon_list = temp_exon_list
        exon_list = {i: exon_list[i] for i in exon_list}

    else:
        # if not flanks or is non coding exons
        exon_list = collections.defaultdict(lambda: collections.defaultdict())
        for i, name in enumerate(exon_names):
            if "N" not in exon_seqs[i]:
                if restrict_size and len(exon_seqs[i]) > 211:
                    exon_list[name.split(".")[0]][int(name.split(".")[-1].split("(")[0])] = exon_seqs[i]
            else:
                exon_list[name.split(".")[0]][int(name.split(".")[-1].split("(")[0])] = exon_seqs[i]
        exon_list = {i: exon_list[i] for i in exon_list}

        intron_list = collections.defaultdict(lambda: [])
        [intron_list[name.split(".")[0]].append(intron_seqs[i]) for i, name in enumerate(intron_names) if name.split(".")[0] in exon_list]
        intron_list = {i: intron_list[i] for i in intron_list}

        temp_exon_list = collections.defaultdict(lambda: [])
        for id in exon_list:
            exons = []
            for exon_id in exon_list[id]:
                exons.append(exon_list[id][exon_id])
            temp_exon_list[id] = exons
        exon_list = temp_exon_list
        exon_list = {i: exon_list[i] for i in exon_list}

    run_intron_length_density_test(motif_file, exon_list, intron_list, output_file, families_file = families_file)



def run_intron_length_density_test(motif_file, exon_list, intron_list, output_file, families_file = None):

    # get the sizes
    exon_sizes = {i: np.median([len(exon) for exon in exon_list[i]]) for i in exon_list}
    intron_sizes = {i: np.median([len(intron) for intron in intron_list[i]]) for i in intron_list}

    # read in the motifs
    motifs = sequo.read_motifs(motif_file)

    stop_motifs = [i for i in motifs if len(re.findall("(?=(TAA|TAG|TGA))", i))]
    non_stop_motifs = [i for i in motifs if i not in stop_motifs]

    motif_densities = {i: seqo.calc_motif_density(exon_list[i], motifs) for i in exon_list}
    stop_motif_densities = {i: seqo.calc_motif_density(exon_list[i], stop_motifs) for i in exon_list}
    non_stop_motif_densities = {i: seqo.calc_motif_density(exon_list[i], non_stop_motifs) for i in exon_list}

    # group by family
    families = gen.read_many_fields(families_file, "\t")
    motif_densities = sequo.group_family_results(motif_densities, families)
    stop_motif_densities = sequo.group_family_results(stop_motif_densities, families)
    non_stop_motif_densities = sequo.group_family_results(non_stop_motif_densities, families)
    exon_sizes = sequo.group_family_results(exon_sizes, families)
    intron_sizes = sequo.group_family_results(intron_sizes, families)

    with open(output_file, "w") as outfile:
        outfile.write("id,exon_size,intron_size,ese_density,stop_ese_density,non_stop_ese_density\n")
        for id in motif_densities:
            if id in intron_sizes:
                outfile.write("{0},{1},{2},{3},{4},{5}\n".format(id, np.median(exon_sizes[id]), np.median(intron_sizes[id]), np.median(motif_densities[id]), np.median(stop_motif_densities[id]), np.median(non_stop_motif_densities[id])))


def calc_seq_hits(coding_exons_fasta, cds_fasta, output_file, motif_file, motif_simulations_directory, required_simulations = None, families_file = None):

    # get exons
    names, seqs = gen.read_fasta(coding_exons_fasta)
    seq_list = collections.defaultdict(lambda: [])
    [seq_list[name.split(".")[0]].append(seqs[i]) for i, name in enumerate(names)]
    seq_list = sequo.pick_random_family_member(families_file, seq_list)
    # get cds
    cds_names, cds_seqs = gen.read_fasta(cds_fasta)
    cds_list = {name.split(".")[0]: cds_seqs[i] for i, name in enumerate(cds_names)}

    exon_list = []
    [exon_list.extend(seq_list[i]) for i in seq_list]
    exon_start_indices = []
    for id in seq_list:
        cds = cds_list[id]
        exon_start_indices.extend([cds.index(exon) for exon in seq_list[id]])

    filelist = {"real": motif_file}
    if required_simulations and required_simulations > 0:
        sim_files = sims = gen.get_filepaths(motif_simulations_directory)[:required_simulations]
        for i, file in enumerate(sim_files):
            filelist["sim_{0}".format(i+1)] = file

    args = [filelist, exon_list, exon_start_indices]
    outputs = simoc.run_simulation_function(list(filelist), args, sequo.calc_hits, sim_run = False)

    with open(output_file, "w") as outfile:
        outfile.write("id,stop_0,stop_1,stop_2,stop_total,non_stop_0,non_stop_1,non_stop_2,non_stop_total,stops_function_0,stops_function_1,stops_function_2,total_stop_motifs,total_non_stop_motifs\n")
        for id in outputs:
            output = outputs[id]
            outline = [id]
            [outline.extend([o[i] for i in sorted(o)] + [sum(o.values())]) for o in output[:2]]
            # outline.append(",")
            outline.extend([output[2][i] for i in sorted(output[2])])
            outline.append(output[3])
            outline.append(output[4])
            outfile.write("{0}\n".format(",".join(gen.stringify(outline))))

def calc_seq_hits_linc(exons_fasta, output_file, motif_file, motif_simulations_directory, required_simulations = None, families_file = None):

    # get exons
    names, seqs = gen.read_fasta(exons_fasta)
    seq_list = collections.defaultdict(lambda: [])
    [seq_list[name.split(".")[0]].append(seqs[i]) for i, name in enumerate(names)]
    seq_list = sequo.pick_random_family_member(families_file, seq_list)
    # seq_list = {i: seq_list[i] for i in seq_list}

    filelist = {"real": motif_file}
    if required_simulations and required_simulations > 0:
        sim_files = sims = gen.get_filepaths(motif_simulations_directory)[:required_simulations]
        for i, file in enumerate(sim_files):
            filelist["sim_{0}".format(i+1)] = file

    args = [filelist, seq_list]
    outputs = simoc.run_simulation_function(list(filelist), args, sequo.calc_hits_lincrna, sim_run = False)

    with open(output_file, "w") as outfile:
        outfile.write("id,stop_hits,non_stop_hits,norm_stop,norm_non_stop,stop_motifs,non_stop_motifs\n")
        for id in outputs:
            output = outputs[id]
            outfile.write("{0}\n".format(",".join(gen.stringify(outputs[id]))))

def process_seq_hits_linc(input_dir, output_file):
    files = gen.get_filepaths(input_dir)

    with open(output_file, "w") as outfile:
        outfile.write("run,stop_total,median_sim_stop_total,normalised_stop,stop_p,adj_stop_p,non_stop_total,median_sim_non_stop_total,normalised_non_stop,non_stop_p,adj_non_stop_p,diff,median_sim_diff,normalised_diff,diff_p,adj_diff_p\n")

        for file_no, file in enumerate(files):

            data = pd.read_csv(file)
            data["diff"] = np.divide(data["norm_stop"], data["norm_non_stop"])
            real = data.loc[data['id'] == 'real']
            sims = data.loc[data['id'] != 'real']

            norm_stops_greater = len(sims[sims["norm_stop"] >= real["norm_stop"].values[0]])
            norm_stops_p = np.divide(norm_stops_greater + 1, len(sims) + 1)
            norm_stops_adj_p = norm_stops_p*len(files)
            norm_non_stops_greater = len(sims[sims["norm_non_stop"] >= real["norm_non_stop"].values[0]])
            norm_non_stops_p = np.divide(norm_non_stops_greater + 1, len(sims) + 1)
            norm_non_stops_adj_p = norm_non_stops_p*len(files)
            diff_greater = len(sims[sims["diff"] >= real["diff"].values[0]])
            diff_p = np.divide(diff_greater + 1, len(sims) + 1)
            diff_adj_p = diff_p*len(files)

            output = [file_no+1]
            output.append(real["norm_stop"].values[0])
            output.append(sims["norm_stop"].median())
            output.append(np.divide(real["norm_stop"].values[0] - sims["norm_stop"].mean(), sims["norm_stop"].mean()))
            output.append(norm_stops_p)
            output.append(norm_stops_adj_p if norm_stops_adj_p < 1 else 1)
            output.append(real["norm_non_stop"].values[0])
            output.append(sims["norm_non_stop"].median())
            output.append(np.divide(real["norm_non_stop"].values[0] - sims["norm_non_stop"].mean(), sims["norm_non_stop"].mean()))
            output.append(norm_non_stops_p)
            output.append(norm_non_stops_adj_p if norm_non_stops_adj_p < 1 else 1)
            output.append(real["diff"].values[0])
            output.append(sims["diff"].median())
            output.append(np.divide(real["diff"].values[0] - sims["diff"].mean(), sims["diff"].mean()))
            output.append(diff_p)
            output.append(diff_adj_p if diff_adj_p < 1 else 1)
            outfile.write("{0}\n".format(",".join(gen.stringify(output))))


def process_seq_hits(input_dir, output_file):

    files = gen.get_filepaths(input_dir)
    with open(output_file, "w") as outfile:
        outfile.write("run,stop_total,median_sim_stop_total,normalised_stop,stop_p,adj_stop_p,non_stop_total,median_sim_non_stop_total,normalised_non_stop,non_stop_p,adj_non_stop_p,diff,median_sim_diff,normalised_diff,diff_p,adj_diff_p\n")

        for file_no, file in enumerate(files):

            data = pd.read_csv(file)
            new_data = pd.DataFrame()
            new_data["id"] = data["id"]
            stop_entries = ["stop_0", "stop_1", "stop_2"]
            norm_stop_entries = ["norm_stop_0", "norm_stop_1", "norm_stop_2"]
            non_stop_entries = ["non_stop_0", "non_stop_1", "non_stop_2"]
            norm_non_stop_entries = ["norm_non_stop_0", "norm_non_stop_1", "norm_non_stop_2"]
            stop_functions = ["stops_function_0", "stops_function_1", "stops_function_2"]
            for i, col in enumerate(stop_entries):
                new_data[norm_stop_entries[i]] = np.divide(data[col], data["total_stop_motifs"])*np.divide(data["total_stop_motifs"], data[stop_functions[i]])
            new_data["total_norm_stop"] = new_data["norm_stop_0"] + new_data["norm_stop_1"] + new_data["norm_stop_2"]
            for i, col in enumerate(non_stop_entries):
                new_data[norm_non_stop_entries[i]] = np.divide(data[col], data["total_non_stop_motifs"])
            new_data["total_norm_non_stop"] = new_data["norm_non_stop_0"] + new_data["norm_non_stop_1"] + new_data["norm_non_stop_2"]
            new_data["diff"] = np.divide(new_data["total_norm_stop"], new_data["total_norm_non_stop"])

            real = new_data.loc[new_data['id'] == 'real']
            sims = new_data.loc[new_data['id'] != 'real']

            norm_stops_greater = len(sims[sims["total_norm_stop"] >= real["total_norm_stop"].values[0]])
            norm_stops_p = np.divide(norm_stops_greater + 1, len(sims) + 1)
            norm_stops_adj_p = norm_stops_p*len(files)
            norm_non_stops_greater = len(sims[sims["total_norm_non_stop"] >= real["total_norm_non_stop"].values[0]])
            norm_non_stops_p = np.divide(norm_non_stops_greater + 1, len(sims) + 1)
            norm_non_stops_adj_p = norm_non_stops_p*len(files)
            diff_greater = len(sims[sims["diff"] >= real["diff"].values[0]])
            diff_p = np.divide(diff_greater + 1, len(sims) + 1)
            diff_adj_p = diff_p*len(files)

            output = [file_no+1]
            output.append(real["total_norm_stop"].values[0])
            output.append(sims["total_norm_stop"].median())
            output.append(np.divide(real["total_norm_stop"].values[0] - sims["total_norm_stop"].mean(), sims["total_norm_stop"].mean()))
            output.append(norm_stops_p)
            output.append(norm_stops_adj_p if norm_stops_adj_p < 1 else 1)
            output.append(real["total_norm_non_stop"].values[0])
            output.append(sims["total_norm_non_stop"].median())
            output.append(np.divide(real["total_norm_non_stop"].values[0] - sims["total_norm_non_stop"].mean(), sims["total_norm_non_stop"].mean()))
            output.append(norm_non_stops_p)
            output.append(norm_non_stops_adj_p if norm_non_stops_adj_p < 1 else 1)
            output.append(real["diff"].values[0])
            output.append(sims["diff"].median())
            output.append(np.divide(real["diff"].values[0] - sims["diff"].mean(), sims["diff"].mean()))
            output.append(diff_p)
            output.append(diff_adj_p if diff_adj_p < 1 else 1)
            outfile.write("{0}\n".format(",".join(gen.stringify(output))))



def calc_overlap_potential(input_fasta, motif_file, output_file, required_simulations = None, controls_directory = None, families_file = None):

    # get the exons
    names, exons = gen.read_fasta(input_fasta)
    exon_list = collections.defaultdict(lambda: [])
    [exon_list[name.split(".")[0]].append(exons[i]) for i, name in enumerate(names)]
    exon_list = sequo.pick_random_family_member(families_file, exon_list)
    exons = []
    [exons.extend(exon_list[i]) for i in exon_list]

    filelist = {"real": motif_file}
    if controls_directory and required_simulations:
        sim_files = gen.get_filepaths(controls_directory)
        for i, file in enumerate(sim_files[:required_simulations]):
            filelist["sim{0}".format(i+1)] = file

    args = [filelist, exons]
    outputs = simoc.run_simulation_function(list(filelist), args, sequo.analyse_overlaps, sim_run = False)

    with open(output_file, "w") as outfile:
        outfile.write("id,stop_overlaps,non_stop_overlaps,diff\n")
        [outfile.write("{0},{1}\n".format(id, ",".join(gen.stringify(outputs[id])))) for id in outputs]

def calc_overlap_diffs(input_fasta, motif_file, output_file, required_simulations = None, controls_directory = None, families_file = None):

    # get the exons
    names, exons = gen.read_fasta(input_fasta)
    exon_list = collections.defaultdict(lambda: [])
    [exon_list[name.split(".")[0]].append(exons[i]) for i, name in enumerate(names)]
    exon_list = sequo.pick_random_family_member(families_file, exon_list)
    exons = []
    [exons.extend(exon_list[i]) for i in exon_list]

    filelist = {"real": motif_file}
    if controls_directory and required_simulations:
        sim_files = gen.get_filepaths(controls_directory)
        for i, file in enumerate(sim_files[:required_simulations]):
            filelist["sim{0}".format(i+1)] = file

    args = [filelist, exons]
    outputs = simoc.run_simulation_function(list(filelist), args, sequo.calc_sequence_overlap_diffs, sim_run = False)

    with open(output_file, "w") as outfile:
        outfile.write("id,single_stop_density,overlap_stop_density,motif_diff,sequence_diff,diff_diff\n")
        [outfile.write("{0},{1}\n".format(id, ",".join(gen.stringify(outputs[id])))) for id in outputs]
