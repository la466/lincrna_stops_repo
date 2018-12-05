import generic as gen
import containers as cont
import file_ops as fo
import seq_ops as seqo
import sequence_ops as sequo
import sim_ops_containers as simoc
import time
import random
import numpy as np
import collections
import zipfile
import os
import multiprocessing as mp

def compare_stop_density(exons_fasta, introns_fasta, output_file, families_file = None):

    exon_names, exon_seqs = gen.read_fasta(exons_fasta)
    exon_list = {name.split("(")[0]: exon_seqs[i] for i, name in enumerate(exon_names)}

    transcript_exon_list = collections.defaultdict(lambda: [])
    [transcript_exon_list[name.split(".")[0]].append(exon_list[name]) for name in exon_list]

    transcript_exon_densities = {name: seqo.calc_seqs_stop_density(transcript_exon_list[name]) for name in transcript_exon_list}

    intron_names, intron_seqs = gen.read_fasta(introns_fasta)
    intron_list = {name.split("(")[0]: intron_seqs[i] for i, name in enumerate(intron_names)}

    transcript_intron_list = collections.defaultdict(lambda: [])
    [transcript_intron_list[name.split(".")[0]].append(intron_list[name]) for name in intron_list]

    transcript_intron_densities =  {name: seqo.calc_seqs_stop_density(transcript_intron_list[name]) for name in transcript_intron_list}
    transcript_intron_densities_min_removed =  {}
    transcript_intron_densities_max_removed =  {}
    for id in transcript_intron_list:
        transcript_intron_densities_min_removed[id] = seqo.calc_intron_seqs_stop_density(transcript_intron_list[id], remove_min = True)
        transcript_intron_densities_max_removed[id] = seqo.calc_intron_seqs_stop_density(transcript_intron_list[id], remove_max = True)

    exons_gc = {name: seqo.calc_gc_seqs_combined(transcript_exon_list[name]) for name in transcript_exon_list}
    introns_gc = {name: seqo.calc_gc_seqs_combined(transcript_intron_list[name]) for name in transcript_intron_list}

    if families_file:
        families = gen.read_many_fields(families_file, "\t")
        transcript_exon_densities = sequo.group_family_results(transcript_exon_densities, families)
        transcript_intron_densities = sequo.group_family_results(transcript_intron_densities, families)
        transcript_intron_densities_min_removed = sequo.group_family_results(transcript_intron_densities_min_removed, families)
        transcript_intron_densities_max_removed = sequo.group_family_results(transcript_intron_densities_max_removed, families)
        exons_gc = sequo.group_family_results(exons_gc, families)
        introns_gc = sequo.group_family_results(introns_gc, families)

    with open(output_file, "w") as outfile:
        outfile.write("id,exon_gc,intron_gc,exon_density,intron_density,intron_density_min_removed,intron_density_max_removed\n")
        for id in transcript_exon_densities:
            if id in transcript_intron_densities:
                args = [id, np.median(exons_gc[id]), np.median(introns_gc[id]), np.median(transcript_exon_densities[id]), np.median(transcript_intron_densities[id]), np.median(transcript_intron_densities_min_removed[id]), np.median(transcript_intron_densities_max_removed[id])]
                outfile.write("{0}\n".format(",".join(gen.stringify(args))))



def coding_exons(input_file, families_file, output_directory):

    output_directory = "{0}/coding_exons".format(output_directory)
    gen.create_output_directories(output_directory)

    # get the coding exons
    coding_exon_names, coding_exon_seqs = gen.read_fasta(input_file)
    coding_exons = {name.split("(")[0]: coding_exon_seqs[i] for i, name in enumerate(coding_exon_names) if len(coding_exon_seqs[i]) >= 210}

    coding_exon_five_prime = {i: coding_exons[i][:69] for i in coding_exons}
    coding_exon_three_prime = {i: coding_exons[i][-69:] for i in coding_exons}

    coding_exon_cores = {}
    for i in coding_exons:
        middle = int(len(coding_exons[i]) / 2)
        coding_exon_cores[i] = coding_exons[i][middle-34:middle+35]

    five_gc = collections.defaultdict(lambda: [])
    three_gc = collections.defaultdict(lambda: [])
    core_gc = collections.defaultdict(lambda: [])

    [five_gc[i.split(".")[0]].append(seqo.calc_seq_gc(coding_exon_five_prime[i])) for i in coding_exon_five_prime]
    [three_gc[i.split(".")[0]].append(seqo.calc_seq_gc(coding_exon_three_prime[i])) for i in coding_exon_three_prime]
    [core_gc[i.split(".")[0]].append(seqo.calc_seq_gc(coding_exon_cores[i])) for i in coding_exon_cores]

    families = gen.read_many_fields(families_file, "\t")

    # group the densities by transcript id
    coding_exon_five_prime_list = collections.defaultdict(lambda: [])
    coding_exon_three_prime_list = collections.defaultdict(lambda: [])
    coding_exon_core_list = collections.defaultdict(lambda: [])

    [coding_exon_five_prime_list[i.split(".")[0]].append(seqo.calc_stop_density(coding_exon_five_prime[i])) for i in coding_exon_five_prime]
    [coding_exon_three_prime_list[i.split(".")[0]].append(seqo.calc_stop_density(coding_exon_three_prime[i])) for i in coding_exon_three_prime]
    [coding_exon_core_list[i.split(".")[0]].append(seqo.calc_stop_density(coding_exon_cores[i])) for i in coding_exon_cores]

    coding_exon_five_prime_list = {i: np.mean(coding_exon_five_prime_list[i]) for i in coding_exon_five_prime_list}
    coding_exon_three_prime_list = {i: np.mean(coding_exon_three_prime_list[i]) for i in coding_exon_three_prime_list}
    coding_exon_core_list = {i: np.mean(coding_exon_core_list[i]) for i in coding_exon_core_list}

    # group the densities by family
    coding_exon_five_prime = sequo.group_family_results(coding_exon_five_prime_list, families)
    coding_exon_three_prime = sequo.group_family_results(coding_exon_three_prime_list, families)
    coding_exon_cores = sequo.group_family_results(coding_exon_core_list, families)

    five_gc = {i: np.mean(five_gc[i]) for i in five_gc}
    three_gc = {i: np.mean(three_gc[i]) for i in three_gc}
    core_gc = {i: np.mean(core_gc[i]) for i in core_gc}



    five_gc = sequo.group_family_results(five_gc, families)
    three_gc = sequo.group_family_results(three_gc, families)
    core_gc = sequo.group_family_results(core_gc, families)

    output_file = "{0}/exon_region_densities.csv".format(output_directory)
    with open(output_file, "w") as outfile:
        outfile.write("id,start_density,core_density,end_density,start_gc,core_gc,end_gc\n")
        for i in sorted(coding_exon_five_prime):
            outfile.write("{0},{1},{2},{3},{4},{5},{6}\n".format(i, np.mean(coding_exon_five_prime[i]), np.mean(coding_exon_cores[i]), np.mean(coding_exon_three_prime[i]), np.mean(five_gc[i]), np.mean(core_gc[i]), np.mean(three_gc[i])))


# def position_ds(cds_fasta, ortholog_cds_fasta, ortholog_transcript_links, output_directory):
#
#     output_directory = "{0}/position_ds".format(output_directory)
#     gen.create_output_directories(output_directory)
#
#     # get a list of the links
#     links = {i[0]: i[1] for i in gen.read_many_fields(ortholog_transcript_links, "\t")}
#
#     cds_names, cds_seqs = gen.read_fasta(cds_fasta)
#     cds_list = {name: cds_seqs[i] for i, name in enumerate(cds_names[:5])}
#
#     ortholog_names, ortholog_seqs = gen.read_fasta(ortholog_cds_fasta)
#     ortholog_cds_list = {links[name]: ortholog_seqs[ortholog_names.index(links[name])] for name in cds_list if links[name] in ortholog_names}
#
#     cds_pairs = {cds_list[id]: ortholog_cds_list[links[id]] for id in cds_list}
#
#     sequo.extract_aligment_sequence_parts(cds_pairs)
#
#
#     # # TODO: get all positions in each cds where if there was a 1 nt mutation at a 4 fold degenerate site, it would create a stop
#     # # check_conservation(human.cds_fasta, ortholog.cds_fasta, orthologs_transcript_links_file, human_ids_after_conservation_file, human_cds_after_ortholog_filter_fasta, max_dS_threshold = 0.2, max_omega_threshold = 0.5, clean_run = clean_run)
#     # # sequo.get_cds_parts_close_to_stop(cds_file, ortholog_fasta, ortholog_transcript_links, 1)
#     # # TODO: concatenate all the sequences together
#     # # TODO: calculate ds for the sequence parts
#     # # TODO: do the same for those needing 2 mutations


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
        exons = [i for i in exon_list if id in i]
        if len(exons):
            print("(W{0}) {1}/{2}: {3}".format(mp.current_process().name.split("-")[-1], i+1, len(id_list), id))
            cds_seq = cds_list[id]
            # get the start index and length of each exon in the cds
            exon_info = [[cds_seq.index(exon_list[i]), len(exon_list[i])] for i in exons]
            # read in the simulated cds seqs
            sim_cds_seqs = gen.read_many_fields(filelist[id], ",")[0]
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
