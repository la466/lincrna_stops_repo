import generic as gen
import seq_ops as seqo
import sequence_ops as sequo
import sim_ops_containers as soc
import main_tests_ops as mto
import ops
import numpy as np
import collections
import re
import os
import pandas as pd
# import scipy.stats as stats
# import matplotlib.pyplot as plt
import conservation
import os
import zipfile
from useful_motif_sets import stops, codon_map
from regex_patterns import codon_pattern
import containers as cont
import itertools as it
from Bio import SeqIO
import scipy.stats




motifs = "source_data/motif_sets/int3.txt"
# motifs = "source_data/motif_sets/RESCUE.txt"
# seqs_file = "clean_run/lincrna/lincRNA.multi_exon.exons.fasta"
# seqs_file = "clean_run/genome_sequences/lincrna/cabili/multi_exons.fasta"

families_file = "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt"
exon_file = "clean_run/genome_sequences/lincrna/cabili/exon_alignments.fasta"
intron_file = "clean_run/genome_sequences/lincrna/cabili/intron_alignments.fasta"

exon_list = collections.defaultdict(lambda: collections.defaultdict())
exon_names, exon_seqs = gen.read_fasta(exon_file)
for i, name in enumerate(exon_names):
    alignments = exon_seqs[i].split(",")
    if len(alignments[0]) > 207:
        exon_list[name.split(".")[0]][int(name.split(".")[-1])] = alignments

intron_list = collections.defaultdict(lambda: collections.defaultdict())
intron_names, intron_seqs = gen.read_fasta(intron_file)
for i, name in enumerate(intron_names):
    alignments = intron_seqs[i].split(",")
    id = name.split(".")[0]
    if len(alignments[0]) and id in exon_list:
        intron_list[id][int(name.split(".")[-1])] = alignments

# now retain only the exons with retained introns
exon_list = {id: exon_list[id] for id in exon_list if id in intron_list}
exon_list = sequo.pick_random_family_member(families_file, exon_list)

with open("temp_files/core_rates.csv", "w") as outfile:
    outfile.write("id,exon_flank,exon_core,intron_core,flank_core,core_core,,stop_ef_rate,stop_ec_rate,stop_ic_rate\n")

    h_exon_flank = []
    m_exon_flank = []
    h_exon_core = []
    m_exon_core = []
    h_intron_core = []
    m_intron_core = []
    hs_exon_flank = []
    ms_exon_flank = []
    hs_exon_core = []
    ms_exon_core = []
    hs_intron_core = []
    ms_intron_core = []

    for i, id in enumerate(exon_list):
        if i < 500:
            human_exon_flanks = []
            human_exon_cores = []
            mac_exon_flanks = []
            mac_exon_cores = []
            human_intron_cores = []
            mac_intron_cores = []

            for intron_id in intron_list[id]:
                human_intron = intron_list[id][intron_id][0]
                mac_intron = intron_list[id][intron_id][1]
                midpoint = int(len(human_intron) / 2)

                # human_core = human_intron[midpoint-33:midpoint+34]
                # mac_core = mac_intron[midpoint-33:midpoint+34]
                human_core = human_intron[20:len(human_intron)-20]
                mac_core = mac_intron[20:len(mac_intron)-20]

                human_intron_cores.append(human_core)
                mac_intron_cores.append(mac_core)

            joined_human_intron_cores = "".join(human_intron_cores)
            joined_mac_intron_cores = "".join(mac_intron_cores)

            if len(joined_human_intron_cores) > 100:
                for exon_id in exon_list[id]:
                    human_exon = exon_list[id][exon_id][0]
                    mac_exon = exon_list[id][exon_id][1]

                    human_flanks = [human_exon[2:69], human_exon[-69:-2]]
                    mac_flanks = [mac_exon[2:69], mac_exon[-69:-2]]

                    human_exon_flanks.extend(human_flanks)
                    mac_exon_flanks.extend(mac_flanks)


                    midpoint = int(len(human_exon) / 2)
                    human_core = human_exon[69:-69]
                    mac_core = mac_exon[69:-69]

                    human_exon_cores.append(human_core)
                    mac_exon_cores.append(mac_core)


            joined_human_exon_flanks = "".join(human_exon_flanks)
            joined_human_exon_cores = "".join(human_exon_cores)
            joined_mac_exon_flanks = "".join(mac_exon_flanks)
            joined_mac_exon_cores = "".join(mac_exon_cores)

            if len(joined_human_exon_flanks) > 100 and len(joined_human_exon_cores) > 100:

                exon_flank_rate = sequo.get_sub_rate(joined_human_exon_flanks, joined_mac_exon_flanks)
                exon_core_rate = sequo.get_sub_rate(joined_human_exon_cores, joined_mac_exon_cores)
                intron_core_rate = sequo.get_sub_rate(joined_human_intron_cores, joined_mac_intron_cores)

                exon_flank_exon_core = np.divide(exon_flank_rate, exon_core_rate)
                exon_core_intron_core = np.divide(exon_core_rate, intron_core_rate)

                h_ef_overlaps = sequo.sequence_overlap_indicies(joined_human_exon_flanks, stops)
                h_ef_non_overlaps = [i for i in range(len(joined_human_exon_flanks)) if i not in h_ef_overlaps]
                h_ef_stops = "".join([joined_human_exon_flanks[i] for i in h_ef_overlaps])
                h_ef_nstops = "".join([joined_human_exon_flanks[i] for i in h_ef_non_overlaps])
                hs_exon_flank.append(h_ef_stops)
                h_exon_flank.append(h_ef_nstops)

                m_ef_stops = "".join([joined_mac_exon_flanks[i] for i in h_ef_overlaps])
                m_ef_nstops = "".join([joined_mac_exon_flanks[i] for i in h_ef_non_overlaps])
                ms_exon_flank.append(m_ef_stops)
                m_exon_flank.append(m_ef_nstops)
                # stop_ec_rate = sequo.get_sub_rate(h_ef_stops, m_ef_stops)

                h_ec_overlaps = sequo.sequence_overlap_indicies(joined_human_exon_cores, stops)
                h_ec_non_overlaps = [i for i in range(len(joined_human_exon_cores)) if i not in h_ec_overlaps]
                h_ec_stops = "".join([joined_human_exon_cores[i] for i in h_ec_overlaps])
                h_ec_nstops = "".join([joined_human_exon_cores[i] for i in h_ec_non_overlaps])
                hs_exon_core.append(h_ec_stops)
                h_exon_core.append(h_ec_nstops)

                m_ec_stops = "".join([joined_mac_exon_cores[i] for i in h_ec_overlaps])
                m_ec_nstops = "".join([joined_mac_exon_cores[i] for i in h_ec_non_overlaps])
                ms_exon_core.append(m_ec_stops)
                m_exon_core.append(m_ec_nstops)
                # stop_ef_rate = sequo.get_sub_rate(h_ec_stops, m_ec_stops)

                h_ic_overlaps = sequo.sequence_overlap_indicies(joined_human_intron_cores, stops)
                h_ec_non_overlaps = [i for i in range(len(joined_human_intron_cores)) if i not in h_ic_overlaps]
                h_ic_stops = "".join([joined_human_intron_cores[i] for i in h_ic_overlaps])
                h_ic_nstops = "".join([joined_human_intron_cores[i] for i in h_ec_non_overlaps])
                hs_intron_core.append(h_ic_stops)
                h_intron_core.append(h_ic_nstops)

                m_ic_stops = "".join([joined_mac_intron_cores[i] for i in h_ic_overlaps])
                m_ic_nstops = "".join([joined_mac_intron_cores[i] for i in h_ec_non_overlaps])
                ms_intron_core.append(m_ic_stops)
                m_intron_core.append(m_ic_nstops)
                # stop_ic_rate = sequo.get_sub_rate(h_ic_stops, m_ic_stops)

                outfile.write("{0},{1},{2},{3},{4},{5}\n".format(id, exon_flank_rate, exon_core_rate, intron_core_rate, exon_flank_exon_core, exon_core_intron_core))

    h_exon_flank = "".join(h_exon_flank)
    m_exon_flank = "".join(m_exon_flank)
    h_exon_core = "".join(h_exon_core)
    m_exon_core = "".join(m_exon_core)
    h_intron_core = "".join(h_intron_core)
    m_intron_core = "".join(m_intron_core)
    hs_exon_flank = "".join(hs_exon_flank)
    ms_exon_flank = "".join(ms_exon_flank)
    hs_exon_core = "".join(hs_exon_core)
    ms_exon_core = "".join(ms_exon_core)
    hs_intron_core = "".join(hs_intron_core)
    ms_intron_core = "".join(ms_intron_core)

    fs = sequo.get_sub_rate(hs_exon_flank, ms_exon_flank)
    ecs = sequo.get_sub_rate(hs_exon_core, ms_exon_core)
    ics = sequo.get_sub_rate(hs_intron_core, ms_intron_core)
    f = sequo.get_sub_rate(h_exon_flank, m_exon_flank)
    ec = sequo.get_sub_rate(h_exon_core, m_exon_core)
    ic = sequo.get_sub_rate(h_intron_core, m_intron_core)

    print(f, fs, np.divide(fs, f))
    print(ec, ecs, np.divide(ecs, ec))
    print(ic, ics, np.divide(ics, ic))



                # if inverse:
                #     overlaps = [i for i in list(range(len(sequence))) if i not in overlaps]
                # chunked_overlaps = chunk_overlaps(overlaps)
                # overlap_motifs = ["".join([sequence[i] for i in chunk]) for chunk in chunked_overlaps]
                # return overlap_motifs
