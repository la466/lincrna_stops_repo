import generic as gen
import file_ops as fo
import conservation as cons
import sim_ops_containers as simopc
import seq_ops as seqo
from regex_patterns import codon_pattern, exon_number_pattern, gene_id_pattern, transcript_id_pattern
import re
import os
import collections
import random
import numpy as np
from Bio.Phylo.PAML import codeml
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import copy
import sys
from useful_motif_sets import nucleotides, stops, codon_map, twofold, fourfold, one_away_codons
from progressbar import ProgressBar
import multiprocessing as mp
import itertools as it
import time
import scipy.stats

pbar = ProgressBar()

def generate_genome_dataset(gtf_file, genome_fasta, dataset_name, dataset_output_directory, id_list = None, filter_by_transcript = None, filter_by_gene = None, filter_one_per_gene = None, clean_run = None, transcript_id_search_pattern = None, gene_id_search_pattern = None, stop_codon_set = None):

    # if we want a clean run, remove any previous output directory
    if clean_run:
        gen.remove_directory(dataset_output_directory)
    gen.create_output_directories(dataset_output_directory)

    # initalise genome object
    genome = Genome_Functions(gtf_file, genome_fasta, dataset_name, dataset_output_directory, clean_run, transcript_id_search_pattern, gene_id_search_pattern, stop_codon_set)
    # get the genome features
    genome.generate_genome_features(id_list, filter_by_transcript, filter_by_gene)
    # load the genome dataset
    genome.load_genome_dataset(id_list)
    # get the cds features
    genome.get_cds_features(filter_by_gene)
    # now build the cds sequences
    genome.build_coding_sequences()
    # filter the sequences for quality
    genome.quality_filter_sequences()
    # filter to one transcript per gene
    if filter_one_per_gene:
        genome.filter_one_transcript_per_gene()
        # have to get the cds_fasta here for the next step
        genome.cds_fasta = genome.filelist["unique_cds_fasta"]
    # get a list of genes from the cleaned transcripts
    genome.list_genes_from_transcripts()

    return genome, genome.filelist


class Genome_Functions(object):

    def __init__(self, gtf_file, genome_fasta, dataset_name, dataset_output_directory, clean_run, transcript_id_search_pattern, gene_id_search_pattern, stop_codon_set):
        self.gtf_file = gtf_file
        self.genome_fasta = genome_fasta
        self.dataset_name = dataset_name
        self.dataset_output_directory = dataset_output_directory
        self.clean_run = clean_run
        if not transcript_id_search_pattern:
            self.transcript_id_pattern = transcript_id_pattern
        else:
            self.transcript_id_pattern = transcript_id_search_pattern
        if not gene_id_search_pattern:
            self.gene_id_pattern = gene_id_pattern
        else:
            self.gene_id_pattern = gene_id_search_pattern
        if not stop_codon_set:
            self.stop_codon_set = stops
        else:
            self.stop_codon_set = stop_codon_set
        self.filelist = {}



    def build_coding_sequences(self):
        """
        Build the coding sequences
        """

        full_cds_fasta = "{0}/{1}.cds.no_filtering_full.fasta".format(self.dataset_output_directory, self.dataset_name)
        if not os.path.isfile(full_cds_fasta) or self.clean_run:
            build_coding_sequences(self.full_cds_features_bed, self.genome_fasta, full_cds_fasta)
        self.full_cds_fasta = full_cds_fasta
        self.filelist["full_cds_fasta"] = full_cds_fasta


    def filter_one_transcript_per_gene(self):
        """
        Filter sequences to only leave one transcript per gene
        """

        unique_cds_fasta = "{0}/{1}.cds.unique_gene_transcripts.step2.fasta".format(self.dataset_output_directory, self.dataset_name)
        if not os.path.isfile(unique_cds_fasta) or self.clean_run:
            filtered_transcript_ids, filtered_transcript_seqs = gen.read_fasta(self.quality_filtered_cds_fasta)
            # # get a list of transcripts from the gtf file and keep only those that we want
            features = gen.read_many_fields(self.dataset_features_bed, "\t")
            transcript_list = extract_transcript_features(features, filtered_transcript_ids)
            transcript_list = {i: transcript_list[i] for i in transcript_list if i in filtered_transcript_ids}
            # keep only the longest transcript if more than one per gene
            unique_gene_transcripts = filter_one_transcript_per_gene(transcript_list)
            transcript_list = unique_gene_transcripts
            unique_gene_cds = {name: filtered_transcript_seqs[i] for i, name in enumerate(filtered_transcript_ids) if name in transcript_list}
            # write the unique transcripts to file
            fo.write_fasta(unique_gene_cds, unique_cds_fasta)
        self.unique_cds_fasta = unique_cds_fasta
        self.filelist["unique_cds_fasta"] = unique_cds_fasta


    def get_cds_features(self, filter_by_gene = None):
        """
        Get a list of CDS features from the feature list
        """

        # get the CDS features
        full_cds_features_bed = "{0}/{1}.cds.full_features.bed".format(self.dataset_output_directory, self.dataset_name)
        if not os.path.isfile(full_cds_features_bed) or self.clean_run:
            cds_features = extract_cds_features(self.dataset_features_bed, self.ids, filter_by_gene = filter_by_gene)
            # write the cds_features to an output_file
            fo.write_features_to_bed(cds_features, full_cds_features_bed)
        self.full_cds_features_bed = full_cds_features_bed
        self.filelist["full_cds_features_bed"] = full_cds_features_bed


    def generate_genome_features(self, input_list = None, filter_by_transcript = None, filter_by_gene = None):
        """
        Get a list of genome features. If wanting to get features with specific IDs, provide
        those in a list.

        Args:
            input_list (list): if set, list of ids to filter features by
            filter_by_transcript (bool): if true, filter by transcript id
            filter_by_gene (bool): if true, filter by gene id
        """

        dataset_features_bed = "{0}/{1}.{2}_features_dataset.bed".format(self.dataset_output_directory, self.gtf_file.split("/")[-1], self.dataset_name)
        # only run if the file doesn't exist or a clean run
        if not os.path.isfile(dataset_features_bed) or self.clean_run:
            generate_genome_features_dataset(self.dataset_name, self.gtf_file, dataset_features_bed, self.transcript_id_pattern, self.gene_id_pattern, input_list = input_list, filter_by_transcript = filter_by_transcript, filter_by_gene = filter_by_gene)
        self.dataset_features_bed = dataset_features_bed
        self.filelist["dataset_features_bed"] = dataset_features_bed


    def list_genes_from_transcripts(self):
        """
        Get a list of transcripts with the associated genes

        Args:
            input_fasta (str): path to input fasta. Set here because we might not filter
                before to one transcript per gene
        """

        unique_transcript_gene_list_file = "{0}/{1}.cds.genes.bed".format(self.dataset_output_directory, self.dataset_name)
        if not os.path.isfile(unique_transcript_gene_list_file) or self.clean_run:
            get_clean_genes(self.cds_fasta, self.full_cds_features_bed, unique_transcript_gene_list_file)
        self.unique_transcript_gene_list_file = unique_transcript_gene_list_file
        self.filelist["unique_transcript_gene_list_file"] = unique_transcript_gene_list_file


    def load_genome_dataset(self, input_list = None):
        """
        Load a genome dataset from file.
        """

        try:
            entries = gen.read_many_fields(self.dataset_features_bed, "\t")
            # self.extracted_features = [[i[0], int(i[1]), int(i[2]), i[3], i[4], i[5], i[6], i[7]] for i in entries if len(i) and i[0][0] != "#"]
            if input_list:
                self.ids = input_list
            else:
                self.ids = list(set([i[3] for i in entries if "#" not in i[0]]))

            print("{0} dataset loaded...".format(self.dataset_name))

        except FileNotFoundError:
            raise FileNotFoundError


    def quality_filter_sequences(self):
        """
        Filter the sequences in the genome and remove any that dont
        pass filtering
        """

        quality_filtered_cds_fasta = "{0}/{1}.cds.quality_filtered.step1.fasta".format(self.dataset_output_directory, self.dataset_name)
        if not os.path.isfile(quality_filtered_cds_fasta) or self.clean_run:
            quality_filter_cds_sequences(self.full_cds_fasta, quality_filtered_cds_fasta, stop_codon_set = self.stop_codon_set)
        self.quality_filtered_cds_fasta = quality_filtered_cds_fasta
        self.cds_fasta = quality_filtered_cds_fasta
        self.filelist["quality_filtered_cds_fasta"] = quality_filtered_cds_fasta



class PAML_Functions(object):

    def __init__(self, input_file = None, output_file = None, tree_file = None, working_dir = None):
        self.input_file = input_file
        self.output_file = output_file
        self.tree_file = tree_file
        if not working_dir:
            working_dir = os.getcwd()
        self.working_dir = working_dir

    def cleanup(self):
        """
        Clean up any output files
        """

        cwd = os.getcwd()
        output_files = ["2NG.dN", "2NG.dS", "2NG.t", "4fold.nuc", "codeml.ctl", "lnf", "rst", "rst1", "rub"]
        output_files = ["{0}/{1}".format(cwd, i) for i in output_files]
        [gen.remove_file(i) for i in output_files]
        if cwd != self.working_dir:
            gen.remove_directory(self.working_dir)


    def generate_tree_file(self):
        """
        Generate a simple tree file for use with Biopython codeml wrapper
        """

        temp_dir = "temp_files"
        gen.create_output_directories(temp_dir)
        tree_file = "{0}/temp_tree.{1}.tree".format(temp_dir, random.random())
        with open(tree_file, "w") as outfile:
            outfile.write("2 1\n(1,2);")
        self.tree_file = tree_file


    def run_codeml(self, input_file = None, output_file = None, tree_file = None, command = None, verbose = False):
        """
        Run Biopython codeml.

        Args:
            input_file (str): if set, use specified input file
            output_file (str): if set, use specified output file
            tree_file (str): if set, use specified tree file
            command (str): if set, use the given command
            verbose (bool): if set, print output to terminal
        """

        # create a tree file if not given
        if not tree_file:
            self.generate_tree_file()

        if not input_file:
            input_file = self.input_file
        if not output_file:
            output_file = self.output_file
        # set up codeml
        cml = codeml.Codeml(alignment = input_file, out_file = output_file, working_dir = self.working_dir, tree = self.tree_file)
        cml.set_options(seqtype = 1, runmode = 0, model = 0, NSsites = [])
        if command:
            cml_dict = cml.run(command = command, verbose = True)
        else:
            cml_dict = cml.run()
        return cml_dict


def build_coding_sequences(cds_features_bed, genome_fasta, output_file):
    """
    Build CDS sequences from the list of features.

    Args:
        cds_features_bed (str): path to file containing cds features
        genome_fasta (str): path to genome fasta file
        output_file (str): path to output fasta file
    """

    print("Building cds...")

    stop_codons = ["TAA", "TAG", "TGA"]

    temp_dir = "temp_dir"
    gen.create_output_directories(temp_dir)

    # create temp file to contain sequences
    temp_file = "{0}/temp_cds_features.fasta".format(temp_dir)
    fo.fasta_from_intervals(cds_features_bed, temp_file, genome_fasta, names=True)

    # now build the sequences
    sequence_parts = collections.defaultdict(lambda: collections.defaultdict())
    sequence_names, seqs = gen.read_fasta(temp_file)
    for i, name in enumerate(sequence_names):
        transcript = name.split(".")[0]
        exon_id = int(name.split(".")[1].split("(")[0])

        seq = seqs[i]
        # set the exon_id to an arbritrarily high number if it is the annotate stop codon
        if len(seq) == 3 and seq in stop_codons:
            exon_id = 9999999
        sequence_parts[transcript][exon_id] = seq

    with open(output_file, "w") as outfile:
        for transcript in sequence_parts:
            sequence = []
            for exon_id in sorted(sequence_parts[transcript]):
                sequence.append(sequence_parts[transcript][exon_id])
            outfile.write(">{0}\n{1}\n".format(transcript, "".join(sequence)))

    gen.remove_directory(temp_dir)


def build_sequences_from_exon_fasta(input_fasta, output_fasta):
    names, seqs = gen.read_fasta(input_fasta)
    outputs = collections.defaultdict(lambda: collections.defaultdict())
    for i, name in enumerate(names):
        id = name.split(".")[0]
        exon = int(name.split(".")[1].split("(")[0])
        outputs[id][exon] = seqs[i]

    with open(output_fasta, "w") as outfile:
        for id in outputs:
            seq = []
            [seq.append(outputs[id][exon]) for exon in sorted(outputs[id])]
            outfile.write(">{0}\n{1}\n".format(id, "".join(seq)))


def calc_codon_sets_ds_wrapper(codon_sets, restricted_sequences, output_file = None, output_directory = None):
    """
    Wrapper to calculate the ds score for the sequence parts where the synonymous
    site can and cannot form on of the codons in a codon set in sequence parts that overlap
    a set of given motifs

    Args:
        motif_set_list (list): list of integers to iterate over, corresponds to motif_sets
        motif_sets (dict): dictionary containing paths to motif set files
        codon_sets (list): a list of list containing codon sets to iterate over
        sequence_alignments (dict): dictionary containing the aligned sequences for each transcript
        output_file (str): if set, output to given file
        output_directory (str): if set, save files to output directory

    Returns:
        outputs (list): list of output file paths
    """

    if not output_file and not output_directory:
        print("Please provide output file or directory...")
        raise Exception

    outputs = []

    # print(mp.current_process().name.split("-")[-1], motif_set_list)
    for i, codon_set in enumerate(codon_sets):
        print("(W{0}) {1}/{2}".format(mp.current_process().name.split("-")[-1], i+1, len(codon_sets)))
        # set up the output filepath
        if output_file:
            output_file = output_file
        if output_directory:
            output_file = "{0}/{1}.txt".format(output_directory, "_".join(sorted(codon_set)))
        outputs.append(output_file)
        # calculcate the ds score for the codon set
        calc_motif_set_codons_ds([codon_set], restricted_sequences, output_file)

    return outputs


def calc_motif_set_codons_ds(codon_sets, sequence_alignments, output_file):
    """
    Given sets of codons, 1) extract all the pieces of the sequence alignments
    where the synonymous site creates the codon in either of the +1 or +2 frames
    and for those sites that dont 2) calculate the ds score for the resulting cases

    Args:
        codon_sets (list): list of lists containing codon sets e.g. [["TAA", "TAG"], ["CAT", "TGG"]]
        sequence_alignments (dict): dict containing sequence alignments
        output_file (str): path to output file to write results to
    """


    codon_hits_ds = {}
    no_codon_hits_ds = {}

    with open(output_file, "w") as outfile:
        outfile.write("codon_set,hits_ds,no_hits_ds,hits_query_count,no_hits_query_count\n")
        # for each of the codon sets provided
        for codon_set in codon_sets:
            hits_sequences = {}
            no_hits_sequences = {}
            # for each of the sequence alignments
            for i, id in enumerate(sequence_alignments):
                # split those sequences where the synonymous site has a hit to the current codon set
                hit_seqs, no_hits_seqs = get_sequences_synonymous_hits_to_codons(sequence_alignments[id], codon_set)
                hits_sequences[i] = hit_seqs
                no_hits_sequences[i] = no_hits_seqs

            # now calculate the ds scores for each set
            hits_alignment_strings = list_alignments_to_strings(hits_sequences)
            no_hits_alignment_strings = list_alignments_to_strings(no_hits_sequences)
            # calculate the ds scores
            hits_ds = cons.calc_ds(hits_alignment_strings)
            no_hits_ds = cons.calc_ds(no_hits_alignment_strings)
            # write to file
            args = ["_".join(sorted(codon_set)), hits_ds, no_hits_ds, int(np.divide(len(hits_alignment_strings[0]), 3)), int(np.divide(len(no_hits_alignment_strings[0]), 3))]
            outfile.write("{0}\n".format(",".join(gen.stringify(args))))


def get_codon_overlaps(overlap_indices, sequence_length, inverse = None, full_set = None):

    kept_indices = []
    for i in range(0, sequence_length, 3):
        # get each codons indices
        local_range = list(range(i, i+3))
        # if one of the codons indices is in the overlap indices
        # and the synonymous site is in the overlaps
        if list(set(local_range) & set(overlap_indices)) and i+2 in overlap_indices:
            # if we are doing the inverse, i.e. all sites that didn't have the original overlap
            if inverse:
                # keep only sites that strictly do not have any overlap with the motifs
                if not list(set(local_range) & set(full_set)):
                    kept_indices.extend(local_range)
            else:
                kept_indices.extend(local_range)
    return kept_indices

def chunk_indices(flat_list):
    chunked_list = []
    # chunk = []
    for i in range(0, len(flat_list), 3):
        chunked_list.append(flat_list[i:i+3])
    #     if index + 1 not in flat_list:
    #         chunked_list.append(chunk)
    #         chunk = []
    # chunked_list.append(chunk)
    # chunked_list = [i for i in chunked_list if len(i)]
    return chunked_list


def sequence_overlap_indicies(sequence, query_motifs, list_set = True):
    """
    Given a sequence, return all indicies that one of the query motifs hit

    Args:
        sequence (str): a sequence
        query_motifs (list): list of motifs to identify in sequence
        list_set (bool): if true, return only unique positions
    Returns:
        overlap_indices (list): list of sequence indicies
    """

    overlap_indices = []
    # compile a search string that has all the query motifs
    regex_pattern = re.compile("(?=({0}))".format("|".join(query_motifs)))
    hits = re.finditer(regex_pattern, sequence)
    for hit in hits:
        overlap_indices.extend(range(hit.start(), hit.start() + len(hit.group(1))))
    if list_set:
        return sorted(list(set(overlap_indices)))
    else:
        return sorted(overlap_indices)

def sequence_motif_overlap(sequence, query_motif):
    """
    Given a sequence, return all indicies that a motif hits

    Args:
        sequence (str): a sequence
        query_motif (list): motif to search for
    Returns:
        overlap_indices (list): list of sequence indicies
    """

    overlap_indices = []
    # compile a search string that has all the query motifs
    regex_pattern = re.compile("(?=({0}))".format(query_motif))
    hits = re.finditer(regex_pattern, sequence)
    for hit in hits:
        overlap_indices.extend(range(hit.start(), hit.start() + len(hit.group(1))))
    return sorted(list(set(overlap_indices)))


def get_sequence_overlap_codons(sequence, overlap_indices):
    """
    Given a sequence, return the codons in the sequence that overlap with one of
    the overlap indices.

    Args:
        sequence (str): sequence
        overlap_indices (list): list of sequence indices
    """

    overlap_codons = []
    non_overlap_codons = []

    # for each codon start
    for i in range(0, len(sequence), 3):
        # get the codon indices
        query_values = list(range(i, i+3))
        # check if there is an overlap
        if query_values[-1] in overlap_indices:
            codon = sequence[i:i+3]
            overlap_codons.append(codon)
        else:
            non_overlap_codons.append(sequence[i:i+3])

    return overlap_codons, non_overlap_codons


def get_sequence_mutation_codons(sequences, overlap_indices):
    """
    Given a sequence, return the codons in the sequence that are one away from a
    stop codon and those that arent

    Args:
        sequence (str): sequence
        overlap_indices (list): list of sequence indices
    """

    one_away_codons = [[], []]
    other_codons = [[], []]

    plus_1_one_away = ["TCA", "TTA", "TCG", "TGG", "TTG"]
    plus_2_one_away = ["AAA", "CAA", "GAA", "AAG", "CAG", "GAG", "AGA", "CGA", "GGA"]

    # for each codon start
    for i in range(0, len(sequences[0]), 3):
        # get the codon indices
        query_values = list(range(i, i+3))
        # check if there is an overlap
        if query_values[-1] in overlap_indices:
            # in the case where only the next nucleotide of a motif is ESE
            if i+3 in overlap_indices and i+4 not in overlap_indices:
                test_codon_seq_1 = sequences[0][i+1:i+4]
                test_codon_seq_2 = sequences[1][i+1:i+4]
                if test_codon_seq_1 in plus_1_one_away or test_codon_seq_2 in plus_1_one_away:
                    one_away_codons[0].append(sequences[0][i:i+3])
                    one_away_codons[1].append(sequences[1][i:i+3])
            elif i+3 in overlap_indices and i+4 in overlap_indices:
                test_codon_1_seq_1 = sequences[0][i+1:i+4]
                test_codon_2_seq_1 = sequences[0][i+2:i+5]
                test_codon_1_seq_2 = sequences[1][i+1:i+4]
                test_codon_2_seq_2 = sequences[1][i+2:i+5]
                if test_codon_1_seq_1 in plus_1_one_away or test_codon_2_seq_1 in plus_1_one_away or test_codon_1_seq_2 in plus_2_one_away or test_codon_2_seq_2 in plus_2_one_away:
                    one_away_codons[0].append(sequences[0][i:i+3])
                    one_away_codons[1].append(sequences[1][i:i+3])
            else:
                other_codons[0].append(sequences[0][i:i+3])
                other_codons[1].append(sequences[1][i:i+3])

    return one_away_codons, other_codons


def get_sequence_overlap_stop_codons(sequences, overlaps):

    stops_overlaps = [[], []]
    non_stops_overlaps = [[], []]

    # for each codon start
    for i in range(0, len(sequences[0]), 3):
        codon_indices = list(range(i, i+3))
        if codon_indices[-1] in overlaps:
            query1 = sequences[0][i:i+6]
            query2 = sequences[1][i:i+6]
            if query1[1:4] in stops or query1[2:5] in stops or query2[1:4] in stops or query2[2:5] in stops:
                stops_overlaps[0].append(sequences[0][i:i+3])
                stops_overlaps[1].append(sequences[1][i:i+3])
            else:
                non_stops_overlaps[0].append(sequences[0][i:i+3])
                non_stops_overlaps[1].append(sequences[1][i:i+3])
    return stops_overlaps, non_stops_overlaps



def get_sequence_non_overlap_stop_codons(sequences, overlaps):

    stops_overlaps = [[], []]
    non_stops_overlaps = [[], []]

    # for each codon start
    for i in range(0, len(sequences[0]), 3):
        codon_indices = list(range(i, i+3))
        if codon_indices[-1] not in overlaps:
            query1 = sequences[0][i:i+6]
            query2 = sequences[1][i:i+6]
            if query1[1:4] in stops or query1[2:5] in stops or query2[1:4] in stops or query2[2:5] in stops:
                stops_overlaps[0].append(sequences[0][i:i+3])
                stops_overlaps[1].append(sequences[1][i:i+3])
            else:
                non_stops_overlaps[0].append(sequences[0][i:i+3])
                non_stops_overlaps[1].append(sequences[1][i:i+3])
    return stops_overlaps, non_stops_overlaps

def extract_motif_codons_from_sequence(sequences, query_list, get_stops = None, query_overlaps = None, only_both_hits = None):
    """
    Given a set of sequences, get the codons that the query motifs overlaps.

    Args:
        sequences (list): list of sequences
        query_list (list)/(str): list of motifs to query, or filepath with overlap indices

    Returns:
        overlap_codons (list): list of codons that overlap motifs
        non_overlap_codons (list): list of codons for each seq not overlapping
    """

    # use this here to pass the given overlap list
    if query_overlaps:
        overlaps = query_list
    else:
        # # get a list of overlaps in any of the sequences
        # overlaps = sorted(list(set(gen.flatten([sequence_overlap_indicies(seq, query_list) for seq in sequences]))))

        sequence1 = sequences[0]
        sequence2 = sequences[1]
        hits = sequence_overlap_indicies(sequence1, query_list)
        chunked_hits = chunk_overlaps(hits)

        #all hits
        overlaps = hits

        # if only_both_hits:
        #     seq1 = sequences[0]
        #     seq2 = sequences[1]
        #
        #     #get the overlaps and chunk into motif overlaps for motifs that are an ese in both cases
        #     retained_hits = []
        #     for motif in query_list:
        #         # predict hits to each motif
        #         hits = sequence_motif_overlap(seq1, motif)
        #         chunked_hits = chunk_overlaps(hits)
        #         # get the motifs
        #         hit_seq1_motifs = ["".join([seq1[i] for i in chunk]) for chunk in chunked_hits]
        #         hit_seq2_motifs = ["".join([seq2[i] for i in chunk]) for chunk in chunked_hits]
        #         # ask whether they are both in the motifs
        #         for j, motif in enumerate(hit_seq1_motifs):
        #             if motif in query_list and hit_seq2_motifs[j] in query_list:
        #                 retained_hits.extend(chunked_hits[j])
        #
        #     # get a set of retained hits and chunk them
        #     overlaps = list(set(retained_hits))

    # now get the codons from these overlaps
    # note this assumes that the sequences are of length 3
    overlap_codons = []
    non_overlap_codons = []
    # stop_overlap_codons = []
    # non_stop_overlap_codons = []


    for sequence in sequences:
        overlapped, non_overlapped = get_sequence_overlap_codons(sequence, overlaps)
        overlap_codons.append(overlapped)
        non_overlap_codons.append(non_overlapped)
        # stop_overlap_codons.append(stop_overlapped)
        # non_stop_overlap_codons.append(non_stop_overlapped)

    if get_stops:
        stop_overlap_codons, non_stop_overlap_codons = get_sequence_overlap_stop_codons(sequences, overlaps)
        return overlap_codons, non_overlap_codons, stop_overlap_codons, non_stop_overlap_codons
    else:
        return overlap_codons, non_overlap_codons


def extract_non_motif_codons_from_sequence(sequences, query_motifs, get_stops = None):
    """
    Given a set of sequences, get the codons that the query motifs overlaps.

    Args:
        sequences (list): list of sequences
        query_motifs (list): list of motifs to query

    Returns:
        overlap_codons (list): list of codons that overlap motifs
        non_overlap_codons (list): list of codons for each seq not overlapping
    """

    # get a list of overlaps in any of the sequences
    overlaps = sorted(list(set(gen.flatten([sequence_overlap_indicies(seq, query_motifs) for seq in sequences]))))
    # now get the codons from these overlaps
    # note this assumes that the sequences are of length 3
    non_overlap_stop_codons, non_overlap_non_stop_codons = get_sequence_non_overlap_stop_codons(sequences, overlaps)
    return non_overlap_stop_codons, non_overlap_non_stop_codons


def extract_motif_codons_mutations_from_sequence(sequences, query_motifs, get_stops = None):
    """
    Given a set of sequences, get the codons that the query motifs overlaps.

    Args:
        sequences (list): list of sequences
        query_motifs (list): list of motifs to query

    Returns:
        overlap_codons (list): list of codons that overlap motifs
        non_overlap_codons (list): list of codons for each seq not overlapping
    """

    # get a list of overlaps in any of the sequences
    overlaps = sorted(list(set(gen.flatten([sequence_overlap_indicies(seq, query_motifs) for seq in sequences]))))
    # now get the codons from these overlaps
    # note this assumes that the sequences are of length 3
    # one_away_codons = []
    # other_codons = []
    #
    # for sequence in sequences:
    one_away, other_codons = get_sequence_mutation_codons(sequences, overlaps)
        # one_away_codons.append(one_away)
        # other_codons.append(other_codons)

    return one_away, other_codons


def extract_motif_codons_from_alignments(alignment_sequences, query_list, get_stops = None, query_overlaps = None, only_both_hits = None):
    """
    Wrapper to extract the codons that contain motif hits from alignments sequences.
    Given a list of sequences, extract the parts that overlap the given motifs

    Args:
        alignment_sequences (dict): dictionary containing sequences
        query_list (list)/(str): list of motifs to query, or a filepath with overlap indices

    Returns:
        overlap_codons (list): list of codons that hit one of the motifs
        non_overlap_codons (list): list of codons that didnt hit one of the motifs
    """

    overlap_codons = []
    non_overlap_codons = []

    if get_stops:
        stop_overlap_codons = []
        non_stop_overlap_codons = []

    for id in alignment_sequences:
        if isinstance(query_list, dict):
            id_query_list = query_list[id]
        else:
            id_query_list = query_list

        if get_stops:

            overlap, non_overlap, stop_overlap, non_stop_overlap = extract_motif_codons_from_sequence(alignment_sequences[id], id_query_list, get_stops = True, query_overlaps = query_overlaps, only_both_hits = only_both_hits)
            overlap_codons.append(overlap)
            non_overlap_codons.append(non_overlap)
            stop_overlap_codons.append(stop_overlap)
            non_stop_overlap_codons.append(non_stop_overlap)

        else:
            overlap, non_overlap = extract_motif_codons_from_sequence(alignment_sequences[id], id_query_list, get_stops = False, query_overlaps = query_overlaps, only_both_hits = only_both_hits)
            overlap_codons.append(overlap)
            non_overlap_codons.append(non_overlap)


    if get_stops:
        return overlap_codons, non_overlap_codons, stop_overlap_codons, non_stop_overlap_codons
    else:
        return overlap_codons, non_overlap_codons

def extract_non_ese_stop_alignments(alignment_sequences, query_motifs):
    """
    Wrapper to extract the codons that contain motif hits from alignments sequences.
    Given a list of sequences, extract the parts that overlap the given motifs

    Args:
        alignment_sequences (dict): dictionary containing sequences
        query_motifs (list): list of motifs to query

    Returns:
        overlap_codons (list): list of codons that hit one of the motifs
        non_overlap_codons (list): list of codons that didnt hit one of the motifs
    """

    non_overlap_stop_codons = []
    non_overlap_non_stop_codons = []

    for id in alignment_sequences:
        stop_overlap, non_stop_overlap = extract_non_motif_codons_from_sequence(alignment_sequences[id], query_motifs)
        non_overlap_stop_codons.append(stop_overlap)
        non_overlap_non_stop_codons.append(non_stop_overlap)
    return non_overlap_stop_codons, non_overlap_non_stop_codons

def extract_motif_codons_mutations_from_alignments(alignment_sequences, query_motifs):
    """
    Wrapper to extract the codons that contain motif hits from alignments sequences.
    Given a list of sequences, extract the parts that overlap the given motifs

    Args:
        alignment_sequences (dict): dictionary containing sequences
        query_motifs (list): list of motifs to query

    Returns:
        # overlap_codons (list): list of codons that hit one of the motifs
        # non_overlap_codons (list): list of codons that didnt hit one of the motifs
    """

    one_away_codons = []
    other_codons = []

    for id in alignment_sequences:
        one_away, others = extract_motif_codons_mutations_from_sequence(alignment_sequences[id], query_motifs, get_stops = False)
        one_away_codons.append(one_away)
        other_codons.append(others)
    return one_away_codons, other_codons

def extract_motif_sequences_from_alignment(alignment_seqs, motif_set, reverse = None, only_fourfold = True, codon_set = None):
    """
    Keep anything that looks like it belongs in the motif set
    from either of the alignment sequences

    Args:
        alignment_seqs (list): list containing [seq1, seq2] of aligned sequences
        motif_set (list): list of motifs to query the alignment sequences.
        If the motif sits in frame, keep all sites, but if the next motif hit
        isn't straight away, add a buffer of NNN to prevent any extra codons being created.
        If the motif sits out of frame, if the first or first two are hits, then it is the last
        codon of the motif. Take the first two nucleotides and add N to the end to remove the
        synonymous site of the codon. If the last two or last nucleotides are hits, it is the
        first codon of the motif. Take the last two nucleotides and add N to keep the synonymous
        site. Then append all codons together.
        only_fourfold (bool): if true, only use fourfold degenerates
    Returns:
        remaining_motif_sequences (list): list containing [seq1, seq1] but only
        sites that overlap one of the motifs
    """

    retained_overlap_sequences = {}
    retained_non_overlap_sequences = {}
    retained_stop_overlap_sequences = {}
    retained_non_stop_overlap_sequences = {}


    for i,id in enumerate(alignment_seqs):
        if i:
            alignment_set = alignment_seqs[id]
            # get a list of all indices of all positions that overlap with something
            # that looks like a motif in the set and all the positions that dont
            overlap_indices = sorted(get_motifs_overlap_indices(alignment_set, motif_set))
            non_overlap_indices = [i for i in range(0, len(alignment_set[0])) if i not in overlap_indices]

            # now get the codon indices that correspond to the overlap sets
            # only want those codons that have an overlap at a synonymous site
            overlap_codon_indices = get_codon_overlaps(overlap_indices, len(alignment_set[0]))
            # for the indicies that don't overlap, we only want to keep those that strictly arent overlaps
            non_overlap_codon_indices = get_codon_overlaps(non_overlap_indices, len(alignment_set[0]), inverse = True, full_set = overlap_indices)

            # now need to check if the codons are retained if we are checking fourfold
            # needs to be done here because we want a codon to be fourfold
            if only_fourfold:
                retained_overlaps = []
                retained_non_overlaps = []
                overlap_codon_indices_chunked = [overlap_codon_indices[i:i+3] for i in list(range(0, len(overlap_codon_indices), 3))]
                for codon_index in overlap_codon_indices_chunked:
                    codon1 = alignment_set[0][codon_index[0]:codon_index[-1]+1]
                    codon2 = alignment_set[1][codon_index[0]:codon_index[-1]+1]
                    if codon1 in fourfold or codon2 in fourfold:
                        retained_overlaps.extend(codon_index)
                non_overlap_codon_indices_chunked = [non_overlap_codon_indices[i:i+3] for i in list(range(0, len(non_overlap_codon_indices), 3))]
                for codon_index in non_overlap_codon_indices_chunked:
                    codon1 = alignment_set[0][codon_index[0]:codon_index[-1]+1]
                    codon2 = alignment_set[1][codon_index[0]:codon_index[-1]+1]
                    if codon1 in fourfold or codon2 in fourfold:
                        retained_non_overlaps.extend(codon_index)

                overlap_codon_indices = retained_overlaps
                non_overlap_codon_indices = retained_non_overlaps

            # now group the indices into chunks, based on whether n + 1 exists
            overlap_codon_indices = chunk_indices(overlap_codon_indices)
            non_overlap_codon_indices = chunk_indices(non_overlap_codon_indices)


            overlap_sequences = [[sequence[i[0]:i[-1]+1] if i != "-" else "---" for i in overlap_codon_indices] for sequence in alignment_set]
            non_overlap_sequences = [[sequence[i[0]:i[-1]+1] if i != "-" else "---" for i in non_overlap_codon_indices] for sequence in alignment_set]
            # join the sequences
            motif_overlap_sequences = ["".join(i) for i in overlap_sequences]
            motif_non_overlap_sequences = ["".join(i) for i in non_overlap_sequences]
            # convert so its just hte synonymous site that differs
            motif_overlap_sequences = ["".join(["GC{0}".format(seq[i+2]) if "-" not in seq[i:i+3] else seq[i:i+3] for i in range(0, len(seq), 3)]) for seq in motif_overlap_sequences]
            motif_non_overlap_sequences = ["".join(["GC{0}".format(seq[i+2]) if "-" not in seq[i:i+3] else seq[i:i+3] for i in range(0, len(seq), 3)]) for seq in motif_non_overlap_sequences]
            # add to dictionary
            retained_overlap_sequences[id] = motif_overlap_sequences
            retained_non_overlap_sequences[id] = motif_non_overlap_sequences

            if codon_set:
                # now get the indices that create a stop from the overlap set
                stop_overlap_codon_indices = []
                non_stop_overlap_codon_indices = []
                for chunk in overlap_codon_indices:
                    codons = [alignment_set[0][chunk[1]:chunk[-1]+2], alignment_set[0][chunk[2]:chunk[-1]+3], alignment_set[1][chunk[1]:chunk[-1]+2], alignment_set[1][chunk[2]:chunk[-1]+3]]
                    if list(set(codon_set) & set(codons)):
                        stop_overlap_codon_indices.append(chunk)
                    else:
                        non_stop_overlap_codon_indices.append(chunk)

                stop_overlap_sequences = [[sequence[i[0]:i[-1]+1] if i != "-" else "---" for i in stop_overlap_codon_indices] for sequence in alignment_set]
                non_stop_overlap_sequences = [[sequence[i[0]:i[-1]+1] if i != "-" else "---" for i in non_stop_overlap_codon_indices] for sequence in alignment_set]
                motif_stop_overlap_sequences = ["".join(i) for i in stop_overlap_sequences]
                motif_non_stop_overlap_sequences = ["".join(i) for i in non_stop_overlap_sequences]
                motif_stop_overlap_sequences = ["".join(["GC{0}".format(seq[i+2]) if "-" not in seq[i:i+3] else seq[i:i+3] for i in range(0, len(seq), 3)]) for seq in motif_stop_overlap_sequences]
                motif_non_stop_overlap_sequences = ["".join(["GC{0}".format(seq[i+2]) if "-" not in seq[i:i+3] else seq[i:i+3] for i in range(0, len(seq), 3)]) for seq in motif_non_stop_overlap_sequences]
                retained_stop_overlap_sequences[id] = motif_stop_overlap_sequences
                retained_non_stop_overlap_sequences[id] = motif_non_stop_overlap_sequences

    if codon_set:
        return retained_overlap_sequences, retained_non_overlap_sequences, retained_stop_overlap_sequences, retained_non_stop_overlap_sequences
    else:
        return retained_overlap_sequences, retained_non_overlap_sequences

                # motif_seq = []
                # #for each codon in the sequences
                # # for seq_pos in range(0, len(seq), 3):
                # for seq_pos in range(0, 60, 3):
                #
                #     local_range = list(range(seq_pos, seq_pos + 3))
                #     next_range = list(range(seq_pos + 3, seq_pos + 6))
                #     codon = seq[seq_pos:seq_pos+3]
                #
                #     print(no, codon)
                #
                #     if len(list(set(local_range) & set(overlap_codon_indices))):
                #         # check if fourfold degenerate
                #         if only_fourfold:
                #             if codon in fourfold:
                #                 motif_seq.append(codon)
                #         else:
                #             motif_seq.append(codon)
                #         # if the next codon doesnt also appear, thats the end of the motif
                #         if not len(list(set(next_range) & set(overlap_codon_indices))):
                #             overlap_kept_sequences[no].append("".join(motif_seq))
                #             motif_seq = []
                #     elif len(list(set(local_range) & set(non_overlap_codon_indices))):
                #         # check if fourfold degenerate
                #         if only_fourfold:
                #             if codon in fourfold:
                #                 motif_seq.append(codon)
                #         else:
                #             motif_seq.append(codon)
                #         # if the next codon doesnt also appear, thats the end of the motif
                #         if not len(list(set(next_range) & set(non_overlap_codon_indices))):
                #             non_overlap_kept_sequences[no].append("".join(motif_seq))
                #             motif_seq = []



            # for no, sequence in enumerate(alignment_set):
            #     motif_seq = []
            #     # for each codon in the sequences
            #     for seq_pos in range(0, len(sequence), 3):
            #         local_range = list(range(seq_pos, seq_pos + 3))
            #         # if a position in the codon is in the overlap list
            #         if len(list(set(local_range) & set(overlap_indices))):
            #             codon_overlap = list(set(local_range) & set(overlap_indices))
            #             # ensure that the synonymous site is in the overlap, otherwise we don't want it
            #             if seq_pos + 2 in indices_to_keep:
            #                 codon = sequence[seq_pos:seq_pos+3]
            #                 # if we want to keep only the fourfold codons
            #                 if only_fourfold and codon in fourfold:
            #                     motif_seq.append(codon)
            #                 else:
            #                     motif_seq.append(codon)
            #
            #             next_range = list(range(seq_pos+3, seq_pos + 6))
            #             if not len(list(set(next_range) & set(overlap_indices))):
            #                 kept_sequences[no].append("".join(motif_seq))
            #                 motif_seq = []

            # print(kept_sequences)


        # kept_sequences = ["CCC".join(i) for i in kept_sequences]
        # print(kept_sequences[0])
        # remaining_motif_sequences[id] = kept_sequences
    # return []
    # return remaining_motif_sequences


def retain_only_fourfold_codons(query_sets):


    retained = []
    if len(query_sets) > 1:
        # set up the blank lists to hold fourfolds
        [retained.append([]) for i in range(len(query_sets))]
        # now loop over each codon
        for codon_index in range(len(query_sets[0])):
            query_codons = [query_sets[i][codon_index] for i in range(len(query_sets))]
            fourfolds = [i for i in query_codons if i in fourfold]
            if "---" in query_codons:
                if len(fourfolds):
                    [retained[i].append(query_sets[i][codon_index]) for i in range(len(query_sets))]
            else:
                query_aminos = list(set([codon_map[codon] for codon in query_codons]))
                if len(fourfolds) == len(query_codons) and len(query_aminos) == 1:
                    [retained[i].append(query_sets[i][codon_index]) for i in range(len(query_sets))]
    else:
        if len(query_sets):
            retained.append([i for i in query_sets[0] if i in fourfold])

    return retained


def convert_fourfold_to_gcn(codon_sets):
    converted = []
    for set in codon_sets:
        new_set = []
        for codon_set in set:
            new_set.append(["GC{0}".format(codon[-1]) if codon != "---" else "---" for codon in codon_set ])
        converted.append(new_set)

    return converted

def strings_to_alignment(strings):
    alignments = [[],[]]
    for set in strings:
        alignments[0].append(set[0])
        alignments[1].append(set[1])
    alignments = ["".join(i) for i in alignments]
    return alignments

def calc_strings_ds(strings):

    if sum([len(i) for i in strings]) != 0:
        alignments = strings_to_alignment(strings)
        ds = cons.calc_ds(alignments)
    else:
        ds = 0
    return ds

def ese_ds_wrapper(ids, motif_set_filepaths, alignments, output_directory, only_fourfold = True):
    """
    Wrapper to calculate the ds scores for motif sets in given alignemnts

    Args:
        ids (list): list of motif file ids
        motif_set_filepaths (dict): filepaths of motif sets
        alignments (dict): sequence alignments
        output_directory (str): path to output_directory
        only_fourfold (bool): if True, only retain fourfold degenerate codons
    """

    outputs = []
    if ids:
        for i, id in enumerate(ids):
            gen.print_parallel_status(i, ids)

            # get the filepath of the motifs
            filepath = motif_set_filepaths[id]
            # set up the temp file
            output_file = "{0}/{1}".format(output_directory, filepath.split("/")[-1])
            outputs.append(output_file)

            # only run if the set doesn't exist
            if output_file not in os.listdir(output_directory):
                # read in the motifs
                motif_set = read_motifs(filepath)
                # now get all the alignment parts that overlap the motif sets
                overlap_codons, non_overlap_codons, overlap_stops, overlap_non_stops = extract_motif_codons_from_alignments(alignments, motif_set, get_stops = True, only_both_hits = True)

                if only_fourfold:

                    overlap_codons = [retain_only_fourfold_codons(i) for i in overlap_codons]
                    non_overlap_codons = [retain_only_fourfold_codons(i) for i in non_overlap_codons]
                    stop_overlap_codons = [retain_only_fourfold_codons(i) for i in overlap_stops]
                    non_stop_overlap_codons = [retain_only_fourfold_codons(i) for i in overlap_non_stops]



                    # overlap_codons = convert_fourfold_to_gcn(overlap_codons)
                    # non_overlap_codons = convert_fourfold_to_gcn(non_overlap_codons)
                    # stop_overlap_codons = convert_fourfold_to_gcn(stop_overlap_codons)
                    # non_stop_overlap_codons = convert_fourfold_to_gcn(non_stop_overlap_codons)



                overlap_strings = []
                non_overlap_strings = []
                stop_overlap_strings = []
                non_stop_overlap_strings = []
                for set in overlap_codons:
                    overlap_strings.append(["".join(i) for i in set])
                for set in non_overlap_codons:
                    non_overlap_strings.append(["".join(i) for i in set])
                for set in stop_overlap_codons:
                    stop_overlap_strings.append(["".join(i) for i in set])
                for set in non_stop_overlap_codons:
                    non_stop_overlap_strings.append(["".join(i) for i in set])

                overlap_ds = calc_strings_ds(overlap_strings)
                non_overlap_ds = calc_strings_ds(non_overlap_strings)
                stop_overlap_ds = calc_strings_ds(stop_overlap_strings)
                non_stop_overlap_ds = calc_strings_ds(non_stop_overlap_strings)

                with open(output_file, "w") as outfile:
                    outfile.write("{0},{1},{2},{3},{4}\n".format(id, overlap_ds, non_overlap_ds, stop_overlap_ds, non_stop_overlap_ds))

    return outputs

def ese_ds_random_overlaps_wrapper(ids, motif_set_filepaths, alignments, output_directory, only_fourfold = True):
    """
    Wrapper to calculate the ds scores for motif sets in given alignemnts

    Args:
        ids (list): list of motif file ids
        motif_set_filepaths (dict): filepaths of motif sets
        alignments (dict): sequence alignments
        output_directory (str): path to output_directory
        only_fourfold (bool): if True, only retain fourfold degenerate codons
    """

    outputs = []
    if ids:
        for i, id in enumerate(ids):
            gen.print_parallel_status(i, ids)

            # get the filepath of the motifs
            filepath = motif_set_filepaths[id]
            # set up the temp file
            output_file = "{0}/{1}".format(output_directory, filepath.split("/")[-1])
            outputs.append(output_file)

            # only run if the set doesn't exist
            if output_file not in os.listdir(output_directory):
                # read in the motifs
                if id == "real":
                    motif_set = read_motifs(filepath)
                    # now get all the alignment parts that overlap the motif sets
                    overlap_codons, non_overlap_codons, overlap_stops, overlap_non_stops = extract_motif_codons_from_alignments(alignments, motif_set, get_stops = True)
                else:
                    overlap_ids, overlaps = gen.read_fasta(filepath)
                    overlap_list = {}
                    for id_count, overlap_id in enumerate(overlap_ids):
                        if overlaps[i] != "na":
                            overlap_list[overlap_id] = [int(pos) for pos in overlaps[i].split(",")]
                        else:
                            overlap_list[overlap_id] = []
                    overlap_codons, non_overlap_codons, overlap_stops, overlap_non_stops = extract_motif_codons_from_alignments(alignments, overlap_list, get_stops = True, query_overlaps = True)


                if only_fourfold:

                    overlap_codons = [retain_only_fourfold_codons(i) for i in overlap_codons]
                    non_overlap_codons = [retain_only_fourfold_codons(i) for i in non_overlap_codons]
                    stop_overlap_codons = [retain_only_fourfold_codons(i) for i in overlap_stops]
                    non_stop_overlap_codons = [retain_only_fourfold_codons(i) for i in overlap_non_stops]



                    # overlap_codons = convert_fourfold_to_gcn(overlap_codons)
                    # non_overlap_codons = convert_fourfold_to_gcn(non_overlap_codons)
                    # stop_overlap_codons = convert_fourfold_to_gcn(stop_overlap_codons)
                    # non_stop_overlap_codons = convert_fourfold_to_gcn(non_stop_overlap_codons)



                overlap_strings = []
                non_overlap_strings = []
                stop_overlap_strings = []
                non_stop_overlap_strings = []
                for set in overlap_codons:
                    overlap_strings.append(["".join(i) for i in set])
                for set in non_overlap_codons:
                    non_overlap_strings.append(["".join(i) for i in set])
                for set in stop_overlap_codons:
                    stop_overlap_strings.append(["".join(i) for i in set])
                for set in non_stop_overlap_codons:
                    non_stop_overlap_strings.append(["".join(i) for i in set])

                overlap_ds = calc_strings_ds(overlap_strings)
                non_overlap_ds = calc_strings_ds(non_overlap_strings)
                stop_overlap_ds = calc_strings_ds(stop_overlap_strings)
                non_stop_overlap_ds = calc_strings_ds(non_stop_overlap_strings)

                with open(output_file, "w") as outfile:
                    outfile.write("{0},{1},{2},{3},{4}\n".format(id, overlap_ds, non_overlap_ds, stop_overlap_ds, non_stop_overlap_ds))

    return outputs

def non_ese_ds_wrapper(ids, motif_set_filepaths, alignments, output_directory, only_fourfold = True):
    """
    Wrapper to calculate the ds scores for motif sets in given alignemnts

    Args:
        ids (list): list of motif file ids
        motif_set_filepaths (dict): filepaths of motif sets
        alignments (dict): sequence alignments
        output_directory (str): path to output_directory
        only_fourfold (bool): if True, only retain fourfold degenerate codons
    """

    outputs = []
    if ids:
        for i, id in enumerate(ids):
            gen.print_parallel_status(i, ids)

            # get the filepath of the motifs
            filepath = motif_set_filepaths[id]
            # set up the temp file
            output_file = "{0}/{1}".format(output_directory, filepath.split("/")[-1])
            outputs.append(output_file)

            # only run if the set doesn't exist
            if output_file not in os.listdir(output_directory):
                # read in the motifs
                motif_set = read_motifs(filepath)
                # now get all the alignment parts that overlap the motif sets
                non_overlap_stop_codons, non_overlap_non_stop_codons = extract_non_ese_stop_alignments(alignments, motif_set)


                if only_fourfold:

                    non_overlap_stop_codons = [retain_only_fourfold_codons(i) for i in non_overlap_stop_codons]
                    non_overlap_non_stop_codons = [retain_only_fourfold_codons(i) for i in non_overlap_non_stop_codons]

                    non_overlap_stop_codons = convert_fourfold_to_gcn(non_overlap_stop_codons)
                    non_overlap_non_stop_codons = convert_fourfold_to_gcn(non_overlap_non_stop_codons)


                non_overlap_stop_strings = []
                non_overlap_non_stop_strings = []

                for set in non_overlap_stop_codons:
                    non_overlap_stop_strings.append(["".join(i) for i in set])
                for set in non_overlap_non_stop_codons:
                    non_overlap_non_stop_strings.append(["".join(i) for i in set])

                overlap_ds = calc_strings_ds(non_overlap_stop_strings)
                non_overlap_ds = calc_strings_ds(non_overlap_non_stop_strings)

                with open(output_file, "w") as outfile:
                    outfile.write("{0},{1},{2}\n".format(id, overlap_ds, non_overlap_ds))

    return outputs


def ese_ds_mutation_wrapper(ids, motif_set_filepaths, alignments, output_directory, only_fourfold = True):
    """
    Wrapper to calculate the ds scores for motif sets in given alignemnts

    Args:
        ids (list): list of motif file ids
        motif_set_filepaths (dict): filepaths of motif sets
        alignments (dict): sequence alignments
        output_directory (str): path to output_directory
        only_fourfold (bool): if True, only retain fourfold degenerate codons
    """

    outputs = []
    if ids:
        for i, id in enumerate(ids):
            gen.print_parallel_status(i, ids)

            # get the filepath of the motifs
            filepath = motif_set_filepaths[id]
            # set up the temp file
            output_file = "{0}/{1}".format(output_directory, filepath.split("/")[-1])
            outputs.append(output_file)

            # only run if the set doesn't exist
            if output_file not in os.listdir(output_directory):
                # read in the motifs
                motif_set = read_motifs(filepath)
                motif_set = [i for i in motif_set if not len(re.findall("(?=(TAA|TAG|TGA))", i))]

                # now get all the alignment parts that overlap the motif sets
                overlap_one_away_codons, overlap_other_codons = extract_motif_codons_mutations_from_alignments(alignments, motif_set)

                if only_fourfold:

                    overlap_one_away_codons = [retain_only_fourfold_codons(i) for i in overlap_one_away_codons]
                    overlap_other_codons = [retain_only_fourfold_codons(i) for i in overlap_other_codons]

                    # overlap_one_away_codons = convert_fourfold_to_gcn(overlap_one_away_codons)
                    # overlap_other_codons = convert_fourfold_to_gcn(overlap_other_codons)

                one_away_strings = []
                other_strings = []

                for set in overlap_one_away_codons:
                    one_away_strings.append(["".join(i) for i in set])
                for set in overlap_other_codons:
                    other_strings.append(["".join(i) for i in set])

                one_away_ds = calc_strings_ds(one_away_strings)
                other_ds = calc_strings_ds(other_strings)

                with open(output_file, "w") as outfile:
                    outfile.write("{0},{1},{2}\n".format(id, one_away_ds, other_ds))

    return outputs


def calc_motif_sets_all_ds_wrapper(motif_set_list, motif_sets, codon_sets, sequence_alignments, output_file = None, output_directory = None):
    """
    Wrapper to calculate the ds score for the sequence parts where the synonymous
    site can and cannot form on of the codons in a codon set in sequence parts that overlap
    a set of given motifs

    Args:
        motif_set_list (list): list of integers to iterate over, corresponds to motif_sets
        motif_sets (dict): dictionary containing paths to motif set files
        codon_sets (list): a list of list containing codon sets to iterate over
        sequence_alignments (dict): dictionary containing the aligned sequences for each transcript
        output_file (str): if set, output to given file
        output_directory (str): if set, save files to output directory

    Returns:
        outputs (list): list of output file paths
    """

    if not output_file and not output_directory:
        print("Please provide output file or directory...")
        raise Exception

    outputs = []

    if motif_set_list:

        # print(mp.current_process().name.split("-")[-1], motif_set_list)
        for i, set_no in enumerate(motif_set_list):
            print("(W{0}) {1}/{2}".format(mp.current_process().name.split("-")[-1], i+1, len(motif_set_list)))
            motif_set = [i[0] for i in gen.read_many_fields(motif_sets[set_no], "\t") if "#" not in i[0] and ">" not in i[0]]
            stop_motif_set = [i for i in motif_set if len(re.findall("(TAA|TAG|TGA)", i))]
            non_stop_motif_set = [i for i in motif_set if i not in stop_motif_set]

            # now get all the parts of the sequences where only the motif sets occur
            retained_overlap_sequence, retained_non_overlap_sequence, retained_stop_overlap_sequences, reatined_non_stop_overlap_sequences = extract_motif_sequences_from_alignment(sequence_alignments, motif_set, codon_set = stops)

            # # now calculate the ds scores for each set
            retained_overlap_alignments = list_alignments_to_strings(retained_overlap_sequence)
            retained_non_overlap_alignments = list_alignments_to_strings(retained_non_overlap_sequence)
            retained_stop_overlap_alignments = list_alignments_to_strings(retained_stop_overlap_sequences)
            retained_non_stop_overlap_alignments = list_alignments_to_strings(reatined_non_stop_overlap_sequences)
            # stops_ese_alignments = list_alignments_to_strings(stops_extracted)
            # non_stops_ese_alignments = list_alignments_to_strings(non_stops_extracted)
            # #
            # # calculate the ds scores
            retained_overlap_ds = cons.calc_ds(retained_overlap_alignments)
            retained_non_overlap_ds = cons.calc_ds(retained_non_overlap_alignments)
            retained_stop_overlap_ds = cons.calc_ds(retained_stop_overlap_alignments)
            retained_non_stop_overlap_ds = cons.calc_ds(retained_non_stop_overlap_alignments)

            # set up the output filepath
            if output_file:
                output_file = output_file
            if output_directory:
                output_file = "{0}/{1}.{2}.txt".format(output_directory, set_no, random.random())
            outputs.append(output_file)

            with open(output_file, "w") as outfile:
                outfile.write("{0},{1},{2},{3},{4}\n".format(set_no, retained_overlap_ds, retained_non_overlap_ds, retained_stop_overlap_ds, retained_non_stop_overlap_ds))


    return outputs


def get_all_mutations(seq):
    mutation_seqs = []
    seq_nts = list(seq)
    for i,nt in enumerate(seq_nts):
        for mut_nt in nucleotides:
            if mut_nt != nt:
                query_nts = copy.deepcopy(seq_nts)
                query_nts[i] = mut_nt
                mut_seq = "".join(query_nts)
                mutation_seqs.append(mut_seq)
    return mutation_seqs


def calc_ds_mutation_wrapper(motif_set_list, motif_sets, codon_sets, sequence_alignments, output_file = None, output_directory = None):
    """
    Wrapper to calculate the ds score for the sequence parts where the synonymous
    site can and cannot form on of the codons in a codon set in sequence parts that overlap
    a set of given motifs

    Args:
        motif_set_list (list): list of integers to iterate over, corresponds to motif_sets
        motif_sets (dict): dictionary containing paths to motif set files
        codon_sets (list): a list of list containing codon sets to iterate over
        sequence_alignments (dict): dictionary containing the aligned sequences for each transcript
        output_file (str): if set, output to given file
        output_directory (str): if set, save files to output directory

    Returns:
        outputs (list): list of output file paths
    """

    if not output_file and not output_directory:
        print("Please provide output file or directory...")
        raise Exception

    outputs = []

    if motif_set_list:

        # print(mp.current_process().name.split("-")[-1], motif_set_list)
        for i, set_no in enumerate(motif_set_list):
            print("(W{0}) {1}/{2}".format(mp.current_process().name.split("-")[-1], i+1, len(motif_set_list)))
            motif_set = [i[0] for i in gen.read_many_fields(motif_sets[set_no], "\t") if "#" not in i[0] and ">" not in i[0]]
            stop_motif_set = [i for i in motif_set if len(re.findall("(TAA|TAG|TGA)", i))]
            non_stop_motif_set = [i for i in motif_set if i not in stop_motif_set]

            non_stop_mutation_set = [get_all_mutations(i) for i in non_stop_motif_set]
            mutation_set_no_stops = [non_stop_motif_set[i] for i, motif_set in enumerate(non_stop_mutation_set) if sum([len(re.findall("(TAA|TAG|TGA)", motif)) for motif in motif_set]) == 0]
            mutation_set_stops = [i for i in non_stop_motif_set if i not in mutation_set_no_stops]

            retained_mutation_no_stops =  extract_motif_sequences_from_alignment(sequence_alignments, mutation_set_no_stops)[0]
            retained_mutation_stops =  extract_motif_sequences_from_alignment(sequence_alignments, mutation_set_stops)[0]

            # # now calculate the ds scores for each set
            stops_alignments = list_alignments_to_strings(retained_mutation_stops)
            non_stops_alignments = list_alignments_to_strings(retained_mutation_no_stops)

            # # calculate the ds scores)
            if sum([len(i) for i in stops_alignments]):
                stops_ds = cons.calc_ds(stops_alignments)
            else:
                stops_ds = "nan"
            if sum([len(i) for i in non_stops_alignments]):
                non_stops_ds = cons.calc_ds(non_stops_alignments)
            else:
                non_stops_ds = "nan"

            # set up the output filepath
            if output_file:
                output_file = output_file
            if output_directory:
                output_file = "{0}/{1}.{2}.txt".format(output_directory, set_no, random.random())
            outputs.append(output_file)

            with open(output_file, "w") as outfile:
                outfile.write("{0},{1},{2}\n".format(set_no, stops_ds, non_stops_ds))

    return outputs


def calc_motif_sets_codons_ds_wrapper(motif_set_list, motif_sets, codon_sets, sequence_alignments, output_file = None, output_directory = None):
    """
    Wrapper to calculate the ds score for the sequence parts where the synonymous
    site can and cannot form on of the codons in a codon set in sequence parts that overlap
    a set of given motifs

    Args:
        motif_set_list (list): list of integers to iterate over, corresponds to motif_sets
        motif_sets (dict): dictionary containing paths to motif set files
        codon_sets (list): a list of list containing codon sets to iterate over
        sequence_alignments (dict): dictionary containing the aligned sequences for each transcript
        output_file (str): if set, output to given file
        output_directory (str): if set, save files to output directory

    Returns:
        outputs (list): list of output file paths
    """

    if not output_file and not output_directory:
        print("Please provide output file or directory...")
        raise Exception

    outputs = []

    # print(mp.current_process().name.split("-")[-1], motif_set_list)
    for i, set_no in enumerate(motif_set_list):
        print("(W{0}) {1}/{2}".format(mp.current_process().name.split("-")[-1], i+1, len(motif_set_list)))

        motif_set = [i[0] for i in gen.read_many_fields(motif_sets[set_no], "\t") if "#" not in i[0] and ">" not in i[0]]

        # now get all the parts of the sequences where only the motif sets occur
        restricted_sequences = extract_motif_sequences_from_alignment(sequence_alignments, motif_set)
        # set up the output filepath
        if output_file:
            output_file = output_file
        if output_directory:
            output_file = "{0}/{1}.{2}.txt".format(output_directory, set_no, random.random())
        outputs.append(output_file)
        # calculcate the ds score for the codon set
        calc_motif_set_codons_ds(codon_sets, restricted_sequences, output_file)

    return outputs


def calc_purine_content(seqs, reverse = None):
    """
    Calculate the purine (or pyrimidine) content of a list of sequences

    Args:
        seqs (list): list of sequences
        reverse (bool): if set, calculate the pyrimidine content instead
    Returns:
        purine_content (float): purine content of sequences
    """

    query_set = ["A", "G"]
    if reverse:
        query_set = ["C", "T"]

    seq_nts = list("".join(seqs))
    purine_content = np.divide(len([i for i in seq_nts if i in query_set]), len(seq_nts))
    return purine_content


def clean_feature_file(bed_file):
    """
    After extracting features, may want to make it in a more usable format

    Args:
        bed_file (str): path to bed file
    """

    temp_file = "{0}.bed".format(random.random())
    entries = gen.read_many_fields(bed_file, "\t")
    with open(temp_file, "w") as outfile:
        for entry in entries:
            entry[3] = "{0}.{1}".format(entry[3], entry[4])
            entry[4] = "."
            outfile.write("{0}\n".format("\t".join(entry)))

    gen.run_process(["mv", temp_file, bed_file])
    gen.remove_file(temp_file)



def extract_gtf_feature(input_file, required_feature, ids_to_keep = None, filter_by_gene = None):
    """
    Get a feature from a gtf file

    Args:
        input_file (str): path to file containing features
        required_feature: (str): feature to keep
        ids_to_keep (list): if set, list containing ids to keep
        filter_by_gene (bool): if true, match gene ids rather than transcript ids

    Returns:
        features_list (dict): dict containing features, sorted according to
            their position in the transcript. dict[transcript_id] = [cds_features]
    """

    print("Getting {0}s...".format(required_feature))

    features = gen.read_many_fields(input_file, "\t")
    features_list = collections.defaultdict(lambda: [])
    # get a list of all that match
    for feature in features:
        if feature[-1] == required_feature:
            # if filtering by gene, check that gene is in the input list, else
            # check the transcript
            if ids_to_keep:
                if filter_by_gene and feature[6] in ids_to_keep:
                    features_list[feature[3]].append(feature)
                elif feature[3] in ids_to_keep:
                    features_list[feature[3]].append(feature)
            # otherwise keep all features
            else:
                features_list[feature[3]].append(feature)

    for id in features_list:
        # now we need to sort the exons in order, reversing if on the minus strand
        strand = features_list[id][0][5]
        if strand == "+":
            features_list[id] = sorted(features_list[id], key = lambda x:x[1])
        elif strand == "-":
            features_list[id] = sorted(features_list[id], key = lambda x:x[2], reverse = True)

    return features_list


def extract_aligment_sequence_parts(cds_pairs, reverse = False):

    # put the muscle executable in tools directory for your os
    muscle_exe = "../tools/muscle3.8.31_i86{0}64".format(sys.platform)
    if not os.path.isfile(muscle_exe):
        print("Could not find the MUSCLE exe {0}...".format(muscle_exe))
        raise Exception

    # setup the muscle alignment
    alignment_functions = cons.Alignment_Functions(muscle_exe)

    kept_sequence1 = []
    kept_sequence2 = []

    for focal_seq in cds_pairs:
        ortholog_seq = cds_pairs[focal_seq]

        # get the alignments for the sequences
        focal_iupac_protein = Seq(focal_seq, IUPAC.unambiguous_dna).translate()
        ortholog_iupac_protein = Seq(ortholog_seq, IUPAC.unambiguous_dna).translate()
        alignment_functions.align_seqs(focal_iupac_protein, ortholog_iupac_protein)
        # extract the alignments
        alignment_functions.extract_alignments()
        # now we want to get the nucleotide sequences for the alignments
        aligned_sequences = alignment_functions.revert_alignment_to_nucleotides(input_seqs = [focal_seq, ortholog_seq])
        # clean up the files
        alignment_functions.cleanup()

        # extract the parts of the sequences that are distance from stop codon
        # get_sequence_parts_near_stop(aligned_sequences)


        keep_only_one_synonymous_from_stop = get_alignment_one_synonymous_from_stop(aligned_sequences, reverse = reverse)
        kept_sequence1.append(keep_only_one_synonymous_from_stop[0])
        kept_sequence2.append(keep_only_one_synonymous_from_stop[1])

    return "".join(kept_sequence1), "".join(kept_sequence2)


def extract_alignments(seq1, seq2, muscle_exe = None):
    """
    Given 2 sequences, align protein sequences and convert back to nucleotide sequences

    Args:
        seq1 (str): sequence 1
        seq2 (str): sequence 2
        muscle_ext (str): path to muscle exe file

    Returns:
        aligned_sequences (list): list containing the aligned sequences
    """

    if not muscle_exe:
        print("Could not find the MUSCLE exe {0}...".format(muscle_exe))
        raise Exception

    # setup the muscle alignment
    alignment_functions = cons.Alignment_Functions(muscle_exe)
    # get the alignments for the sequences
    seq1_iupac_protein = Seq(seq1, IUPAC.unambiguous_dna).translate()
    seq2_iupac_protein = Seq(seq2, IUPAC.unambiguous_dna).translate()
    alignment_functions.align_seqs(seq1_iupac_protein, seq2_iupac_protein)
    # extract the alignments
    alignment_functions.extract_alignments()
    # now we want to get the nucleotide sequences for the alignments
    aligned_sequences = alignment_functions.revert_alignment_to_nucleotides(input_seqs = [seq1, seq2])
    # clean up the files
    alignment_functions.cleanup()

    return aligned_sequences


def extract_alignments_from_file(cds_fasta, ortholog_fasta, ortholog_transcript_links, output_file):
    """
    Given a file of CDSs, extract the alignment as given by a link file with
    the orthologous sequence

    Args:
        cds_fasta (str): path to cds fasta file
        ortholog_fasta (str): path to ortholog fasta file
        ortholog_transcript_links (str): path to file defining the ortholog sequence to use (lowest dS)
        output_file (str): path to output file
    """

    print("Extracting sequence alignments to file...")

    # get a list of the links
    links = {i[0]: i[1] for i in gen.read_many_fields(ortholog_transcript_links, "\t")}
    # get a list of cds
    cds_names, cds_seqs = gen.read_fasta(cds_fasta)
    cds_list = {name: cds_seqs[i] for i, name in enumerate(cds_names)}
    # get a list of the orthologs
    ortholog_names, ortholog_seqs = gen.read_fasta(ortholog_fasta)
    ortholog_cds_list = {links[name]: ortholog_seqs[ortholog_names.index(links[name])] for name in cds_list if links[name] in ortholog_names}
    # pair the sequences
    cds_pairs = {id: [cds_list[id], ortholog_cds_list[links[id]]] for id in cds_list}

    # put the muscle executable in tools directory for your os
    muscle_exe = "../tools/muscle3.8.31_i86{0}64".format(sys.platform)

    with open(output_file, "w") as outfile:
        for id in pbar(cds_pairs):
            alignments = extract_alignments(cds_pairs[id][0], cds_pairs[id][1], muscle_exe = muscle_exe)
            outfile.write(">{0}\n{1},{2}\n".format(id, alignments[0], alignments[1]))


def extract_gtf_features(input_list, gtf_file_path, transcript_id_search_pattern, gene_id_search_pattern, filter_by_transcript = None, filter_by_gene = None):
    """
    Given a .gtf file, filter the entries based on an input list

    Args:
        input_list (list): list containing ids for which to filter from
        gtf_file_path (str): path to gtf file
        filter_by_transcript (bool): if true, filter by transcript id
        filter_by_gene (bool): if true, filter by gene id

    Returns:
        outputs (list): list containing features that passed the feature filtering
    """

    print("Getting features from .gtf file...")

    outputs = []

    # only run if necessary
    if len(input_list):
        # read the gtf file and remove meta data
        entries = [i for i in gen.read_many_fields(gtf_file_path, "\t") if "#" not in i[0]]

        for entry in entries:
            entry_info = entry[8]
            # look for transcript id
            transcript_search_results = re.search(transcript_id_search_pattern, entry_info)
            gene_search_results = re.search(gene_id_search_pattern, entry_info)
            try:
                # try to get transcript id and gene id
                current_transcript_id = transcript_search_results.group(0)
                current_gene_id = gene_search_results.group(0)
                # if filtering by ids, do that here
                passed_filter = True
                if filter_by_transcript:
                    passed_filter = current_transcript_id in input_list
                if filter_by_gene:
                    passed_filter = current_gene_id in input_list
                # if it passes the filter
                if passed_filter:
                    exon_number_search_results = re.search(exon_number_pattern, entry_info)
                    try:
                        # try to get the exon number
                        exon_number = exon_number_search_results.group(1)
                    except AttributeError:
                        exon_number = "nan"
                    # minus 1 for the start because gtf are in base 0, but want base 1
                    outputs.append([entry[0], str(int(entry[3])-1), entry[4], current_transcript_id, exon_number, entry[6], current_gene_id, entry[2]])
            # failed to find transcript or gene id
            except AttributeError:
                pass

    return outputs


def extract_cds_features(input_file, input_list, filter_by_gene = None):
    """
    Get all the CDS features from a file

    Args:
        input_file (str): path to file containing features
        input_list (list): list containing ids to keep
        filter_by_gene (bool): if true, match gene ids rather than transcript ids

    Returns:
        cds_features_list (dict): dict containing cds features, sorted according to
            their position in the transcript. dict[transcript_id] = [cds_features]
    """

    print("Getting cds features...")


    features = gen.read_many_fields(input_file, "\t")
    # get the features labelled as stop codons
    stop_codons = extract_stop_codon_features(features, input_list, filter_by_gene = filter_by_gene)
    start_codons = extract_start_codon_features(features, input_list, filter_by_gene = filter_by_gene)
    cds_features_list = collections.defaultdict(lambda: [])
    # get a list of all the exon parts that contribute to the cds
    # [cds_features_list[feature[3]].append(feature) for feature in self.features if feature[-1] == "CDS" and feature[3] in self.ids]
    for feature in features:
        if feature[-1] == "CDS":
            # if filtering by gene, check that gene is in the input list, else
            # check the transcript
            if filter_by_gene and feature[6] in input_list:
                cds_features_list[feature[3]].append(feature)
            elif feature[3] in input_list:
                cds_features_list[feature[3]].append(feature)

    for id in cds_features_list:
        # get the start codon coordinates if they exist
        start_codon = start_codons[id]
        start_codon_coordinates = [[i[1], i[2]] for i in start_codon]
        # get the stop codon coordinates if they exist
        stop_codon = stop_codons[id]
        stop_codon_coordinates = [[i[1], i[2]] for i in stop_codon]
        # get the cds coordinates if they arent in the stop codon list
        cds_features_list[id] = [i for i in cds_features_list[id] if [i[1], i[2]] not in stop_codon_coordinates and [i[1], i[2]] not in start_codon_coordinates]
        # now we need to sort the exons in order, reversing if on the minus strand
        strand = cds_features_list[id][0][5]
        if strand == "+":
            cds_features_list[id] = sorted(cds_features_list[id], key = lambda x:x[1])
        elif strand == "-":
            cds_features_list[id] = sorted(cds_features_list[id], key = lambda x:x[2], reverse = True)
        # add the start codon if it exists
        for start_codon in start_codons[id]:
            cds_features_list[id] = [start_codon] + cds_features_list[id]
        # add the stop codon if it exists
        for stop_codon in stop_codons[id]:
            cds_features_list[id].append(stop_codon)

    return cds_features_list



def extract_multi_exons_entries_to_bed(input_bed, output_bed = None):

    entries = gen.read_many_fields(input_bed, "\t")
    if not output_bed:
        output_bed_file = "temp_files/{0}.bed".format(random.random())
    else:
        output_bed_file = output_bed

    with open(output_bed_file, "w") as outfile:
        for entry in entries:
            sizes = [int(i) for i in entry[10].split(",") if i]
            starts = [int(i) for i in entry[11].split(",") if i]
            for i, exon in enumerate(sizes):
                start = int(entry[1]) + starts[i]
                end = start + sizes[i]
                bed_entry = copy.deepcopy(entry)
                bed_entry[1] = start
                bed_entry[2] = end
                bed_entry[3] = "{0}.{1}".format(bed_entry[3], i+1)
                outfile.write("{0}\n".format("\t".join(gen.stringify(bed_entry[:6]))))

    if not output_bed:
        gen.run_process(["mv", output_bed_file, input_bed])
        gen.remove_file(output_bed_file)


def keep_fourfold_indices(indices, sequences):

    kept_indices = []

    length = len(sequences[0])
    # for each of the codons in the sequences
    for i in range(0, length, 3):
        # get the indices of the codon
        local_range = list(range(i, i+3))
        # if any of those indices are in ese hits
        if list(set(local_range) & set(indices)):
            codons = [i[local_range[0]:local_range[-1]+1] for i in sequences]
            local_fourfold = [i for i in codons if i in fourfold]
            if len(local_fourfold):
                kept_indices.extend([i for i in local_range if i in indices])

    return kept_indices





def extract_stop_codon_features(input_features, input_list, filter_by_gene = None):
    """
    Get all the features that match stop codon

    Args:
        input_features (list): list containing features from gtf file
        input_list (list): list containing ids to keep

    Returns:
        feature_list (dict): dictionary containing features corresponding to
            stop codons. dict[transcript_id] = feature
    """

    print("Getting stop_codon features...")

    feature_list = collections.defaultdict(lambda: [])
    for feature in input_features:
        if feature[-1] == "stop_codon":
            # if filtering by gene, check that gene is in the input list, else
            # check the transcript
            if filter_by_gene and feature[6] in input_list:
                feature_list[feature[3]].append(feature)
            elif feature[3] in input_list:
                feature_list[feature[3]].append(feature)
    return feature_list


def extract_start_codon_features(input_features, input_list, filter_by_gene = None):
    """
    Get all the features that match start codon

    Args:
        input_features (list): list containing features from gtf file
        input_list (list): list containing ids to keep

    Returns:
        feature_list (dict): dictionary containing features corresponding to
            start codons. dict[transcript_id] = feature
    """

    print("Getting start_codon features...")

    feature_list = collections.defaultdict(lambda: [])
    for feature in input_features:
        if feature[-1] == "start_codon":
            # if filtering by gene, check that gene is in the input list, else
            # check the transcript
            if filter_by_gene and feature[6] in input_list:
                feature_list[feature[3]].append(feature)
            elif feature[3] in input_list:
                feature_list[feature[3]].append(feature)
    return feature_list


def extract_transcript_features(features_list, id_list):
    """
    Get the list of transcript features

    Args:
        features_list (list): list of features from the gtf file that were extracted
        id_list (list): list of ids to search for

    Returns:
        feature_list (dict): dictionary containing transcripts.
            dict[transcript_id] = feature
    """

    print("Getting transcript features...")

    feature_list = collections.defaultdict(lambda: [])
    for feature in features_list:
        if feature[-1] == "transcript" and feature[3] in id_list:
            feature_list[feature[3]].append(feature)
    return feature_list


def filter_bed_file(input_bed, filter_columns, filter_values, inclusive, output_file = None):
    """
    Given a bed file, filter columns based on values

    Args:
        input_bed (str): path to input bed file
        filter_columns (list): list containing columns to filter by
        filter_values (list): list of lists containing values that the column must correspond to
        inclusive (list): list of bool values saying whether the filter value should be included or not
        output_file (str): if set, use as output path. If not, assume bed file is to be filtered in place
    """

    # create temp file if no output file is given
    if not output_file:
        file_to_write = "temp_files/{0}.bed".format(random.random())
    else:
        file_to_write = output_file

    # read entries
    entries = gen.read_many_fields(input_bed, "\t")
    # for each of the columns, filter by values
    for i, column in enumerate(filter_columns):
        values = filter_values[i]
        inc = inclusive[i]
        if inc:
            entries = [entry for entry in entries if entry[column] in values]
        else:
            entries = [entry for entry in entries if entry[column] not in values]

    # write the remaining entries to file
    with open(file_to_write, "w") as outfile:
        [outfile.write("{0}\n".format("\t".join(entry))) for entry in entries]

    # remove the temp file if created
    if not output_file:
        gen.run_process(["mv", file_to_write, input_bed])
        gen.remove_file(file_to_write)


def filter_by_exon_number(input_bed, input_fasta, single_exon_bed, multi_exon_bed, single_exon_fasta, multi_exon_fasta):
    """
    Given a bed file of CDS entries and a fasta file containing the sequences, filter to group the sequences
    into single and multi exon cases.

    Args:
        input_bed (str): path to bed file
        input_fasta (str): path to fasta file
        single_exon_bed (str): path to file containing bed entries for single exon genes
        multi_exon_bed (str): path to file containing bed entries for multi exon genes
        single_exon_fasta (str): path to file containing single exon gene sequences
        multi_exon_fasta (str): path to file containing multi exon gene sequences
    """

    # get a list of bed entries and hold by transcript id
    bed_entries = gen.read_many_fields(input_bed, "\t")
    info = collections.defaultdict(lambda: [])
    [info[i[3].split(".")[0]].append(i) for i in bed_entries if i[-1] == "CDS"]
    # now get a list of sequences
    names, seqs = gen.read_fasta(input_fasta)
    seq_list = {name: seqs[i] for i, name in enumerate(names)}

    with open(single_exon_bed, "w") as outfile1:
        with open(multi_exon_bed, "w") as outfile2:
            with open(single_exon_fasta, "w") as outfile3:
                with open(multi_exon_fasta, "w") as outfile4:
                    for id in seq_list:
                        if len(info[id]) > 1:
                            [outfile2.write("{0}\n".format("\t".join(i))) for i in info[id]]
                            outfile4.write(">{0}\n{1}\n".format(id, seq_list[id]))
                        else:
                            [outfile1.write("{0}\n".format("\t".join(i))) for i in info[id]]
                            outfile3.write(">{0}\n{1}\n".format(id, seq_list[id]))


def filter_one_transcript_per_gene(transcript_list):
    """
    If there is more than one transcript per gene, keep only the longest.

    Args:
        transcript_list (dict): dict of transcripts

    Returns:
        longest_transcripts (dict): dict containing the longest transcript per gene
    """

    print("Filtering to one transcript per gene...")

    ids = []
    # get a list of genes with the associated transcripts
    gene_list = collections.defaultdict(lambda: [])
    for transcript in transcript_list:
        transcript_info = transcript_list[transcript][0]
        gene_list[transcript_info[6]].append(transcript_info)
        ids.append(transcript_info[6])

    # now get the longest one
    longest_transcripts = {}
    for gene_id in gene_list:
        if len(gene_list[gene_id]) > 1:
            lengths = [int(i[2]) - int(i[1]) for i in gene_list[gene_id]]
            longest_transcripts[gene_list[gene_id][lengths.index(max(lengths))][3]] = gene_list[gene_id][lengths.index(max(lengths))]
        else:
            longest_transcripts[gene_list[gene_id][0][3]] = gene_list[gene_id][0]
    return longest_transcripts


def get_alignment_one_synonymous_from_stop(aligned_sequences):

    seq1, seq2 = aligned_sequences[0], aligned_sequences[1]
    kept_seq1, kept_seq2 = [], []
    removed_seq1, removed_seq2 = [],[]

    # for each codon in the alignment
    for i in range(0, len(seq1)-3, 3):
        one_away_stop_list = []
        # get the codons in the +1 reading frame
        query_codons_f1 = [seq1[i+1:i+4], seq2[i+1:i+4]]
        query_codons_f1 = [i for i in query_codons_f1 if "-" not in i]
        # get the codons in the +2 reading frame
        query_codons_f2 = [seq1[i+2:i+5], seq2[i+2:i+5]]
        query_codons_f2 = [i for i in query_codons_f2 if "-" not in i]
        # now get all the codons that are 1 away from the query codon in both frames, and keep if its a stop
        one_away_stop_list.extend([[j for j in i if j in stops] for i in [one_away_codons[i][1] for i in query_codons_f1]])
        one_away_stop_list.extend([[j for j in i if j in stops] for i in [one_away_codons[i][0] for i in query_codons_f2]])
        # get the number of stops for the set of query codons
        one_away_stop_count = sum([len(i) for i in one_away_stop_list])
        # add codons to list if one stop exists in either frame
        if one_away_stop_count > 0:
            kept_seq1.append(seq1[i:i+3])
            kept_seq2.append(seq2[i:i+3])
        else:
            removed_seq1.append(seq1[i:i+3])
            removed_seq2.append(seq2[i:i+3])

    return [["".join(kept_seq1), "".join(kept_seq2)], ["".join(removed_seq1), "".join(removed_seq2)]]

# def get_alignment_synonymous_stop(aligned_sequences):
#
#     seq1, seq2 = aligned_sequences[0], aligned_sequences[1]
#     kept_seq1, kept_seq2 = [], []
#     removed_seq1, removed_seq2 = [],[]
#
#     # for each codon in the alignment
#     for i in range(0, len(seq1)-3, 3):
#         one_away_stop_list = []
#         # get the codons in the +1 reading frame
#         query_codons_f1 = [seq1[i+1:i+4], seq2[i+1:i+4]]
#         query_codons_f1 = [i for i in query_codons_f1 if "-" not in i]
#         # get the codons in the +2 reading frame
#         query_codons_f2 = [seq1[i+2:i+5], seq2[i+2:i+5]]
#         query_codons_f2 = [i for i in query_codons_f2 if "-" not in i]
#         # now get all the codons that are 1 away from the query codon in both frames, and keep if its a stop
#         one_away_stop_list.extend([[j for j in i if j in stops] for i in [one_away_codons[i][1] for i in query_codons_f1]])
#         one_away_stop_list.extend([[j for j in i if j in stops] for i in [one_away_codons[i][0] for i in query_codons_f2]])
#         # get the number of stops for the set of query codons
#         one_away_stop_count = sum([len(i) for i in one_away_stop_list])
#         # add codons to list if one stop exists in either frame
#         if one_away_stop_count > 0:
#             kept_seq1.append(seq1[i:i+3])
#             kept_seq2.append(seq2[i:i+3])
#         else:
#             removed_seq1.append(seq1[i:i+3])
#             removed_seq2.append(seq2[i:i+3])
#
#     return [["".join(kept_seq1), "".join(kept_seq2)], ["".join(removed_seq1), "".join(removed_seq2)]]


def generate_genome_features_dataset(dataset_name, dataset_gtf_file, dataset_features_output_bed, transcript_id_search_pattern, gene_id_search_pattern, input_list = None, filter_by_transcript = None, filter_by_gene = None):
    """
    Create a dataset containing all the relevant genome features from the
    given list.

    Args:
        dataset_name (str): name for the dataset
        dataset_gtf_file (str): path to genome gtf file
        dataset_features_output_bed (str): path to bed file that contains the extacted features
        input_list (list): if set, use as a list to filter transcripts from
        filter_by_transcript (bool): if true, filter by transcript id
        filter_by_gene (bool): if true, filter by gene id
    """

    print("Generating genome dataset...")

    # make sure that if we want filtering, we have the ids
    if filter_by_transcript or filter_by_gene:
        if not input_list:
            print("Input list required...")
            raise Exception
    if not input_list:
        input_list = ["foo"]

    # do the filtering
    args = [dataset_gtf_file, transcript_id_search_pattern, gene_id_search_pattern, filter_by_transcript, filter_by_gene]
    # reduce the number of workers here otherwise it doesnt like it
    features = gen.run_parallel_function(input_list, args, extract_gtf_features, parallel = True, workers = int(os.cpu_count() / 2) - 1)
    # output the features to the bed file
    with open(dataset_features_output_bed, "w") as outfile:
        outfile.write("# {0} lines in {1}\n".format(dataset_name, dataset_gtf_file.split("/")[-1]))
        [outfile.write("{0}\n".format("\t".join(feature))) for feature in features]

    print("{0} dataset created...".format(dataset_name))


def get_clean_genes(input_fasta, input_bed, output_bed):
    """
    Get a list of gene names for sequences in fasta

    Args:
        input_fasta (str): path to fasta file
        input_bed (str): path to bed file
        output_bed (str): path to output file
    """

    ids = gen.read_fasta(input_fasta)[0]
    entries = gen.read_many_fields(input_bed, "\t")
    clean_genes = {i[3].split(".")[0]: i[6] for i in entries if i[3].split(".")[0] in ids and i[3].split(".")[0]}
    with open(output_bed, "w") as outfile:
        [outfile.write("{0}\t{1}\n".format(i, clean_genes[i])) for i in clean_genes]


def get_coding_exon_coordinates(full_bed, cds_bed, output_file):
    """
    Given a list of exons that make up the cds, and a list of all exons,
    filter to only include fully coding exonsself.

    Args:
        full_bed (str): path to file containing all exons
        cds_bed (str): path to file containing cds exons
        output_file (str): path to output_file
    """

    print("Getting coding exon coordinates...")

    gen.create_output_directories("temp_files")
    temp_file1 = "temp_files/{0}.bed".format(random.random())
    temp_file2 = "temp_files/{0}.bed".format(random.random())
    # intersect the bed file to get all 100% hits
    intersect_coding_exons(full_bed, cds_bed, temp_file1)
    # now remove any potential lingering terminal exons
    remove_terminal_exons(full_bed, temp_file1, temp_file2)
    # sort and clean up output
    sort_bed_file(temp_file2, output_file)
    gen.remove_file(temp_file1)
    gen.remove_file(temp_file2)


def get_intron_coordinates(input_bed, output_bed):
    """
    Given a bed file of exon coordinates, extract the intron coordinates

    Args:
        input_bed (str): path to bed file containing exon coordinates
        output_bed (str): path to output file
    """

    exons = gen.read_many_fields(input_bed, "\t")
    exon_list = collections.defaultdict(lambda: collections.defaultdict())
    for i, exon in enumerate(exons):
        id = exon[3].split(".")[0]
        exon_no = int(exon[3].split(".")[1])
        exon_list[id][exon_no] = exon

    with open(output_bed, "w") as outfile:
        for id in exon_list:
            for exon_no in exon_list[id]:
                # check if there is another exon afterwards
                if exon_no + 1 in exon_list[id]:
                    entry = copy.deepcopy(exon_list[id][exon_no])
                    strand = exon_list[id][exon_no][5]
                    if strand == "-":
                        temp_start = int(exon_list[id][exon_no+1][2])
                        temp_end = int(exon_list[id][exon_no][1])
                        if temp_start > temp_end:
                            start = temp_end
                            end = temp_start
                        else:
                            start = temp_start
                            end = temp_end
                        entry[1] = start
                        entry[2] = end
                        entry[3] = "{0}.{1}-{2}".format(id, exon_no, exon_no+1)
                    else:
                        entry[1] = exon_list[id][exon_no][2]
                        entry[2] = exon_list[id][exon_no+1][1]
                        entry[3] = "{0}.{1}-{2}".format(id, exon_no, exon_no+1)
                    outfile.write("{0}\n".format("\t".join(gen.stringify(entry))))

def get_motifs_overlaps(seqs, motif_set):

    overlap_indices = []
    for i, motif in enumerate(motif_set):
        print("{0}/{1}".format(i+1, len(motif_set)))
        for seq in seqs:
            hits = re.finditer('(?=({0}))'.format(motif), seq)
            [overlap_indices.extend(list(range(i.start(), i.start() + len(motif)))) for i in hits]
    # overlap_indices = sorted(list(set(overlap_indices)))

    return overlap_indices



def get_motifs_overlap_indices(seqs, motif_set, reverse = None):
    """
    For a set of aligned sequences, get a list of all indices that overlap something
    in the motif set

    Args:
        seqs (list): list of seqeunces
        motif_set (list): list of motifs to query

    Returns:
        overlap_indices (list): list of indices that overlap
    """

    overlap_indices = []
    for motif in motif_set:
        for seq in seqs:
            hits = re.finditer('(?=({0}))'.format(motif), seq)
            [overlap_indices.extend(list(range(i.start(), i.start() + len(motif)))) for i in hits]
    overlap_indices = sorted(list(set(overlap_indices)))

    # if reverse, get all sites that arent in the overlaps
    if reverse:
        non_overlaps = [i for i in range(0, len(seqs[0])) if i not in overlap_indices]
        overlap_indices = non_overlaps

    return overlap_indices


def get_non_coding_exon_coordinates(full_bed, reduced_bed, output_file):
    """
    Get all non coding exons

    Args:
        full_bed (str): path to bed file containing full exon entries
        reduced_bed (str): path to bed file containing cds exon entries
        output_file (str): path to output file
    """

    print("Getting non coding exon coordinates...")

    gen.create_output_directories("temp_files")
    temp_file = "temp_files/{0}.bed".format(random.random())
    # intersect the files, keeping only those entries in the full bed that
    # arent in the reduced bed
    intersect_non_coding_exons(full_bed, reduced_bed, temp_file)
    sort_bed_file(temp_file, output_file)


def get_ortholog_transcript_pairings(input_file1, input_file2, ortholog_pairs_file, ortholog_fasta, output_file):
    """
    Give two files containing transcript-gene pairings, link the transcripts in the first
    file to the the transcripts in the second.

    Args:
        input_file1 (str): path to the first transcript-gene pairings file
        input_file2 (str): path to the second transcript-gene pairings file
        ortholog_pairs_file (str): path to the file containing the ensembl biomart ortholog pairings
        ortholog_fasta (str): path to fasta file containing filtered orthlog seqs
        output_file (str): path to output file
    """

    # get the names of the sequnces in the fasta
    fasta_names = gen.read_fasta(ortholog_fasta)[0]
    # read in the transcript-gene links
    transcripts1 = gen.read_many_fields(input_file1, "\t")
    transcripts2 = gen.read_many_fields(input_file2, "\t")
    # get the list of genes that are queried
    queried_genes = [i[1] for i in transcripts1]
    # get the links between transcripts and genes
    # for the first file get transcript-genes
    transcript_gene_links1 = collections.defaultdict(lambda: [])
    [transcript_gene_links1[i[0]].append(i[1]) for i in transcripts1]
    # for the second file get gene-transcripts
    transcript_gene_links2 = collections.defaultdict(lambda: [])
    [transcript_gene_links2[i[1]].append(i[0]) for i in transcripts2]
    # read in the ortholog pairings
    ortholog_pairings = collections.defaultdict(lambda: [])
    [ortholog_pairings[i[0]].append(i[1]) for i in gen.read_many_fields(ortholog_pairs_file, "\t") if i[0] in queried_genes]
    # now links them together
    transcript_links = collections.defaultdict(lambda: [])
    for transcript_id in transcript_gene_links1:
        genes = transcript_gene_links1[transcript_id]
        for gene in genes:
            orthologs = ortholog_pairings[gene]
            for ortholog in orthologs:
                if transcript_gene_links2[ortholog]:
                    [transcript_links[transcript_id].append(i) for i in transcript_gene_links2[ortholog] if i in fasta_names]
    # write to output file
    with open(output_file, "w") as outfile:
        for transcript_id in transcript_links:
            [outfile.write("{0}\t{1}\n".format(transcript_id, j)) for j in transcript_links[transcript_id]]


def get_orthologous_pairs(input_bed, input_pairs_file, output_bed):
    """
    Given a list of genes, get the ortholog from a file only if it exists

    Args:
        input_bed (str): path to file containing transcript and gene id links
        input_pairs_file (str): path to file containing specific ortholog links
        output_bed (str): path to output file
    """

    print("Filtering orthologs...")

    gene_ids = [i[1] for i in gen.read_many_fields(input_bed, "\t")]
    ortholog_entries = gen.read_many_fields(input_pairs_file, ",")
    entries_kept = []
    with open(output_bed, "w") as outfile:
        for entry in ortholog_entries:
            # ensure the gene is in the gene ids required and that there is an ortholog
            if entry[1] in gene_ids and entry[4] and entry[1] not in entries_kept:
                outfile.write("{0}\t{1}\n".format(entry[1], entry[4]))
                entries_kept.append(entry[1])
    return entries_kept


def get_purine_matched_sets(focal_set, sets, reverse = None):
    """
    Given a set of sequences, retain those only with the same purine (or pyrimidine)
    content of the focal set

    Args:
        focal_set (list): list of seqeunces
        sets (list): list of lists of test sequences
        reverse (bool): if set, calculate pyrimidine content

    Returns:
        kept_sets (list): list of lists keeping retained sets
    """

    kept_sets = []
    purine_content = calc_purine_content(focal_set, reverse = reverse)
    [kept_sets.append(i) for i in sets if calc_purine_content(i) == purine_content]
    return kept_sets


def get_sequences_synonymous_hits_to_codons(seqs, codon_set):
    """
    Given a sequence, return 1) all the parts of the sequence where the synonymous site in any frame makes a codon
    in the codon_set, 2) all cases where it isnt

    Args:
        seq (str): sequence string
        codon_set (list): list of codons to query

    Returns:
        sorted_seqs (list): [hits, no hits]
    """

    hits = [[],[]]
    no_hits = [[],[]]
    for i in range(0, len(seqs[0])-3, 3):
        if seqs[0][i+1:i+4] in codon_set or seqs[0][i+2:i+5] in codon_set or seqs[1][i+1:i+4] in codon_set or seqs[1][i+2:i+5] in codon_set:
            hits[0].append(seqs[0][i:i+3])
            hits[1].append(seqs[1][i:i+3])
        else:
            no_hits[0].append(seqs[0][i:i+3])
            no_hits[1].append(seqs[1][i:i+3])

    hits = ["".join(i) for i in hits]
    no_hits = ["".join(i) for i in no_hits]

    return [hits, no_hits]


def get_transcript_and_orthologs(input_file1, input_file2, ortholog_transcript_links):
    """
    Given the sequences for a genome and orthologous seqeunces, pair the sequences in a
    dictionary to use elsewhere.

    Args:
        input_file1 (str): path to fasta containing first set of sequences
        input_file2 (str): path to fasta containing second set of sequences
        ortholog_transcript_links (str): path to file containing the links

    Returns:
        sequences (dict): dictionary containing the linked sequences.
            dict[transcript_id] = [[genome1_seq], [..ortholog_seqs]]
    """

    # get the links between transcripts
    links = collections.defaultdict(lambda: [])
    [links[i[0]].append(i[1]) for i in gen.read_many_fields(ortholog_transcript_links, "\t")]
    # read in the sequences
    names, seqs = gen.read_fasta(input_file1)
    ortholog_names, ortholog_seqs = gen.read_fasta(input_file2)
    # now create the output linking the sequences
    sequences = collections.defaultdict(lambda: [[],{}])
    for transcript_id in links:
        sequences[transcript_id][0].append(seqs[names.index(transcript_id)])
        for i in links[transcript_id]:
            sequences[transcript_id][1][i] = ortholog_seqs[ortholog_names.index(i)]
    return sequences


def group_family_results(result_list, families, return_groups = None):
    """
    Group the results of sequences in paralagous family together

    Args:
        result_list (dict): dictionary of results with transcript ids a keys
        families (list): list of lists containing families

    Returns:
        outputs (dict): dictionary containing list of outputs sorted by families,
            dict[id] = [result]
    """

    family_ids = []
    [family_ids.extend(i) for i in families]

    outputs = collections.defaultdict(lambda: [])
    returned_groups = collections.defaultdict(lambda: [])
    for id in result_list:
        if id in family_ids:
            list_id = [i for i, lst in enumerate(families) if id in lst][0]
            outputs["group_{0}".format(list_id)].append(result_list[id])
            if return_groups:
                returned_groups["group_{0}".format(list_id)].append(id)
        else:
            outputs[id].append(result_list[id])
            returned_groups[id].append(id)
    if return_groups:
        return outputs, returned_groups
    else:
        return outputs

def return_median_family_id(result_list, families):
    """
    Group the results of sequences in paralagous family together. Return the
    id of the result with the median value

    Args:
        result_list (dict): dictionary of results with transcript ids a keys
        families (list): list of lists containing families

    Returns:
        outputs (dict): dictionary containing list of outputs sorted by families,
            dict[id] = id
    """

    # get a list of all ids that are part of the families file
    family_ids = []
    [family_ids.extend(i) for i in families]

    seen_ids = []
    outputs = collections.defaultdict(lambda: [[],[]])

    for id in result_list:
        if id in family_ids:
            list_id = [i for i, lst in enumerate(families) if id in lst][0]
            outputs["group_{0}".format(list_id)][0].append(result_list[id])
            outputs["group_{0}".format(list_id)][1].append(id)
        else:
            outputs[id][0].append(result_list[id])
            outputs[id][1].append(id)

    output_ids = {}
    for id in outputs:
        median = np.nanmedian(outputs[id][0])
        if np.isfinite(median):
            if median in outputs[id][0]:
                median_id = outputs[id][1][outputs[id][0].index(median)]
            else:
                sorted_values = sorted(outputs[id][0])
                mid = int(np.divide(len(sorted_values), 2))
                values = [sorted_values[mid-1], sorted_values[+1]]
                chosen_value = np.random.choice(values)
                median_id = outputs[id][1][outputs[id][0].index(chosen_value)]
        else:
            median_id = outputs[id][1][0]
        output_ids[id] = median_id

    return output_ids



def return_median_family_result(result_list, families):
    """
    Group the results of sequences in paralagous family together

    Args:
        result_list (dict): dictionary of results with transcript ids a keys
        families (list): list of lists containing families

    Returns:
        outputs (dict): dictionary containing list of outputs sorted by families,
            dict[id] = [result]
    """

    grouped_results = group_family_results(result_list, families)
    outputs = {i: np.median(grouped_results[i]) for i in grouped_results}

    return outputs


def intersect_coding_exons(full_bed, reduced_bed, output_file):
    """
    Intersect the two files, only keep those entries that have 100% hits for both

    Args:
        full_bed (str): path to file containing all exons
        reduced_bed (str): path to file containing filtered cds exons
        output_file (str): path to output file
    """

    args = ["bedtools", "intersect", "-a", full_bed, "-b", reduced_bed, "-wo", "-f", "1", "-r", "-s"]
    gen.run_process(args, file_for_output = output_file)


def intersect_non_coding_exons(full_bed, reduced_bed, output_file):
    """
    Given a bed file and reudced file, return entries in the full bed file
    that do not appear in the reduced file

    Args:
        full_bed (str): path to bed file containing full entries
        reduced_bed (str): path to bed file containing reduced entries
        output_file (str): path to output file
    """

    args = ["bedtools", "intersect", "-a", full_bed, "-b", reduced_bed, "-wo", "-r", "-s", "-v"]
    gen.run_process(args, file_for_output = output_file)


def list_alignments_to_strings(seq_list):
    """
    For a dictionary containing alignments for each transcript,
    i.e. dict[id] = [seq1, seq2], join each id to form just 2 alignment strings

    Args:
        seq_list (dict): dictionary containing list of alignments, dict[id] = [seq1, seq2]

    Returns:
        alignments (list): list containing [joined_seqs1, joined_seqs2]
    """

    alignments = [[],[]]

    [alignments[0].append(seq_list[i][0]) for i in sorted(seq_list) if len(seq_list[i][0])]
    [alignments[1].append(seq_list[i][1]) for i in sorted(seq_list) if len(seq_list[i][1])]
    alignments = ["".join(i) for i in alignments]
    return alignments


def list_transcript_ids_from_features(gtf_file_path, exclude_pseudogenes=True, full_chr=False, check_chr = True, transcript_id_search_pattern = transcript_id_pattern, ):
    """
    Given a gtf file, return a list of unique transcript ids associated with
    the specific features

    Args:
        gtf_file_path (str): path to input .gtf file
        exclude_pseudogenes (bool): if true, exclude entries labelled as pseudogenes
        full_chr (bool): if true, chromosome is represented as "chrX"

    Returns:
        unique_ids (list): list of unique transcript ids
    """

    print("Getting transcript ids...")

    ids = []

    # get entries, exclude pseudogenes if necessary
    if exclude_pseudogenes:
        gtf_entries = gen.run_process(["grep", "gene_biotype \"protein_coding\"", gtf_file_path])
        entries = [i for i in gtf_entries.split("\n") if len(i)]
    else:
        entries = ["".join(gen.read_many_fields(gtf_file_path, "\t"))]

    # get the index to check for the chromosome
    chr_index = 0
    if full_chr:
        chr_index = 3

    for entry in entries:
        # ensure this is a correct chromosome
        if check_chr and entry[chr_index] in "123456789XY" and entry[chr_index+1] in "0123456789XY\t" or not check_chr:
            # search for the transcript id
            id = re.search(transcript_id_search_pattern, entry)
            # if the id exists, add to list
            if id:
                ids.append(id.group(0))

    # now get a unique set
    unique_ids = sorted(list(set(ids)))

    return unique_ids


def pick_random_family_member(families_file, seqs_dict, output_file = None):
    """
    Given a bed file containing families, pick a random member from the family
    to keep and return the filtered dictionary.

    Args:
        families_file (str): path to bed file containing families
        seqs_dict (dict): dictionary with transcript ids as keys

    Returns:
        seqs_dict (dict): dictionary but only containing one entry per family
    """

    if families_file:
        np.random.seed()
        families = gen.read_many_fields(families_file, "\t")
        new_families = []
        for fam in families:
            fam = [i for i in fam if i in seqs_dict]
            if len(fam):
                new_families.append(fam)
        families = new_families
        # get a random member from each family
        query_ids = [i for i in seqs_dict]
        queried_ids = []
        ids_to_keep = []
        chosen_family_members = []
        for i, id in enumerate(query_ids):
            if id not in queried_ids:
                group = [i for i in families if id in i]
                if len(group):
                    group = group[0]
                    queried_ids.extend(group)
                    choice = np.random.choice(group)
                    ids_to_keep.append(choice)
                    chosen_family_members.append(choice)
                else:
                    queried_ids.append(id)
                    ids_to_keep.append(id)

        seqs_dict = {i: seqs_dict[i] for i in ids_to_keep}

        if output_file:
            with open(output_file, "w") as outfile:
                [outfile.write("{0}\n".format(i)) for i in sorted(chosen_family_members)]

    return seqs_dict


def quality_filter_cds_sequences(input_fasta, output_fasta, stop_codon_set = None):
    """
    Quality filter coding sequences

    Args:
        input_fasta (str): path to fasta file containing input CDS sequences
        output_fasta (str): path to output file to put sequences that pass filtering
    """

    print("Filtering cds...")

    if not stop_codon_set:
        stop_codon_set = stops

    # copile regex searches
    actg_regex = re.compile("[^ACTG]")
    codon_regex = re.compile(".{3}")

    # read the sequences
    names, seqs = gen.read_fasta(input_fasta)

    print("{0} sequences prior to filtering...".format(len(seqs)))

    # filter the sequences
    with open(output_fasta, "w") as outfile:
        pass_count = 0
        for i, name in enumerate(names):
            seq = seqs[i]
            passed = True
            # check to see if the first codon is ATG
            if passed and seq[:3] != "ATG":
                passed = False
            # check to see if the last codon is a stop codon
            if passed and seq[-3:] not in stop_codon_set:
                passed = False
            # check to see if sequence is a length that is a
            # multiple of 3
            if passed and len(seq) % 3 != 0:
                passed = False
            # check if there are any non ACTG characters in string
            non_actg = re.subn(actg_regex, '!', seq)[1]
            if passed and non_actg != 0:
                passed = False
            # check if there are any in frame stop codons
            codons = re.findall(codon_regex, seq[3:-3])
            inframe_stops = [codon for codon in codons if codon in stop_codon_set]
            if passed and len(inframe_stops):
                passed = False
            # only if passed all the filters write to file
            if passed:
                outfile.write(">{0}\n{1}\n".format(name, seq))
                pass_count += 1

    print("{0} sequences after filtering...".format(pass_count))


def remove_terminal_exons(full_bed, intersect_bed, output_file):
    """
    Remove cases where the last exon might be the terminal exon or there are
    incorrectly annotated items.

    Args:
        full_bed (str): path to bed file containing all exons
        intersect_bed (str): path to bed file containing an intersect output
        output_file (str): path to output file
    """

    # get a list of exons, ensure it is a correct entry and get the last one
    full_exons = [i for i in gen.read_many_fields(full_bed, "\t") if len(i) > 3]
    exon_list = collections.defaultdict(lambda: [])
    [exon_list[i[3].split(".")[0]].append(int(i[3].split(".")[1])) for i in full_exons]
    transcript_terminal_exons = {id: max(exon_list[id]) for id in exon_list}
    # read in the intersect cases
    intersect_cases = gen.read_many_fields(intersect_bed, "\t")
    # now filter out any remaining terminal exons
    with open(output_file, "w") as outfile:
        for intersect in intersect_cases:
            intersect_id = intersect[9].split(".")[0]
            intersect_exon_id = int(intersect[9].split(".")[1])
            if intersect_id in transcript_terminal_exons:
                original_id = intersect[3].split(".")[0]
                original_exon_id = int(intersect[3].split(".")[1])

                # ensure a match to the same transcript
                if intersect_id == original_id and intersect_exon_id == original_exon_id:
                    # ensure this is not the first exon
                    if original_exon_id != 1:
                        final_exon = transcript_terminal_exons[original_id]
                        # if the intersect case is not the final case in the full exons file
                        if original_exon_id != final_exon:
                            outfile.write("{0}\n".format("\t".join(gen.stringify(intersect[:6]))))


def sort_bed_file(input_bed, output_bed):
    """
    Sort bed file

    Args:
        input_bed (str): path to bed file to be sorted
        output_bed (str): path to sorted output file
    """

    sort_output = "temp_files/{0}.bed".format(random.random())
    # sort the bed files
    gen.run_process(["sort-bed", input_bed], file_for_output = sort_output)
    # move to the required output file
    gen.run_process(["mv", sort_output, output_bed])
    gen.remove_file(sort_output)



# def get_potential_stops(codon):
#
#     stops = []
#     nts_list = list(codon)
#     for i, codon_nt in enumerate(nts_list):
#         for nt in nucleotides:
#             if nt != codon_nt:
#                 query_list = copy.deepcopy(nts_list)
#                 query_list[i] = nt
#                 print(query_list)
#                 if "".join(query_list) in stops:
#                     stops.append("".join(query_list))
#     print(stops)
#     return stops

def get_one_away_codons():

    codon_list = [i for i in sorted(codon_map)]
    codon_list.extend(stops)
    one_away_codons = collections.defaultdict(lambda: [])

    for codon in codon_list:
        nt_list = list(codon)
        for i, codon_nt in enumerate(nt_list):
            for nt in nucleotides:
                if nt != codon_nt:
                    query_nts = copy.deepcopy(nt_list)
                    query_nts[i] = nt
                    one_away_codons[codon].append("".join(query_nts))

    return one_away_codons


def get_one_away_indicies(seq, one_away_codons):

    kept = []
    # for the first sequence
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        # only do it for the codons that may have one that codes for it after
        if i+3 < len(seq):
            query_codon1 = seq[i+1:i+4]
            query_codon2 = seq[i+2:i+5]
            potential_stops1 = one_away_codons[query_codon1]
            potential_stops2 = one_away_codons[query_codon2]
            # only keep those that might overlap stop and are strictly twofold or fourfold degenerates
            if len(set(potential_stops1) & set(stops)) or len(set(potential_stops2) & set(stops)) and codon in twofold or codon in fourfold:
                kept.append(i)
    return kept


def keep_only_potential_stops(seq1, seq2):

    one_away_codons = get_one_away_codons()
    kept_seq1 = get_one_away_indicies(seq1, one_away_codons)
    kept_seq2 = get_one_away_indicies(seq2, one_away_codons)
    kept_all = sorted(list(set(kept_seq1 + kept_seq2)))

    kept_seq1 = [seq1[i:i+3] for i in kept_all]
    kept_seq2 = [seq2[i:i+3] for i in kept_all]

    print("\n")
    print(seq1)
    print(seq2)
    print("\n")
    print(kept_seq1)
    print(kept_seq2)

    kept_seq1 = "".join(kept_seq1)
    kept_seq2 = "".join(kept_seq2)

    return [kept_seq1, kept_seq2]


def replace_motifs_in_seq(seq, motifs):

    motif_hits = []


    nt_hits = []
    for motif in motifs:
        hits = re.finditer("(?=({0}))".format(motif), seq)
        hits = [i for i in hits]
        for i in hits:
            nt_hits.extend(list(range(i.start(), i.start() + len(motif))))
    # clean up hits
    nt_hits = sorted(list(set(nt_hits)))
    seq_nts = list(seq)
    for hit in nt_hits:
        seq_nts[hit] = "X"
    replaced_seq = "".join(seq_nts)
    return replaced_seq


def replace_motifs_in_seqs(seq_list, motif_set = None):

    if not motif_set:
        print("Please provide a motif set...")
        raise Exception

    replaced_seq_list = {}
    for id in seq_list:
        seq = seq_list[id]
        removed_motifs_seq = replace_motifs_in_seq(seq, motif_set)
        replaced_seq_list[id] = removed_motifs_seq

    return replaced_seq_list



def generate_random_motifs(motif_file, output_dir, required, stop_restricted = None, exclude_motif_set = None, match_gc = None):

    # get all the motifs in the provided motif file
    motif_set = [i[0] for i in gen.read_many_fields(motif_file, "\t") if "#" not in i[0] and ">" not in i[0]]
    # get the possible lengths of the motifs in the provided motif file
    lengths = list(set([len(i) for i in motif_set]))

    # generate all possible choices of motifs with lengths in lengths
    choices = {i: ["".join(j) for j in it.product(nucleotides, repeat = i)] for i in lengths}

    # if the motifs needs to be restricted by stops
    if stop_restricted:
        choices = {i: keep_only_no_stops_in_one_reading_frame(choices[i]) for i in choices}

    # set up the arguments
    required_sets = list(range(required))
    args = [motif_set, choices, output_dir]
    kwargs_dict = {"exclude_motif_set": exclude_motif_set}
    output_files = simopc.run_simulation_function(required_sets, args, generate_motif_set, kwargs_dict = kwargs_dict, sim_run = False)

    # check that none of the sets are the same
    repeated_files = check_repeated_motif_sets(output_files)
    if len(repeated_files):
        repeated = True
        while repeated:
            [gen.remove_file(file) for file in repeated_files]
            required_sets = list(range(len(repeated_files)))
            output_files = output_files = simopc.run_simulation_function(required_sets, args, generate_motif_set, kwargs_dict = kwargs_dict, sim_run = False)
            # now check for repeats, if there are, repeat again
            repeated_files = check_repeated_motif_sets(output_files)
            if not len(repeated_files):
                repeated = False


def check_repeated_motif_sets(filelist):
    """
    Checks to see if the contents of a motif set file match another in the list
    """
    sets = []
    repeat_files = []
    for file in filelist:
        # sort the motifs
        motif_set = sorted([i[0] for i in gen.read_many_fields(file, "\t")])
        # if that motif set isnt in the list, add to list
        if motif_set not in sets:
            sets.append(motif_set)
        # if the motif set is in the list, add to the list to return
        else:
            repeat_files.append(file)
    return repeat_files


def generate_motif_set(required_sets, motif_set, motif_choices, output_dir, exclude_motif_set = None):

    # sort the motifs for later use
    motif_set = sorted(motif_set)

    # if we want to exclude the motifs found in the motif set, remove from the list
    if exclude_motif_set:
        motif_choices = {i: [motif for motif in motif_choices[i] if motif not in motif_set] for i in motif_choices}

    output_files = []

    for i, set in enumerate(required_sets):
        gen.print_parallel_status(i, required_sets)
        np.random.seed()
        set_generated = False
        # whilst we havent generated a full set
        while not set_generated:
            random_set = []
            # for motif in the set
            for motif in motif_set:
                chosen = False
                # choose a random motif that isnt in the set
                while not chosen:
                    choice = np.random.choice(motif_choices[len(motif)])
                    # ensure the motif hasnt already been chosen
                    if choice not in random_set:
                        random_set.append(choice)
                        chosen = True
            # sort the random set and ensure it is not the real motif set
            random_set = sorted(random_set)
            if random_set != sorted(motif_set):
                set_generated = True

        output_file = "{0}/{1}.txt".format(output_dir, random.random())
        output_files.append(output_file)
        with open(output_file, "w") as outfile:
            [outfile.write("{0}\n".format(i)) for i in sorted(random_set)]

    return output_files


def keep_only_no_stops_in_one_reading_frame(seq_list):
    kept_seqs = []
    for i in seq_list:
        if len(list(set(stops) & set(re.findall(codon_pattern, i)))) == 0 or len(list(set(stops) & set(re.findall(codon_pattern, i[1:])))) == 0 or len(list(set(stops) & set(re.findall(codon_pattern, i[2:])))) == 0:
            kept_seqs.append(i)
    return kept_seqs


def read_motifs(input_file):
    return [i[0] for i in gen.read_many_fields(input_file, "\t") if i[0][0] not in [">", "#"]]


def calc_nucleotide_content(seqs):
    content = {}
    nts = list("".join(seqs))
    for nt in sorted(nucleotides):
        content[nt] = np.divide(nts.count(nt), len(nts))
    return content


def randomise_ese_overlaps_wrapper(alignment_file, exons_fasta, ese_file, output_directory, required = None):
    """


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

    if required:
        args = [alignments, ese_file, output_directory]
        simopc.run_simulation_function(list(range(required)), args, randomise_ese_overlaps, sim_run = False)


def randomise_ese_overlaps(simulations, alignments, motif_file, output_directory):

    outputs = []
    if len(simulations):

        for i, simulation in enumerate(simulations):
            gen.print_parallel_status(i, simulations)
            np.random.seed()
            motifs = read_motifs(motif_file)

            output_file = "{0}/{1}.txt".format(output_directory, random.random())
            outputs.append(output_file)

            # for each sequence, get the motif overlaps, then generate random overlaps of the same sizes
            with open(output_file, "w") as outfile:
                for i, id in enumerate(alignments):
                    sequences = alignments[id]
                    overlaps = sorted(list(set(gen.flatten([sequence_overlap_indicies(seq, motifs) for seq in sequences]))))
                    random_overlaps = generate_random_overlaps(overlaps, len(sequences[0]))
                    if len(random_overlaps):
                        outfile.write(">{0}\n{1}\n".format(id, ",".join(gen.stringify(random_overlaps))))
                    else:
                        outfile.write(">{0}\nna\n".format(id))

    return outputs

def generate_random_overlaps(overlap_list, sequence_length):

    random_overlaps = []

    overlap_sets = []
    new_set = []
    for i in overlap_list:
        new_set.append(i)
        if i + 1 not in overlap_list:
            overlap_sets.append(new_set)
            new_set = []

    if len(overlap_list):
        choices = [i for i in list(range(sequence_length - max([len(i) for i in overlap_sets]))) if i not in overlap_list]

        for set in overlap_sets:
            start = np.random.choice(choices)
            end = start + len(set)
            set_range = list(range(start, end))
            if len([i for i in set_range if i in overlap_list]) == 0 and len([i for i in set_range if i in random_overlaps]) == 0:
                random_overlaps.extend(set_range)

        random_overlaps = sorted(random_overlaps)
    return random_overlaps


def retrieve_alignments(input_fasta1, input_fasta2, links_file, output_file):
    """
    Given a list of transcript ID pairs, extract the two sequences that are linked

    Args:
        input_fasta1 (str): path to fasta 1
        input_fasta2 (str): path to fasta 2
        links_file (str): path to the file containing the links
        output_file (str): path to the output file
    """
    # get the first set of sequences
    input_names1, input_seqs1 = gen.read_fasta(input_fasta1)
    inputs1 = {name: input_seqs1[i] for i, name in enumerate(input_names1)}
    # get the second set of sequences
    input_names2, input_seqs2 = gen.read_fasta(input_fasta2)
    inputs2 = {name: input_seqs2[i] for i, name in enumerate(input_names2)}
    # read in the links
    links = gen.read_many_fields(links_file, "\t")
    # now link them together and write to file
    with open(output_file, "w") as outfile:
        for link in links:
            outfile.write(">{0}\n{1},{2}\n".format(link[0], inputs1[link[0]], inputs2[link[1]]))

def calculate_motifs_stop_ds(motif_set_ids, alignments, motif_sets, output_directory):
    """
    Calculate the ds scores for motifs and for the stops in those motifs
    """
    outputs = []
    # the there is a motif set id
    if motif_set_ids:
        for i, id in enumerate(motif_set_ids):
            gen.print_parallel_status(i, motif_set_ids)
            # now set the set of motifs
            motif_set = read_motifs(motif_sets[id])
            # now extract the motif hits for each sequence

            overlap_list = [[], []]
            overlap_stop_list = [[], []]
            overlap_non_stop_list = [[], []]
            non_overlap_list = [[], []]
            non_overlap_stop_list = [[], []]
            non_overlap_non_stop_list = [[], []]

            for alignment_id in alignments:
                alignment_sequences = alignments[alignment_id]
                # predict overlaps to the human sequence
                overlaps = sequence_overlap_indicies(alignment_sequences[0], motif_set)
                overlap_codons, overlap_stops, overlap_non_stops, non_overlaps, non_overlap_stops, non_overlap_non_stops = extract_alignment_overlaps(overlaps, alignment_sequences, stop_overlaps = True)
                # now retain only fourfold degenerate codons
                overlap_codons = retain_fourfolds(overlap_codons)
                overlap_stops = retain_fourfolds(overlap_stops)
                overlap_non_stops = retain_fourfolds(overlap_non_stops)
                non_overlaps = retain_fourfolds(non_overlaps)
                non_overlap_stops = retain_fourfolds(non_overlap_stops)
                non_overlap_non_stops = retain_fourfolds(non_overlap_non_stops)

                # now convert to strings
                overlap_strings = list_to_strings(overlap_codons)
                overlap_stops_strings = list_to_strings(overlap_stops)
                overlap_non_stops_strings = list_to_strings(overlap_non_stops)
                non_overlaps_strings = list_to_strings(non_overlaps)
                non_overlap_stops_strings = list_to_strings(non_overlap_stops)
                non_overlap_non_stops_strings = list_to_strings(non_overlap_non_stops)

                for i in range(len(overlap_strings)):
                    overlap_list[i].append(overlap_strings[i])
                    overlap_stop_list[i].append(overlap_stops_strings[i])
                    overlap_non_stop_list[i].append(overlap_non_stops_strings[i])
                    non_overlap_list[i].append(non_overlaps_strings[i])
                    non_overlap_stop_list[i].append(non_overlap_stops_strings[i])
                    non_overlap_non_stop_list[i].append(non_overlap_non_stops_strings[i])

            # now convert to strings
            overlap_list = list_to_strings(overlap_list)
            overlap_stop_list = list_to_strings(overlap_stop_list)
            overlap_non_stop_list = list_to_strings(overlap_non_stop_list)
            non_overlap_list = list_to_strings(non_overlap_list)
            non_overlap_stop_list = list_to_strings(non_overlap_stop_list)
            non_overlap_non_stop_list = list_to_strings(non_overlap_non_stop_list)

            # calc ds scores
            overlap_ds = cons.calc_ds(overlap_list)
            overlap_stops_ds = cons.calc_ds(overlap_stop_list)
            overlap_non_stops_ds = cons.calc_ds(overlap_non_stop_list)
            non_overlaps_ds = cons.calc_ds(non_overlap_list)
            non_overlap_stops_ds = cons.calc_ds(non_overlap_stop_list)
            non_overlap_non_stops_ds = cons.calc_ds(non_overlap_non_stop_list)

            output_file = "{0}/run_{1}.txt".format(output_directory, id)
            outputs.append(output_file)
            with open(output_file, "w") as outfile:
                if id != "real":
                    id += 1
                output = [overlap_ds, overlap_stops_ds, overlap_non_stops_ds, non_overlaps_ds, non_overlap_stops_ds, non_overlap_non_stops_ds]
                outfile.write("{0},{1}\n".format(id, ",".join(gen.stringify(output))))


    return outputs


def list_to_strings(codon_sets):
    concats = []
    for set in codon_sets:
        concats.append("".join(set))
    return concats


def retain_fourfolds(codon_sets):

    retained = [[], []]
    for i in range(len(codon_sets[0])):
        codon1 = codon_sets[0][i]
        codon2 = codon_sets[1][i]
        if codon1 in fourfold and codon2 in fourfold and codon_map[codon1] == codon_map[codon2]:
            retained[0].append(codon1)
            retained[1].append(codon2)
        elif codon1 in fourfold and codon2 == "---":
            retained[0].append(codon1)
            retained[1].append(codon2)
    return retained


def extract_alignment_overlaps(overlaps, sequences, stop_overlaps = None):
    """
    Given a set of overlap positions, retrieve the parts of the sequences that
    these overlaps correspond to.

    Args:
        overlaps (list): list of indices
        sequences (list): list containing the alignment pair
    """

    overlap_set, non_overlaps = [[], []], [[], []]
    if stop_overlaps:
        overlap_stops, overlap_non_stops, non_overlap_stops, non_overlap_non_stops = [[], []], [[], []],[[], []], [[], []]

    for i in range(0, len(sequences[0]), 3):
        # if the synonymous site overlaps motif
        if i+2 in overlaps:
            overlap_set[0].append(sequences[0][i:i+3])
            overlap_set[1].append(sequences[1][i:i+3])

            # here we have two frames where it might overlap a stop
            # +1 frame, we want the +1 codon to overlap a stop and all sites to be in the motif
            if sequences[0][i+1:i+4] in stops and len(list(set(range(i+1, i+4)) & set(overlaps))) == 3:
                overlap_stops[0].append(sequences[0][i:i+3])
                overlap_stops[1].append(sequences[1][i:i+3])
            # +2 frame, we want the +2 codon to overlap a stop and all sites to be in the motif
            elif sequences[0][i+2:i+5] in stops and len(list(set(range(i+2, i+5)) & set(overlaps))) == 3:
                overlap_stops[0].append(sequences[0][i:i+3])
                overlap_stops[1].append(sequences[1][i:i+3])
            else:
                overlap_non_stops[0].append(sequences[0][i:i+3])
                overlap_non_stops[1].append(sequences[1][i:i+3])
        else:
            # non overlap codons are those with strictly no overlap
            codon_range = list(range(i, i+3))
            if not len(list(set(codon_range) & set(overlaps))):
                non_overlaps[0].append(sequences[0][i:i+3])
                non_overlaps[1].append(sequences[1][i:i+3])
                # if the +1 frame is a stop, and the first nt of the next codon is not motif
                if sequences[0][i+1:i+4] in stops and i+3 not in overlaps:
                    non_overlap_stops[0].append(sequences[0][i:i+3])
                    non_overlap_stops[1].append(sequences[1][i:i+3])
                # if the +2 frame is a stop, and the first two nts of the next codon is not motif
                elif sequences[0][i+2:i+5] in stops and i+3 not in overlaps and i+4 not in overlaps:
                    non_overlap_stops[0].append(sequences[0][i:i+3])
                    non_overlap_stops[1].append(sequences[1][i:i+3])
                else:
                    non_overlap_non_stops[0].append(sequences[0][i:i+3])
                    non_overlap_non_stops[1].append(sequences[1][i:i+3])


    if stop_overlaps:
        return overlap_set, overlap_stops, overlap_non_stops, non_overlaps, non_overlap_stops, non_overlap_non_stops
    else:
        return overlap_set, non_overlaps


def chunk_overlaps(overlaps):
    overlaps = sorted(list(set(overlaps)))
    new_overlaps = []
    new_set = []
    for i in overlaps:
        new_set.append(i)
        if i+1 not in overlaps:
            new_overlaps.append(new_set)
            new_set = []
    return new_overlaps

def return_overlap_motifs(sequence, motifs, inverse = None):
    overlaps = sequence_overlap_indicies(sequence, motifs)
    if inverse:
        overlaps = [i for i in list(range(len(sequence))) if i not in overlaps]
    chunked_overlaps = chunk_overlaps(overlaps)
    overlap_motifs = ["".join([sequence[i] for i in chunk]) for chunk in chunked_overlaps]
    return overlap_motifs


def motif_excess_test(filelist, motifs):

    outputs = []

    if len(filelist):
        for i, file in enumerate(filelist):
            gen.print_parallel_status(i, filelist)
            # read the sequences
            sequences = gen.read_many_fields(file, ",")[0]
            # get the sequence regions
            core_motif_hits = []
            flank_motif_hits = []
            for seq in sequences:
                if len(seq) > 211:
                    flank_sequences = [seq[2:69], seq[-69:2]]
                    midpoint = int(len(seq) / 2)
                    core_sequence = seq[midpoint-33:midpoint+34]
                    core_motif_hits.extend(return_overlap_motifs(core_sequence, motifs))
                    [flank_motif_hits.extend(return_overlap_motifs(seq, motifs)) for seq in flank_sequences]
            # now calculate the density of stops in the set of motifs
            core_stop_density = seqo.calc_motif_density(core_motif_hits, stops)
            flank_stop_density = seqo.calc_motif_density(flank_motif_hits, stops)
            # calculate the excess
            excess = np.divide(flank_stop_density - core_stop_density, core_stop_density) * 100
            # add to the outputs
            outputs.append([file.split("/")[-1].split(".")[0], core_stop_density, flank_stop_density, excess])

    return outputs


def get_sub_rate(seq1, seq2, p = None):
    """
    Given two sequences, calculate the rate at which nucleotides differ between
    between the two sequences.

    Args:
        seq1 (str): sequence 1
        seq2 (str): sequence 2

    Returns:
        sub_rate (list):
    """

    # print(count)
    subs = len([i for i in range(len(seq1)) if seq1[i] in nucleotides and seq2[i] in nucleotides and seq1[i] != seq2[i]])
    total_nts = len([i for i in range(len(seq1)) if seq1[i] in nucleotides and seq2[i] in nucleotides])
    sub_rate = np.divide(subs, total_nts)

    if p:
        print(subs, total_nts, np.divide(subs, total_nts))
    return sub_rate


def get_dinucleotide_substitutions(ids, sequence_list, motifs):
    """
    Given a list of sequence ids and sequences, predict hits, then calculate the
    number of each dinucleotide and subs for each dinucleotide for both the hits
    and non hits

    Args:
        ids (list): list of sequence ids
        sequence_list (dict): sequences
        motifs (list): list of motifs

    Return:
        outputs (list): list of counts and subs for hits and non hits
    """

    outputs = []

    motif_dint_count = collections.defaultdict(lambda: 0)
    motif_dint_sub = collections.defaultdict(lambda: 0)
    non_motif_dint_count = collections.defaultdict(lambda: 0)
    non_motif_dint_sub = collections.defaultdict(lambda: 0)

    if len(ids):
        for i, id in enumerate(ids):
            gen.print_parallel_status(i, ids)

            sequences = sequence_list[id][0]
            human = sequences[0]
            mac = sequences[1]
            # get predicted hits to motifs and chunk to get motifs
            hits = sequence_overlap_indicies(human, motifs)
            non_hits = [i for i in range(len(human)) if i not in hits]
            chunked_hits = chunk_overlaps(hits)
            chunked_non_hits = chunk_overlaps(non_hits)
            # get the motif sequences
            hit_human_motifs = ["".join([human[i] for i in chunk]) for chunk in chunked_hits]
            hit_mac_motifs = ["".join([mac[i] for i in chunk]) for chunk in chunked_hits]
            non_hit_human_motifs = ["".join([human[i] for i in chunk]) for chunk in chunked_non_hits]
            non_hit_mac_motifs = ["".join([mac[i] for i in chunk]) for chunk in chunked_non_hits]

            # for each motif hit in the hit motifs, go through each position,
            # get the dinucleotide and see whether it is different
            for i, hit in enumerate(hit_human_motifs):
                human_motif = hit
                mac_motif = hit_mac_motifs[i]

                for pos in range(len(human_motif)-1):
                    if human_motif[pos] in nucleotides and mac_motif[pos] in nucleotides:
                        motif_dint_count[human_motif[pos:pos+2]] += 1
                        if human_motif[pos:pos+2] != mac_motif[pos:pos+2]:
                            motif_dint_sub[human_motif[pos:pos+2]] += 1
            # for each non hit in the non hit motifs, go through each position,
            # get the dinucleotide and see whether it is different
            for i, non_hit in enumerate(non_hit_human_motifs):
                human_motif = non_hit
                mac_motif = non_hit_mac_motifs[i]
                for pos in range(len(human_motif)-1):
                    if human_motif[pos] in nucleotides and mac_motif[pos] in nucleotides:
                        non_motif_dint_count[human_motif[pos:pos+2]] += 1
                        if human_motif[pos:pos+2] != mac_motif[pos:pos+2]:
                            non_motif_dint_sub[human_motif[pos:pos+2]] += 1

    #unpickle
    motif_dint_count = unpickle(motif_dint_count)
    motif_dint_sub = unpickle(motif_dint_sub)
    non_motif_dint_count = unpickle(non_motif_dint_count)
    non_motif_dint_sub = unpickle(non_motif_dint_sub)
    output = [motif_dint_count, motif_dint_sub, non_motif_dint_count, non_motif_dint_sub]
    return output

def unpickle(list):
    return {i: list[i] for i in list}

def calc_dinucleotide_substitution_rates(subs, counts):
    """
    Calculate the substitution rates given the sub counts and total counts.

    Args:
        subs (dict): counts of substitutions for each dint
        counts (dict): total counts for each dint

    Returns:
        sub_rates (dict): the sub rates
    """
    sub_rates = {}
    for dint in sorted(counts):
        sub_rates[dint] = np.divide(subs[dint], counts[dint])
    return sub_rates


def calc_dint_chisquare(subs, counts, group1, group2):

    total_counts = sum(counts.values())
    total_subs = sum(subs.values())
    group1_subs = sum([subs[i] for i in subs if i in group1])
    group2_subs = sum([subs[i] for i in subs if i in group2])
    group1_totals = sum([counts[i] for i in counts if i in group1])
    group2_totals = sum([counts[i] for i in counts if i in group2])
    expected_subs_group1 = np.divide(group1_totals, total_counts)*total_subs
    expected_subs_group2 = np.divide(group2_totals, total_counts)*total_subs
    chisquare = scipy.stats.chisquare([group1_subs, group2_subs], f_exp = [expected_subs_group1, expected_subs_group2])
    output = [group1_totals, group1_subs, expected_subs_group1, group2_totals, group2_subs, expected_subs_group2, total_counts, total_subs, chisquare]
    return output


def calc_overlaps(runs, motifs, exon_list, families_file):
    """
    Wrapper to ask whether there are a disproportionate number of motifs that overlap
    that are not made up of motifs containing stop codons.

    Args:
        runs (list): list of runs to iterate over
        motifs (list): list of motifs
        exon_list (dict): dictionary of exon sequences
        families_file (str): path to famililies file

    Returns:
        outputs (list): list containing results
    """

    outputs = []

    if len(runs):
        stop_motifs = [i for i in motifs if len(re.findall("(?=(TAA|TAG|TGA))", i))]
        non_stop_motifs = [i for i in motifs if i not in stop_motifs]

        for i, run in enumerate(runs):
            np.random.seed()
            gen.print_parallel_status(i, runs)
            # get a list of random family members
            data_exon_list = sequo.pick_random_family_member(families_file, exon_list)
            temp_exon_list = []
            for id in data_exon_list:
                for seq in data_exon_list[id]:
                    temp_exon_list.append(seq)
            data_exon_list = temp_exon_list

            stop_motif_hits = gen.flatten([sequo.return_overlap_motifs(seq, stop_motifs) for seq in data_exon_list])
            non_stop_motif_hits = gen.flatten([sequo.return_overlap_motifs(seq, non_stop_motifs) for seq in data_exon_list])
            total_hits = len(stop_motif_hits) + len(non_stop_motif_hits)

            # anything longer than 6 nucleotides will be an overlap
            stop_overlap_hits = [i for i in stop_motif_hits if len(i) > 6]
            non_stop_overlap_hits = [i for i in non_stop_motif_hits if len(i) > 6]
            total_overlapping_hits = len(stop_overlap_hits) + len(non_stop_overlap_hits)

            # get the proporition of total hits that are overlapping
            proportion_overlapping = np.divide(total_overlapping_hits, total_hits)

            # now calcaulte the expected number of overlaps
            expected_stop_overlaps = proportion_overlapping * len(stop_motif_hits)
            expected_non_stop_overlaps = proportion_overlapping * len(non_stop_motif_hits)

            # do the chisquare test
            chi = chisquare([len(stop_overlap_hits), len(non_stop_overlap_hits)], f_exp=[expected_stop_overlaps, expected_non_stop_overlaps])
            outputs.append([expected_stop_overlaps, len(stop_overlap_hits), expected_non_stop_overlaps, len(non_stop_overlap_hits), np.mean([len(i) for i in stop_overlap_hits]), np.mean([len(i) for i in non_stop_overlap_hits]), chi])

    return outputs


def calc_motif_overlap_density(runs, motifs, exon_list, families_file, sim = None):
    """
    Wrapper to ask whether overlapping hits to motifs containing more stop codons
    than hits that dont overlap.

    Args:
        runs (list): list of runs to iterate over
        motifs (list): list of motifs
        exon_list (dict): dictionary of exon sequences
        families_file (str): path to famililies file
        sim (bool): if true, pick a random set of motifs to use as the overlaps

    Returns:
        outputs (list): list containing results
    """

    outputs = []
    if len(runs):
        np.random.seed()
        for i, run in enumerate(runs):
            gen.print_parallel_status(i, runs)
            # pick a random family member
            data_exon_list = sequo.pick_random_family_member(families_file, exon_list)

            temp_exon_list = []
            for id in data_exon_list:
                for seq in data_exon_list[id]:
                    temp_exon_list.append(seq)
            data_exon_list = temp_exon_list

            motif_hits = gen.flatten([sequo.return_overlap_motifs(seq, motifs) for seq in data_exon_list])

            if sim:
                #if a simulation, pick a random set of hits to be the overlaps
                overlap_hit_length = len([i for i in motif_hits if len(i) > 6])
                overlap_hits = np.random.choice(motif_hits, overlap_hit_length, replace = False)
                non_overlap_hits = motif_hits
            else:
                # otherwise, use real overlaps
                overlap_hits = [i for i in motif_hits if len(i) > 6]
                non_overlap_hits = [i for i in motif_hits if len(i) == 6]

            # calculat the stop codon density in those hits
            stop_non_overlaps = seqo.calc_motif_density(overlap_hits, stops)
            stop_overlaps = seqo.calc_motif_density(non_overlap_hits, stops)

            outputs.append([stop_non_overlaps, stop_overlaps])

    return outputs


def calc_subs(motif1, motif2):
    count = 0
    for i, nt in enumerate(list(motif1)):
        if motif2[i] != nt:
            count += 1
    return count

def get_subs(motifs):
    checked_pairs = []
    subs = []
    for motif1 in motifs:
        for motif2 in motifs:
            if motif1 != motif2:
                if [motif1, motif2] not in checked_pairs and [motif2, motif1] not in checked_pairs:
                    subs.append(calc_subs(motif1, motif2))
                    checked_pairs.append([motif1, motif2])
    return subs


def get_function_frames(motifs):

    non_function_frames = {0: [], 1: [], 2: []}

    for motif in motifs:
        if motif[0:3] in stops or motif[3:6] in stops:
            non_function_frames[0].append(motif)
        elif motif[1:4] in stops:
            non_function_frames[2].append(motif)
        elif motif[2:5] in stops:
            non_function_frames[1].append(motif)

    function_frames = {i: len(motifs) - len(non_function_frames[i]) for i in non_function_frames}
    return function_frames

def calc_hits(file_ids, filelist, exon_list, exon_start_indices):

    outputs = {}

    if len(file_ids):
        for file_no, id in enumerate(file_ids):
            gen.print_parallel_status(file_no, file_ids)
            motifs = read_motifs(filelist[id])
            stop_motifs = [i for i in motifs if len(re.findall("(?=TAA|TAG|TGA)", i)) > 0]
            non_stop_motifs = [i for i in motifs if i not in stop_motifs]

            function_frames = get_function_frames(stop_motifs)

            stop_hits = collections.defaultdict(lambda: [])
            non_stop_hits = collections.defaultdict(lambda: [])

            for i, exon in enumerate(exon_list):
                for frame in [0,1,2]:
                    hits = seqo.calc_motif_count_frame([exon], motifs, frame = frame)
                    true_frame = (exon_start_indices[i] + frame) % 3
                    stop_hits[true_frame].extend([i for i in hits if i in stop_motifs])
                    non_stop_hits[true_frame].extend([i for i in hits if i in non_stop_motifs])

            stop_hits = {i: len(stop_hits[i]) for i in stop_hits}
            non_stop_hits = {i: len(non_stop_hits[i]) for i in non_stop_hits}
            outputs[id] = [stop_hits, non_stop_hits, function_frames, len(stop_motifs), len(non_stop_motifs)]

    return outputs



def calc_hits_lincrna(file_ids, filelist, exon_list, total_bp):

    outputs = {}

    exons = []
    [exons.extend(exon_list[i]) for i in exon_list]

    if len(file_ids):
        for file_no, id in enumerate(file_ids):
            gen.print_parallel_status(file_no, file_ids)
            motifs = read_motifs(filelist[id])
            motif_search = re.compile("(?=({0}))".format("|".join(motifs)))

            stop_motifs = [i for i in motifs if len(re.findall("(?=TAA|TAG|TGA)", i)) > 0]
            non_stop_motifs = [i for i in motifs if i not in stop_motifs]

            all_hits = []
            for transcript_id in exon_list:
                exons = exon_list[transcript_id]
                for exon in exons:
                    matches = re.finditer(motif_search, exon)
                    hit_motifs = [match.group(1) for match in matches]
                    all_hits.extend(hit_motifs)

            stop_hits = len([i for i in all_hits if i in stop_motifs])
            non_stop_hits = len([i for i in all_hits if i in non_stop_motifs])

            # print(len(non_stop_hits))


            # print(id, stop_hits, total_bp)

            stop_hits = 1000*np.divide(stop_hits, total_bp)
            non_stop_hits = 1000*np.divide(non_stop_hits, total_bp)
            all_hits = 1000*np.divide(len(all_hits), total_bp)


            stop_normalised = np.divide(stop_hits, len(stop_motifs))
            non_stop_normalised = np.divide(non_stop_hits, len(non_stop_motifs))
            all_hits_normalised = np.divide(all_hits, len(motifs))
            outputs[id] = [id, stop_hits, non_stop_hits, stop_normalised, non_stop_normalised, len(stop_motifs), len(non_stop_motifs), all_hits, all_hits_normalised]

    return outputs



def analyse_overlaps(file_ids, filelist, exons):

    outputs = {}
    if len(file_ids):
        for i, id in enumerate(file_ids):
            gen.print_parallel_status(i, file_ids)
            motifs, stop_motifs, non_stop_motifs = get_stop_non_stop_motifs(filelist[id])
            overlaps = calc_overlap_potential(motifs)
            sim_stop, sim_non_stop, sim_diff = calc_overlap_means(overlaps, stop_motifs, non_stop_motifs)
            outputs[id] = [sim_stop, sim_non_stop, sim_diff]
    return outputs



def check_for_overlap(seq, i, motifs):
    test_motifs = [
        seq[i-6:i],
        seq[i-5:i+1],
        seq[i-4:i+2],
        seq[i-3:i+3],
        seq[i-2:i+4],
        seq[i-1:i+5],
        seq[i+1:i+7],
        seq[i+2:i+8],
        seq[i+3:i+9],
        seq[i+4:i+10],
        seq[i+5:i+11],
        seq[i+6:i+12]
    ]

    overlaps = [motif for motif in test_motifs if motif in motifs]
    if len(overlaps):
        return True
    else:
        return False


def calc_sequence_overlap_diffs(file_ids, filelist, exons):

    exons = [i for i in exons if len(i) > 211]
    outputs = {}
    if len(file_ids):
        for i, id in enumerate(file_ids):
            gen.print_parallel_status(i, file_ids)
            motifs = read_motifs(filelist[id])

            overlap_hits = []
            non_overlap_hits = []
            for exon in exons:
                for pos in range(0, len(exon)-6):
                    focal_motif = exon[pos:pos+6]
                    if focal_motif in motifs:
                        overlap = check_for_overlap(exon, pos, motifs)
                        if overlap:
                            overlap_hits.append(focal_motif)
                        else:
                            non_overlap_hits.append(focal_motif)



            # calculate the stop codon density in those hits
            stop_overlaps = seqo.calc_motif_density(overlap_hits, stops)
            stop_non_overlaps = seqo.calc_motif_density(non_overlap_hits, stops)

            overlap_stops = len([i for i in overlap_hits if len(re.findall("(?=(TAA|TAG|TGA))", i)) > 0])
            non_overlap_stops = len([i for i in non_overlap_hits if len(re.findall("(?=(TAA|TAG|TGA))", i)) > 0])
            # overlap_non_stops = len([i for i in overlap_hits if len(re.findall("(?=(TAA|TAG|TGA))", i)) == 0])
            # non_overlap_non_stops = len([i for i in non_overlap_hits if len(re.findall("(?=(TAA|TAG|TGA))", i)) == 0])

            # print(np.divide(overlap_non_stops, len(overlap_hits)))
            # print(np.divide(non_overlap_non_stops, len(non_overlap_hits)))

            expected_overlap = np.divide(non_overlap_stops, len(non_overlap_hits))*len(overlap_hits)
            chi = scipy.stats.chisquare([overlap_stops, len(overlap_hits) - overlap_stops], f_exp = [expected_overlap, len(overlap_hits) - expected_overlap])

            print(overlap_stops, len(overlap_hits))
            print(non_overlap_stops, len(non_overlap_hits))
            print(np.divide(overlap_stops, len(overlap_hits)))
            print(np.divide(non_overlap_stops, len(non_overlap_hits)))
            print(chi)
            # print(id, stop_overlaps, stop_non_overlaps, chi)

    return outputs


def calc_overlap_potential(motifs):
    overlaps = collections.defaultdict(lambda: [])
    for motif in sorted(motifs):
        for overlap_motif in motifs:

            if overlap_motif[-2:] == motif[:2]:
                overlaps[motif].append(overlap_motif)
            if overlap_motif[-3:] == motif[:3]:
                overlaps[motif].append(overlap_motif)
            if overlap_motif[-4:] == motif[:4]:
                overlaps[motif].append(overlap_motif)
            if overlap_motif[-5:] == motif[:5]:
                overlaps[motif].append(overlap_motif)
            if motif[1:] == overlap_motif[:5]:
                overlaps[motif].append(overlap_motif)
            if motif[2:] == overlap_motif[:4]:
                overlaps[motif].append(overlap_motif)
            if motif[3:] == overlap_motif[:3]:
                overlaps[motif].append(overlap_motif)
            if motif[4:] == overlap_motif[:2]:
                overlaps[motif].append(overlap_motif)

    overlaps = {i: len(list(set(overlaps[i]))) for i in overlaps}
    return overlaps


def calc_overlap_means(overlaps, stop_motifs, non_stop_motifs):
    stop_overlaps = {i: overlaps[i] for i in overlaps if i in stop_motifs}
    non_stop_overlaps = {i: overlaps[i] for i in overlaps if i in non_stop_motifs}
    stop_mean = np.mean([stop_overlaps[i] for i in stop_overlaps])
    non_stop_mean = np.mean([non_stop_overlaps[i] for i in non_stop_overlaps])
    diff = non_stop_mean - stop_mean

    return stop_mean, non_stop_mean, diff

def get_stop_non_stop_motifs(motif_file):
    motifs = read_motifs(motif_file)
    stop_motifs = [i for i in motifs if len(re.findall("(?=(TAA|TAG|TGA))", i)) > 0]
    non_stop_motifs = [i for i in motifs if i not in stop_motifs]
    return motifs, stop_motifs, non_stop_motifs


def calc_upstream_sim_densities(sims, sim_seqs):
    outputs = []
    if len(sims):
        for sim in sims:
            sim_seqs = randomise_seqs(sim_seqs)
            outputs.append(seqo.calc_motif_density(sim_seqs, stops))
    return outputs

def randomise_seqs(seqs):
    new_seqs = []
    for seq in seqs:
        nts = list(seq)
        np.random.shuffle(nts)
        new_seqs.append("".join(nts))
    return new_seqs

def mask_seq(seq, motifs, mask_character = "C"):
    """
    Given a seq and some motifs, mask out

    Args:
        seq (str): sequence to mask
        motifs (list): list of motifs to mask
        mask_character (str): if set, the character to replace with

    Returns:
        masked_seq (str): masked sequence
    """
    overlaps = sequence_overlap_indicies(seq, motifs)
    masked_seq = "".join([nt if i not in overlaps else mask_character for i, nt in enumerate(list(seq))])
    return masked_seq


def locate_random_motifs(iterations, seq, motifs, output_directory):
    """
    Given a sequence, locate a set amount of random motifs

    Args:
        iterations (list): list of items to iterate over
        seq (str): the sequence to pick motifs frmo
        motifs (list): list of motifs to simulate
        output_directory (str): path to output directory of outputs
    """
    # get a list of positions in the string
    choices = list(range(len(seq)))
    output_files = []

    choice_length = list(set([len(i) for i in motifs]))[0]
    # for each of the iterations
    if len(iterations):
        for i, iteration in enumerate(iterations):
            np.random.seed()
            gen.print_parallel_status(i, iterations)

            test_set = []
            generated = False
            # if not generated
            while not generated:
                required = len(motifs) - len(test_set)
                # choose a random position
                chosen = np.random.choice(choices, len(motifs) - len(test_set))
                # get the motif and check whether it is not already in set
                test_set.extend(list(set([seq[i:i+choice_length] for i in chosen if "X" not in seq[i:i+choice_length] and seq[i:i+choice_length] not in test_set])))
                if len(test_set) == len(motifs):
                    generated = True

            # write motifs to file
            output_file = "{0}/{1}.txt".format(output_directory, random.random())
            with open(output_file, "w") as outfile:
                [outfile.write("{0}\n".format(i)) for i in test_set]

def calc_motif_hit_purine(ids, seq_list, motifs):
    """
    Calculate the purine content in and out of motif hits

    Args:
        ids (list): list of sequence ids
        seq_list (dict): sequences
        motifs (list): list of motifs

    Returns:
        outputs (dict): purine contents for each seq
    """
    outputs = {}

    if len(ids):
        for i, id in enumerate(ids):
            gen.print_parallel_status(i, ids)

            motif_sequence = []
            non_motif_sequence = []
            flank_sequence_motif = []
            flank_sequence_non_motif = []
            for sequence in exons[id]:
                overlaps = sequo.sequence_overlap_indicies(sequence, motifs)
                motif_sequence.append("".join([sequence[i] for i in overlaps]))
                non_motif_sequence.append("".join([sequence[i] for i in range(len(sequence)) if i not in overlaps]))

                flank_sequence = "".join([sequence[:50], sequence[-50:]])
                flank_sequence_motif.append("".join([flank_sequence[i] for i in overlaps if i in range(len(flank_sequence))]))
                flank_sequence_non_motif.append("".join([flank_sequence[i] for i in range(len(flank_sequence)) if i not in overlaps]))

            outputs[id] = [sequo.calc_purine_content("".join(motif_sequence)), sequo.calc_purine_content("".join(non_ese_sequence)), sequo.calc_purine_content("".join(flank_sequence_ese)), sequo.calc_purine_content("".join(flank_sequence_non_ese))]

    return outputs

def remove_motifs(ids, seqs, motifs):
    """
    Given a list of sequence ids, remove the motifs and return the remaining
    sequence.

    Args:
        ids (list): list of seq ids
        seqs (dict): sequences
        motifs (list): list of motifs to remove

    Returns:
        kept (dict): the remaining bits of sequence
    """

    kept = {}
    if len(ids):
        for i, id in enumerate(ids):
            gen.print_parallel_status(i, ids)
            kept_seqs = []
            for seq in seqs[id]:
                overlaps = get_motifs_overlap_indices([seq], motifs)
                kept_seq = "".join([nt for i, nt in enumerate(list(seq)) if i not in overlaps])
                kept_seqs.append(kept_seq)
            kept[id] = kept_seqs
    return kept
