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
from conservation import Alignment_Functions
import copy
import sys
from useful_motif_sets import nucleotides, stops, codon_map, twofold, fourfold, one_away_codons
from progressbar import ProgressBar
import multiprocessing as mp
import itertools as it

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

def extract_motif_sequences_from_alignment(alignment_seqs, motif_set, reverse = None, only_fourfold = True):
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
            overlap_indices = get_motifs_overlap_indices(alignment_set, motif_set)
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


            # now get the indices that create a stop from the overlap set
            stop_overlap_codon_indices = []
            non_stop_overlap_codon_indices = []
            for chunk in overlap_codon_indices:
                codons = [alignment_set[0][chunk[1]:chunk[-1]+2], alignment_set[0][chunk[2]:chunk[-1]+3], alignment_set[1][chunk[1]:chunk[-1]+2], alignment_set[1][chunk[2]:chunk[-1]+3]]
                if list(set(stops) & set(codons)):
                    stop_overlap_codon_indices.append(chunk)
                else:
                    non_stop_overlap_codon_indices.append(chunk)



            overlap_sequences = [[sequence[i[0]:i[-1]+1] if i != "-" else "---" for i in overlap_codon_indices] for sequence in alignment_set]
            non_overlap_sequences = [[sequence[i[0]:i[-1]+1] if i != "-" else "---" for i in non_overlap_codon_indices] for sequence in alignment_set]
            stop_overlap_sequences = [[sequence[i[0]:i[-1]+1] if i != "-" else "---" for i in stop_overlap_codon_indices] for sequence in alignment_set]
            non_stop_overlap_sequences = [[sequence[i[0]:i[-1]+1] if i != "-" else "---" for i in non_stop_overlap_codon_indices] for sequence in alignment_set]

            motif_overlap_sequences = ["".join(i) for i in overlap_sequences]
            motif_non_overlap_sequences = ["".join(i) for i in non_overlap_sequences]
            motif_stop_overlap_sequences = ["".join(i) for i in stop_overlap_sequences]
            motif_non_stop_overlap_sequences = ["".join(i) for i in non_stop_overlap_sequences]

            motif_overlap_sequences = ["".join(["GC{0}".format(seq[i+2]) if "-" not in seq[i:i+3] else seq[i:i+3] for i in range(0, len(seq), 3)]) for seq in motif_overlap_sequences]
            motif_non_overlap_sequences = ["".join(["GC{0}".format(seq[i+2]) if "-" not in seq[i:i+3] else seq[i:i+3] for i in range(0, len(seq), 3)]) for seq in motif_non_overlap_sequences]
            motif_stop_overlap_sequences = ["".join(["GC{0}".format(seq[i+2]) if "-" not in seq[i:i+3] else seq[i:i+3] for i in range(0, len(seq), 3)]) for seq in motif_stop_overlap_sequences]
            motif_non_stop_overlap_sequences = ["".join(["GC{0}".format(seq[i+2]) if "-" not in seq[i:i+3] else seq[i:i+3] for i in range(0, len(seq), 3)]) for seq in motif_non_stop_overlap_sequences]

            retained_overlap_sequences[id] = motif_overlap_sequences
            retained_non_overlap_sequences[id] = motif_non_overlap_sequences
            retained_stop_overlap_sequences[id] = motif_stop_overlap_sequences
            retained_non_stop_overlap_sequences[id] = motif_non_stop_overlap_sequences

    return retained_overlap_sequences, retained_non_overlap_sequences, retained_stop_overlap_sequences, retained_non_stop_overlap_sequences

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
            retained_overlap_sequence, retained_non_overlap_sequence, retained_stop_overlap_sequences, reatined_non_stop_overlap_sequences = extract_motif_sequences_from_alignment(sequence_alignments, motif_set)

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
            mutation_set_stops = [len([i for i in motif_set if len(re.findall("(TAA|TAG|TGA)", i))]) for motif_set in non_stop_mutation_set]

            non_stop_non_stop_mutation = [motif for i, motif in enumerate(non_stop_motif_set) if mutation_set_stops[i] == 0]
            non_stop_stop_mutation = [motif for i, motif in enumerate(non_stop_motif_set) if mutation_set_stops[i] != 0]


            # # now get all the parts of the sequences where only the motif sets occur
            # full_extracted_sequences = extract_motif_sequences_from_alignment(sequence_alignments, motif_set)
            # non_motif_extracted_sequences = extract_motif_sequences_from_alignment(sequence_alignments, motif_set, reverse = True)
            no_stop_mutation = extract_motif_sequences_from_alignment(sequence_alignments, non_stop_non_stop_mutation)
            stop_mutation = extract_motif_sequences_from_alignment(sequence_alignments, non_stop_stop_mutation)

            # # now calculate the ds scores for each set
            # all_ese_alignments = list_alignments_to_strings(full_extracted_sequences)
            # non_motif_alignments = list_alignments_to_strings(non_motif_extracted_sequences)
            stops_alignments = list_alignments_to_strings(stop_mutation)
            non_stops_alignments = list_alignments_to_strings(no_stop_mutation)

            # # calculate the ds scores
            # all_ese_ds = cons.calc_ds(all_ese_alignments)
            # non_ese_ds = cons.calc_ds(non_motif_alignments)
            stops_ds = cons.calc_ds(stops_alignments)
            non_stops_ds = cons.calc_ds(non_stops_alignments)


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
    alignment_functions = Alignment_Functions(muscle_exe)

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
    alignment_functions = Alignment_Functions(muscle_exe)
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


def group_family_results(result_list, families):
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
    for id in result_list:
        if id in family_ids:
            list_id = [i for i, lst in enumerate(families) if id in lst][0]
            outputs["group_{0}".format(list_id)].append(result_list[id])
        else:
            outputs[id].append(result_list[id])

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
    [alignments[0].append(seq_list[i][0]) for i in sorted(seq_list)]
    [alignments[1].append(seq_list[i][1]) for i in sorted(seq_list)]
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
            [outfile.write("{0}\n".format(i)) for i in chosen_family_members]

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
