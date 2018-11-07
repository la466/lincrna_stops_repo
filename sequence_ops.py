import generic as gen
import file_ops as fo
from regex_patterns import exon_number_pattern, gene_id_pattern, transcript_id_pattern
import re
import os
import collections

class Genome_Functions(object):

    def __init__(self, features_file_path, genome_file_path):
        self.features_file_path = features_file_path
        self.genome_file_path = genome_file_path


    def build_cds(self, cds_features_bed, output_file, clean_run = None):
        """
        Build CDS sequences from the list of features.
        """

        if not os.path.isfile(output_file) or clean_run:
            print("Building cds...")

            stop_codons = ["TAA", "TAG", "TGA"]

            temp_dir = "temp_dir"
            gen.create_output_directories(temp_dir)

            # create temp file to contain sequences
            temp_file = "{0}/temp_cds_features.fasta".format(temp_dir)
            fo.fasta_from_intervals(cds_features_bed, temp_file, self.genome_file_path, names=True)

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


    def generate_genome_dataset(self, dataset_name, dataset_features_output_bed, input_list = None, filter_transcripts = True, clean_run = None):
        """
        Create a dataset containing all the relevant genome features from the
        given list.

        Args:
            dataset_name (str): name for the dataset
            input_list (list): if set, use as a list to filter transcripts from
            filter_transcripts (bool): if true, filter transcripts from list
        """

        print("Generating genome dataset...")

        # set some variables
        self.dataset_name = dataset_name
        self.features = []

        # only run if the output file doesnt exist, or it is a clean run
        if not os.path.isfile(dataset_features_output_bed) or clean_run:

            if filter_transcripts:
                if input_list:
                    self.ids = input_list
                else:
                    print("Input list required...")
                    raise Exception
            else:
                self.ids = "foo"

            # do the filtering
            args = [self.features_file_path]
            # reduce the number of workers here otherwise it doesnt like it
            features = gen.run_parallel_function(self.ids, args, filter_gtf_features, parallel=True, workers = int(os.cpu_count()/2) - 1)
            # output the features to the bed file
            with open(dataset_features_output_bed, "w") as outfile:
                outfile.write("# {0} lines in {1}\n".format(dataset_name, self.features_file_path.split("/")[-1]))
                [outfile.write("{0}\n".format("\t".join(feature))) for feature in features]

            print("{0} dataset created...".format(dataset_name))


    def get_cds_features(self):
        """
        Get all the CDS features
        """

        print("Getting cds features...")

        stop_codons = self.get_stop_codon_features()
        cds_features_list = collections.defaultdict(lambda: [])
        # get a list of all the exon parts that contribute to the cds
        # [cds_features_list[feature[3]].append(feature) for feature in self.features if feature[-1] == "CDS" and feature[3] in self.ids]
        for feature in self.features:
            if feature[-1] == "CDS" and feature[3] in self.ids:
                print(feature)
                cds_features_list[feature[3]].append(feature)

        for id in cds_features_list:
            # get the stop codon coordinates if they exist
            stop_codon = stop_codons[id]
            stop_codon_coordinates = [[i[1], i[2]] for i in stop_codon]
            # get the cds coordinates if they arent in the stop codon list
            cds_features_list[id] = [i for i in cds_features_list[id] if [i[1], i[2]] not in stop_codon_coordinates]
            # now we need to sort the exons in order, reversing if on the minus strand
            strand = cds_features_list[id][0][5]
            if strand == "+":
                cds_features_list[id] = sorted(cds_features_list[id], key = lambda x:x[1])
            elif strand == "-":
                cds_features_list[id] = sorted(cds_features_list[id], key = lambda x:x[2], reverse = True)
            # add the stop codon if it exists
            for stop_codon in stop_codons[id]:
                cds_features_list[id].append(stop_codon)

        return cds_features_list


    def get_stop_codon_features(self):
        """
        Get all the features that match stop codon
        """

        print("Getting stop_codon features...")

        feature_list = collections.defaultdict(lambda: [])
        for feature in self.features:
            if feature[-1] == "stop_codon":
                feature_list[feature[3]].append(feature)
        return feature_list


    def get_transcript_features(self):
        """
        Get the list of transcript features
        """

        print("Getting transcript features...")

        feature_list = collections.defaultdict(lambda: [])
        for feature in self.features:
            if feature[-1] == "transcript" and feature[3] in self.ids:
                feature_list[feature[3]].append(feature)
        return feature_list


    def load_genome_dataset(self, dataset_name, dataset_file_path, input_list = None):
        """
        Load a genome dataset from file.
        """

        try:
            entries = gen.read_many_fields(dataset_file_path, "\t")
            self.features = [[i[0], int(i[1]), int(i[2]), i[3], i[4], i[5], i[6], i[7]] for i in entries if len(i) and i[0][0] != "#"]
            self.dataset_name = dataset_name
            self.dataset_features_file_path = dataset_file_path
            if input_list:
                self.ids = input_list
            else:
                self.ids = list(set([i[3] for i in self.features]))

            print("{0} dataset loaded...".format(dataset_name))

        except FileNotFoundError:
            raise FileNotFoundError


def filter_cds(input_fasta, output_fasta, clean_run = None):
    """
    Quality filter coding sequences

    Args:
        input_fasta (str): path to fasta file containing input CDS sequences
        output_fasta (str): path to output file to put sequences that pass filtering
    """

    if not os.path.isfile(output_fasta) or clean_run:
        print("Filtering cds...")

        # copile regex searches
        actg_regex = re.compile("[^ACTG]")
        codon_regex = re.compile(".{3}")

        stop_codons = ["TAA", "TAG", "TGA"]

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
                if passed and seq[-3:] not in stop_codons:
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
                inframe_stops = [codon for codon in codons if codon in stop_codons]
                if passed and len(inframe_stops):
                    passed = False
                # only if passed all the filters write to file
                if passed:
                    outfile.write(">{0}\n{1}\n".format(name, seq))
                    pass_count += 1

        print("{0} sequences after filtering...".format(pass_count))


def filter_gtf_features(input_list, gtf_file_path, filter_transcripts = True):
    """
    Given a .gtf file, filter the entries based on an input list
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
            transcript_search_results = re.search(transcript_id_pattern, entry_info)
            try:
                current_transcript_id = transcript_search_results.group(0)
                # if filtering by ids, do that here
                passed_filter = True
                if filter_transcripts:
                    passed_filter = current_transcript_id in input_list
                # if it passes the filter
                if passed_filter:
                    # look for gene id
                    gene_search_results = re.search(gene_id_pattern, entry_info)
                    try:
                        current_gene_id = gene_search_results.group(0)
                    except AttributeError:
                        pass
                    exon_number_results = re.search(exon_number_pattern, entry_info)
                    # try to get the exon number
                    try:
                        exon_number = exon_number_results.group(1)
                    except AttributeError:
                        exon_number = "nan"
                    # minus 1 for the start because gtf are in base 0, but want base 1
                    outputs.append([entry[0], str(int(entry[3])-1), entry[4], current_transcript_id, exon_number, entry[6], current_gene_id, entry[2]])
            # failed to find ID
            except AttributeError:
                pass

    return outputs


def filter_one_transcript_per_gene(transcript_list):
    """
    If there is more than one transcript per gene, keep only the longest.
    """

    print("Filtering to one transcript per gene...")

    # get a list of genes with the associated transcripts
    gene_list = collections.defaultdict(lambda: [])
    for transcript in transcript_list:
        transcript_info = transcript_list[transcript][0]
        gene_list[transcript_info[6]].append(transcript_info)
    # now get the longest one
    longest_transcripts = {}
    for gene_id in gene_list:
        if len(gene_list[gene_id]) > 1:
            lengths = [int(i[2]) - int(i[1]) for i in gene_list[gene_id]]
            longest_transcripts[gene_list[gene_id][lengths.index(max(lengths))][3]] = gene_list[gene_id][lengths.index(max(lengths))]
        else:
            longest_transcripts[gene_list[gene_id][0][3]] = gene_list[gene_id][0]
    return longest_transcripts


def list_transcript_ids_from_features(gtf_file_path, exclude_pseudogenes=True, full_chr=False):
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
        if entry[chr_index] in "123456789XY" and entry[chr_index+1] in "0123456789XY\t":
            # search for the transcript id
            id = re.search(transcript_id_pattern, entry)
            # if the id exists, add to list
            if id:
                ids.append(id.group(0))

    # now get a unique set
    unique_ids = sorted(list(set(ids)))

    return unique_ids


def write_features_to_bed(feature_list, output_file, clean_run = None):
    """
    Write a set of features to bed file

    Args:
        feature_list (dict): dictionary containing features with transcripts as keys
        output_file (str): path to output file
        clean_run (bool): if true, always write file
    """

    if not os.path.isfile(output_file) or clean_run:
        with open(output_file, "w") as outfile:
            for feature in feature_list:
                for item in feature_list[feature]:
                    item[3] = "{0}.{1}".format(item[3], item[4])
                    item[4] = "."
                    outfile.write("{0}\n".format("\t".join(gen.stringify(item))))
