import generic as gen
import file_ops as fo
from regex_patterns import exon_number_pattern, gene_id_pattern, transcript_id_pattern
import re
import os
import collections


def generate_genome_dataset(gtf_file, genome_fasta, dataset_name, dataset_output_directory, id_list = None, filter_by_transcript = None, filter_by_gene = None, clean_run = None):

    # if we want a clean run, remove any previous output directory
    if clean_run:
        gen.remove_directory(dataset_output_directory)
    gen.create_output_directories(dataset_output_directory)

    # initalise genome object
    genome = Genome_Functions(gtf_file, genome_fasta, dataset_name, dataset_output_directory, clean_run)
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
    genome.filter_one_transcript_per_gene()
    # get a list of genes from the cleaned transcripts
    genome.list_genes_from_transcripts()

    return genome.filelist


class Genome_Functions(object):

    def __init__(self, gtf_file, genome_fasta, dataset_name, dataset_output_directory, clean_run):
        self.gtf_file = gtf_file
        self.genome_fasta = genome_fasta
        self.dataset_name = dataset_name
        self.dataset_output_directory = dataset_output_directory
        self.clean_run = clean_run
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
            generate_genome_features_dataset(self.dataset_name, self.gtf_file, dataset_features_bed, input_list = input_list, filter_by_transcript = filter_by_transcript, filter_by_gene = filter_by_gene)
        self.dataset_features_bed = dataset_features_bed
        self.filelist["dataset_features_bed"] = dataset_features_bed


    def list_genes_from_transcripts(self):
        """
        Get a list of transcripts with the associated genes
        """

        unique_transcript_gene_list_file = "{0}/{1}.cds.genes.bed".format(self.dataset_output_directory, self.dataset_name)
        if not os.path.isfile(unique_transcript_gene_list_file) or self.clean_run:
            get_clean_genes(self.unique_cds_fasta, self.full_cds_features_bed, unique_transcript_gene_list_file)
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
                self.ids = list(set([i[3] for i in entries]))

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
            quality_filter_cds_sequences(self.full_cds_fasta, quality_filtered_cds_fasta)
        self.quality_filtered_cds_fasta = quality_filtered_cds_fasta
        self.filelist["quality_filtered_cds_fasta"] = quality_filtered_cds_fasta


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


def extract_gtf_features(input_list, gtf_file_path, filter_by_transcript = None, filter_by_gene = None):
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
            transcript_search_results = re.search(transcript_id_pattern, entry_info)
            gene_search_results = re.search(gene_id_pattern, entry_info)
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


def filter_one_transcript_per_gene(transcript_list):
    """
    If there is more than one transcript per gene, keep only the longest.

    Args:
        transcript_list (dict): dict of transcripts

    Returns:
        longest_transcripts (dict): dict containing the longest transcript per gene
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


def generate_genome_features_dataset(dataset_name, dataset_gtf_file, dataset_features_output_bed, input_list = None, filter_by_transcript = None, filter_by_gene = None):
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
    args = [dataset_gtf_file, filter_by_transcript, filter_by_gene]
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


def quality_filter_cds_sequences(input_fasta, output_fasta):
    """
    Quality filter coding sequences

    Args:
        input_fasta (str): path to fasta file containing input CDS sequences
        output_fasta (str): path to output file to put sequences that pass filtering
    """

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
