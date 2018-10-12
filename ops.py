import generic as gen
import file_ops as fo
import numpy as np
import itertools as it
import os
import re
import collections


def build_sequences(input_fasta, input_stops_fasta, output_fasta):
    """
    Build sequences from a provided fasta

    Args:
        input_fasta (str): path to input fasta
        input_stops_fasta (str): path to fasta containing the stop codons
        output_fasta (str): path for output fasta containing built sequences
    """

    names, seqs = gen.read_fasta(input_fasta)
    stop_names, stop_seqs = gen.read_fasta(input_stops_fasta)

    # create a dictionary containing all the parts
    sequence_parts_list = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict()))
    for i, name in enumerate(names):
        splits = name.split('(')
        name_splits = splits[0].split('.')
        transcript_id = name_splits[0]
        exon_id = int(name_splits[1])
        sequence_parts_list[transcript_id][exon_id] = seqs[i]

    # now build the sequences out
    sequence_dict = collections.defaultdict(lambda: [])
    for transcript_id in sequence_parts_list:
        for exon_id in sorted(sequence_parts_list[transcript_id]):
            sequence_dict[transcript_id].append(sequence_parts_list[transcript_id][exon_id])

    # add the stop codons
    for i, stop_name in enumerate(stop_names):
        transcript_id = stop_name.split(".")[0]
        sequence_dict[transcript_id].append(stop_seqs[i])

    #output to file
    with open(output_fasta, "w") as outfile:
        for transcript_id in sorted(sequence_dict):
            outfile.write(">{0}\n{1}\n".format(transcript_id, "".join(sequence_dict[transcript_id])))


def get_non_coding_exons(genome_fasta, gtf, coding_exons_bed, non_coding_exons_bed, coding_exons_fasta, non_coding_exons_fasta, output_directory):
    """
    Extract the coding and non coding exons.

    Args:
        genome_fasta (str): path to genome fasta
        gtf (str): path to genome gtf
        coding_exons_bed (str): path to coding exons bed output file
        non_coding_exons_bed (str): path to non coding exons bed output file
        coding_exons_fasta (str): path to coding exons fasta output file
        non_coding_exons_fasta (str): path to non coding exons fasta output file
        output_directory (str): path to the output directory
    """

    # get the required features
    full_exons_bed = "{0}/full_exons.bed".format(output_directory)
    full_CDS_bed = "{0}/full_CDS.bed".format(output_directory)
    full_stops_bed = "{0}/full_stop_codons.bed".format(output_directory)
    if not os.path.isfile(full_exons_bed) or not os.path.isfile(full_CDS_bed) or not os.path.isfile(full_stops_bed):
        print("Extracting features from .gtf file...")
        extract_gtf_features(gtf, ["exon"], full_exons_bed)
        extract_gtf_features(gtf, ["CDS"], full_CDS_bed)
        extract_gtf_features(gtf, ["stop_codon"], full_stops_bed)

    # get the fasta sequences of the features
    full_exons_fasta = "{0}/full_exons.fasta".format(output_directory)
    full_CDS_fasta = "{0}/full_CDS.fasta".format(output_directory)
    full_stops_fasta = "{0}/full_stops.fasta".format(output_directory)
    if not os.path.isfile(full_exons_fasta) or not os.path.isfile(full_CDS_fasta) or not os.path.isfile(full_stops_fasta):
        print("Converting bed to fasta...")
        fo.fasta_from_intervals(full_exons_bed, full_exons_fasta, genome_fasta, names=True)
        fo.fasta_from_intervals(full_CDS_bed, full_CDS_fasta, genome_fasta, names=True)
        fo.fasta_from_intervals(full_stops_bed, full_stops_fasta, genome_fasta, names=True)

    # build the coding sequences
    cds_sequences_fasta = "{0}/cds_sequences.fasta".format(output_directory)
    if not os.path.isfile(cds_sequences_fasta):
        build_sequences(full_CDS_fasta, full_stops_fasta, cds_sequences_fasta)

    # filter the coding sequences
    filtered_cds_sequences_fasta = "{0}/filtered_cds_sequences.fasta".format(output_directory)
    filter_coding_sequences(cds_sequences_fasta, filtered_cds_sequences_fasta)


def extract_features(gtf_file, features, output_file, full_chr_name=None, clean_chrom_only = False):
    """
    Given a GTF file, extract coordinates for specific features and write to .bed.

    Args:
        gtf_file (str): path to gtf file
        features (list): list of feature names to be extracted from the file
        output_file (str): path to the output file
        full_chr_name (bool): if true, add 'chr' to start or name entry
        clean_chrom_only (bool): if true, only allow chromosome names in [1,..,23,X,Y]
    """

    feature_list = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.UserList())))))
    if len(features) > 1:
        list_feature = True
    else:
        list_feature = False

    if clean_chrom_only:
        allowed = [str(i) for i in range(1, 23)] + ["X", "Y"]

    lines = gen.read_many_fields(gtf_file, "\t")
    #ensure protein coding and not entry meta
    lines = [line for line in lines if not line[0].startswith("#") and "pseudogene" not in line[-1]]
    #compile regex to find genes, transcripts and exons
    gene_regex = re.compile("(?<=gene_id \")ENSG[0-9]*")
    trans_regex = re.compile("(?<=transcript_id \")ENST[0-9]*")
    exon_no_regex = re.compile("(?<=exon_number \")[0-9]*")

    for line in lines:
        #if the feature identifier is in the list
        if line[2] in features:
            #get entry meta
            meta = line[-1]
            gene = re.search(gene_regex, meta).group(0)
            trans = re.search(trans_regex, meta).group(0)
            #I added the try ... except in case you want to extract, say, transcript
            #features where the exon number information wouldn't be present.
            try:
                exon_no = re.search(exon_no_regex, meta).group(0)
            except AttributeError:
                exon_no = 0
            chr_no = line[0]
            add = True
            if clean_chrom_only:
                if chr_no not in allowed:
                    add = False
            if add:
                feature_list[chr_no][gene][trans][exon_no][line[2]].append([line[3], line[4], line[6]])
    #output features sorted by feature, chr, gene, transcript id
    with open(output_file, 'w') as output:
        for chr_no in sorted(feature_list):
            for gene in sorted(feature_list[chr_no]):
                for trans in sorted(feature_list[chr_no][gene]):
                    for exon in sorted(feature_list[chr_no][gene][trans]):
                        for feature in sorted(feature_list[chr_no][gene][trans][exon]):
                            for item in feature_list[chr_no][gene][trans][exon][feature]:
                                if not list_feature:
                                    feature = feature
                                if full_chr_name:
                                    chr_name = "chr{0}".format(chr_no)
                                else:
                                    chr_name = str(chr_no)
                                #output and convert to base 0
                                output.write('\t'.join([chr_name, str(int(item[0])-1), item[1], '{0}.{1}.{2}'.format(trans, exon, gene), feature, item[2]]) + '\n')


def extract_gtf_features(gtf, features, bed):
    """
    Extract given features coordinates from a .gtf file and write to .bed

    Args:
        gtf (str): path to gtf file
        features (list): list of features to extract
        bed (str): path to output bed file
    """

    feature_list = []
    #iterate over the desired features
    for feature in features:
        #extract feature from GTF
        gtf_features = gen.run_process(["grep", "\t{0}\t".format(feature), gtf])
        #filter down to only protein-coding ones
        gtf_features = gen.run_process(["grep", "transcript_biotype \"protein_coding\""], input_to_pipe = gtf_features)
        #split lines
        gtf_features = [i.split("\t") for i in gtf_features.split("\n")]
        #add feature to list
        [i.append(feature) for i in gtf_features]
        #append to list
        feature_list.extend(gtf_features)
    #format as .bed. Switch to base 0.
    gtf_features = [["chr{0}".format(i[0]), int(i[3]) - 1, i[4], i[8], i[-1], i[6]] for i in feature_list if len(i) >= 3]

    #pre-compile regex
    trans_regex = re.compile("(?<=transcript_id \")ENST[0-9]*")
    exon_no_regex = re.compile("(?<=exon_number \")[0-9]*")
    #extract transcript IDs and feature numbers
    for pos, feature in enumerate(gtf_features):
        to_parse = feature[3]
        trans = re.search(trans_regex, to_parse).group(0)
        exon_no = re.search(exon_no_regex, to_parse).group(0)
        gtf_features[pos][3] = "{0}.{1}".format(trans, exon_no)
    #write to bed
    with open(bed, "w") as file:
        for feature in gtf_features:
            file.write("{0}\n".format("\t".join([str(i) for i in feature])))


def filter_coding_sequences(input_fasta, output_fasta):
    """
    Quality filter coding sequences

    Args:
        input_fasta (str): path to fasta file containing input CDS sequences
        output_fasta (str): path to output file to put sequences that pass filtering
    """

    # copile regex searches
    actg_regex = re.compile("[^ACTG]")
    codon_regex = re.compile(".{3}")

    stop_codons = ["TAA", "TAG", "TGA"]

    # read the sequences
    names, seqs = gen.read_fasta(input_fasta)

    print("{0} sequences prior to filtering".format(len(seqs)))

    s = []
    stop = []
    le = []
    non = []
    inf = []

    with open(output_fasta, "w") as outfile:
        pass_count = 0
        for i, name in enumerate(names):
            seq = seqs[i]
            passed = True
            while passed:
                # check to see if the first codon is ATG
                if seq[:3] != "ATG":
                    s.append(name)
                    passed = False
                # check to see if the last codon is a stop codon
                if seq[-3:] not in stop_codons:
                    stop.append(name)
                    passed = False
                # check to see if sequence is a length that is a
                # multiple of 3
                if len(seq) % 3 != 0:
                    le.append(name)
                    passed = False
                # check if there are any non ACTG characters in string
                non_actg = re.subn(actg_regex, '!', seq)[1]
                if non_actg != 0:
                    non.append(name)
                    passed = False
                # check if there are any in frame stop codons
                codons = re.findall(codon_regex, seq[3:-3])
                inframe_stops = [codon for codon in codons if codon in stop_codons]
                if len(inframe_stops):
                    inf.append(name)
                    passed = False

                outfile.write(">{0}\n{1}\n".format(name, seq))
                pass_count += 1

    print("{0} seqeunces after filtering...".format(pass_count))
