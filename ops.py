import generic as gen
import file_ops as fo
import numpy as np
import itertools as it
import os
import re
import collections
import random

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


def get_exon_coding(exon_bed, cds_bed, non_coding_bed_output):
    """
    Check whether exons are coding or not.

    Args:
        exon_bed (str): path to the full exon bed file
        cds_bed (str): path to the more limited bed file
        non_coding_bed_output (str): path to the non coding bed output file
    """

    # get the non coding exons
    remove_bed_intersects(exon_bed, cds_bed, non_coding_bed_output)


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
    gene_regex = re.compile("(?<=gene_id \")ENSG[0-9]*")
    #extract transcript IDs and feature numbers
    for pos, feature in enumerate(gtf_features):
        to_parse = feature[3]
        trans = re.search(trans_regex, to_parse).group(0)
        exon_no = re.search(exon_no_regex, to_parse).group(0)
        gene_id = re.search(gene_regex, to_parse).group(0)
        gtf_features[pos][3] = "{0}.{1}".format(trans, exon_no)
        gtf_features[pos].append(gene_id)
    #write to bed
    with open(bed, "w") as file:
        for feature in gtf_features:
            file.write("{0}\n".format("\t".join([str(i) for i in feature])))


def get_unique_transcripts(input_bed, input_fasta, output_fasta):
    """
    Generate a file containing transcripts where there is only one
    transcript per gene remaining.

    Args:
        input_bed (str): path to bed file containing transcript features
        input_fasta (str): path to fasta file containing transcript sequences
        output_fasta (str): path to output fasta file containing unique transcripts
    """

    # get the gene and transcript links
    gene_transcript_links = link_transcripts_to_genes(input_bed)
    # get the unique transcripts
    uniquify_transcripts(input_fasta, gene_transcript_links, output_fasta)


def filter_bed_from_fasta(input_bed, input_fasta, output_file):
    """
    Given a fasta file, filter the bed file to only include entries that
    also appear in the fasta file.

    Args:
        input_bed: path to the input bed file
        input_fasta: path to the input fasta file
        output_file: path to the output file
    """

    bed_lines = gen.read_many_fields(input_bed, "\t")
    names, seqs = gen.read_fasta(input_fasta)
    with open(output_file, "w") as outfile:
        for line in bed_lines:
            transcript = line[3].split(".")[0]
            if transcript in names:
                outfile.write("{0}\n".format("\t".join(line)))


def filter_bed_transcript_id(input_bed1, input_bed2, output_file):
    """
    Filter bed file 1 based on the trancript ids found in bed file 2.

    Args:
        input_bed1 (str): path to the first bed file to be filtered
        input_bed2 (str): path to the second bed file that will be used for filtering
        output_file (str): path to the output file
    """

    bed_lines1 = gen.read_many_fields(input_bed1, "\t")
    bed_lines2 = gen.read_many_fields(input_bed2, "\t")
    # get the transcript ids from the second bed file
    transcript_ids = [line[3].split('.')[0] for line in bed_lines2]
    with open(output_file, "w") as outfile:
        [outfile.write("{0}\n".format("\t".join(line))) for line in bed_lines1 if line[3].split(".")[0] in transcript_ids]


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

    print("{0} seqeunces after filtering...".format(pass_count))


def intersect_bed(bed_file1, bed_file2, overlap = False, overlap_rec = False, write_both = False, sort = False, output_file = None,
                             force_strand = False, no_name_check = False, no_dups = True, chrom = None, intersect = False, hit_count = False, bed_path = None, intersect_bam=None,
                  write_zero = False, write_bed = False, subtract=None, return_non_overlaps=None):
    '''Use bedtools/bedops to intersect coordinates from two bed files.
    Return those lines in bed file 1 that overlap with intervals in bed file 2.
    OPTIONS
    output_file: write output to this file
    use_bedops: use bedops rather than bedtools. Certain options are only valid with one of the two, see below.
    overlap: minimum overlap required as a fraction of the intervals in bed file 1 (EX: 0.8 means that the
    overlap has to be at least 80% of the intervals in bed file 1).
    overlap_rec: require that the overlap as a fraction of the intervals in file 2 be at least as high as
    the threshold indicated in -f.
    write_both: if True, return not only the interval from bed file 1 but, tagged onto the end, also the
    interval from bed file 2 that it overlaps (only
    valid when using bedtools).
    sort: sort bed files before taking the intersection
    force_strand: check that the feature and the bed interval are on the same strand (only valid with bedtools)
    no_name_check: if set to False, checks whether the chromosome names are the same in the too bed files (only valid with bedtools)
    no_dups: if True, only returns each interval once. If set to false, intervals in bed file 1 that overlap several intervals in
    bed file 2 will be returned several times (as many times as there are overlaps with different elements in bed file 2)
    chrom: limit search to a specific chromosome (only valid with bedops, can help in terms of efficiency)
    intersect: rather than returning the entire interval, only return the part of the interval that overlaps an interval in bed file 2.
    hit_count: for each element in bed file 1, return the number of elements it overlaps in bed file 2 (only valid with bedtools)
    intersect_bam: intersect a bam file with a bed file. Requires bam file to be called first
    write_zero: like write_both but also write A intervals that don't overlap with any B intervals,
    write_bed: when intersecting a bam file, write output as bed.'''
    gen.create_directory("temp_data/")
    temp_file_name = "temp_data/temp_bed_file{0}.bed".format(random.random())
    #have it write the output to a temporary file
    bedtools_output = run_bedtools(bed_file1, bed_file2, force_strand, write_both, chrom, overlap, sort, no_name_check, no_dups, output_file = temp_file_name, intersect = intersect, hit_number = hit_count, bed_path = bed_path, intersect_bam = intersect_bam, write_zero = write_zero, overlap_rec = overlap_rec, write_bed = write_bed, subtract = subtract, return_non_overlaps = return_non_overlaps)
    #move it to a permanent location only if you want to keep it
    if output_file:
        gen.run_process(["mv", temp_file_name, output_file])
    else:
        bedtools_output = gen.read_many_fields(temp_file_name, "\t")
    gen.remove_file(temp_file_name)
    return(bedtools_output)


def link_transcripts_to_genes(bed_file):
    """
    Take a bed file and link the transcript ids to gene ids

    Args:
        bed_file (str): path to bed file containing information

    Returns:
        transcript_gene_links (dict): dict[gene_id] = [transcript_id]
    """

    lines = gen.read_many_fields(bed_file, "\t")
    data = [[i[3].split(".")[0], i[-1]] for i in lines if len(lines) > 2]
    transcript_gene_links = collections.defaultdict(lambda: [])
    for transcript in data:
        # only include the transcript if its not already been included
        if transcript[0] not in transcript_gene_links[transcript[1]]:
            transcript_gene_links[transcript[1]].append(transcript[0])

    return transcript_gene_links


def remove_bed_intersects(bed_file1, bed_file2, output_file):
    intersect_bed(bed_file1, bed_file2, overlap_rec = True, output_file = output_file, force_strand = True, return_non_overlaps = True, no_dups=False)


def remove_bed_overlaps(input_file, output_file):
    '''
    Given a bed file, only leave non-overlapping elements, regardless of the strand of the overlap.
    Adapted from function written by RS

    Args:
        input_file (str): path to the input file to remove overlaps
        output_file (str): path to the output file
    '''
    #check how many columns there are in the bedfile
    with open(input_file) as file:
        line = file.readline()
        column_number = line.count("\t") + 1
    # merge overlapping intervals and have it count how many of the elements from the original file contribute to each
    # interval in the new file
    # note that bedops takes the column numbers in base 1
    if column_number > 3:
        columns = ",".join([str(i) for i in range(4, column_number + 1)] + ["1"])
        operations = ",".join(["distinct" for i in range(4, column_number + 1)] + ["count"])
    else:
        columns = "1"
        operations = "count"
    merge_result = gen.run_process(["bedtools", "merge", "-i", input_file, "-c", columns, "-o", operations])
    #only leave those intervals that do not result from a merge and delete counts column
    merge_result = merge_result.split("\n")
    with open(output_file, "w") as outfile:
        for line in merge_result:
            if line[-2:] == "\t1":
                outfile.write("{0}\n".format(line[:-2]))


def run_bedtools(A_file, B_file, force_strand = False, write_both = False, chrom = None, overlap = None, sort = False, no_name_check = False, no_dups = True, hit_number = False, output_file = None, intersect = False, bed_path = None, overlap_rec = None, intersect_bam = None, write_zero = None, write_bed = False, subtract=False, return_non_overlaps=None):
    '''
    See intersect_bed for details.
    '''
    if write_both:
        write_option = "-wo"
    elif hit_number:
        write_option = "-c"
    elif write_zero:
        write_option = "-wao"
    else:
        write_option = "-wa"
    if sort:
        sort_bed(A_file, A_file)
        sort_bed(B_file, B_file)
    if subtract:
        arg = "subtractBed"
    else:
        arg = "intersectBed"
    bedtools_args = [arg, "-a", A_file,"-b", B_file, write_option]
    if intersect:
        del bedtools_args[-1]
    if overlap:
        bedtools_args.extend(["-f", str(overlap)])
    if overlap_rec:
        bedtools_args.append("-r")
    if force_strand:
        bedtools_args.append("-s")
    if no_name_check:
        bedtools_args.append("-nonamecheck")
    if no_dups:
        bedtools_args.append("-u")
    if return_non_overlaps:
        bedtools_args.append("-v")
    if chrom:
        print("Bedtools cannot be restricted to a single chromosome. Use bedops!")
        raise Exception
    if hit_number and no_dups:
        print("When counting hits, each interval in the first bed file is only reported once by default. Set no_dups to False!")
        raise(Exception)
    if bed_path:
        bedtools_args[0] = "{0}{1}".format(bed_path, bedtools_args[0])
    if intersect_bam:
        if A_file[-4:] != ".bam":
            print("Bam file must be called first")
            raise Exception
        if B_file[-4:] != ".bed":
            print("Bed file must be called second")
            raise Exception
        bedtools_args = ["intersectBed", write_option, "-abam", A_file, "-b", B_file]
        if write_bed:
            bedtools_args.append("-bed")
    bedtools_output = gen.run_process(bedtools_args, file_for_output = output_file)
    return(bedtools_output)


def sort_bed(input_file, output_file):
    """
    Sort a bed file.

    Args:
        input_file (str): path to the input file
        output_file (str): path to the output file
    """

    # Do like this so we can sort a file and keep the same name
    gen.create_output_directories("temp_data")
    temp_file_name = "temp_data/temp_sorted_bed{0}.bed".format(random.random())
    gen.run_process(["sort-bed", input_file], file_for_output = temp_file_name)
    gen.run_process(["mv", temp_file_name, output_file])
    gen.remove_file(temp_file_name)


def uniquify_transcripts(input_fasta, transcript_gene_links, output_fasta):
    """
    Given a fasta and the links between genes and transcripts, filter to only
    leave one transcript per gene. Choose the longest transcript.

    Args:
        input_fasta (str): path to the input fasta
        transcript_gene_links (dict): the dictionary containing trancript links to genes, keys are gene ids
        output_fasta (str): path to the output fasta
    """

    names, seqs = gen.read_fasta(input_fasta)

    with open(output_fasta, "w") as outfile:
        for gene_id in sorted(transcript_gene_links):
            # get the lengths of the transcripts associated with the gene
            sequence_lengths = [len(seqs[names.index(transcript_id)]) for transcript_id in transcript_gene_links[gene_id]]
            # get the transcript id of the longest transcript
            max_length_transcript = transcript_gene_links[gene_id][sequence_lengths.index(max(sequence_lengths))]
            # write the longest transcript to file
            outfile.write(">{0}\n{1}\n".format(max_length_transcript, seqs[names.index(max_length_transcript)]))
