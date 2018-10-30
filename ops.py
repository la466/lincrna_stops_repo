import generic as gen
import seq_ops as seqo
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

def calc_motif_gc(motif_set):
    motif = list("".join(motif_set))
    return np.divide(sum([motif.count(i) for i in ["G", "C"]]), len(motif))

def calc_motif_density(seq, motif_set):
    codon_regex = re.compile(".{3}")
    c1 = re.findall(codon_regex, seq[1:])
    c2 = re.findall(codon_regex, seq[2:])
    return (np.sum([c1.count(codon) for codon in motif_set]) + np.sum([c2.count(codon) for codon in motif_set]))*3


def calc_motif_densities(motif_sets, cds_fasta, temp_dir):

    temp_files = []
    cds_seqs = gen.read_fasta(cds_fasta)[1]

    if motif_sets:
        for i, set in enumerate(motif_sets):
            print("{0}/{1} {2}".format(i+1, len(motif_sets), "_".join(set)))
            temp_file = "{0}/{1}.bed".format(temp_dir, "_".join(set))
            temp_files.append(temp_file)
            if not os.path.isfile(temp_file):
                stop_count = 0
                nt_count = 0
                for seq in cds_seqs:
                    stop_count += calc_motif_density(seq[:-3], set)
                    nt_count += len(seq[:-3])
                with open(temp_file, "w") as outfile:
                    outfile.write("{0},{1}\n".format(calc_motif_gc(set), np.divide(stop_count, nt_count)))
    return temp_files



def check_exon_files(input_bed1, input_bed2):
    """
    Do a sanity check to make sure there are no coding exons in the
    non coding exons file and vice versa.

    Args:
        input_bed1 (str): path to the first bed file
        input_bed2 (str): path to the second bed file

    Returns:

    """

    bed_lines1 = gen.read_many_fields(input_bed1, "\t")
    bed_lines2 = gen.read_many_fields(input_bed2, "\t")
    transcripts1 = [line[3] for line in bed_lines1]
    transcripts2 = [line[3] for line in bed_lines2]
    # get any overlap
    overlap = list(set(transcripts1) & set(transcripts2))
    if len(overlap):
        print("Something's gone wrong. Coding exons and non coding exons are present in both files...")
        raise Exception
    return True


def get_coding_exons(exons_file, cds_file, output_file, remove_overlapping = False):
    """
    Given a bed file of exon coordinates and a bed file of CDS coordinates,
    write a new bed file that only contains those exon coordinates form the former file that are
    1) fully coding 2) internal
    NB! Assumes that all the coordinates are from non-overlapping transcripts. If this is not the case,
    set remove_overlaps to True and it'll remove overlapping intervals.
    Modified from LA and RS.

    Args:
        exons_file (str): path to the bed file containing exon coordinates
        cds_file (str): path to bed file containing the cds coordinates
        output_file (str): path to output file
        remove_overlapping (bool): if true, remove overlapping intervals
    """

    if remove_overlapping:
        sort_bed(exons_file, exons_file)
        remove_bed_overlaps(exons_file, exons_file)
    #filter out anything that isn't fully coding
    #you have to write_both because you want to make sure that they
    #haven't been kept because of an overlap to a transcript that doesn't appear in the exons file
    gen.create_output_directories("temp_data")
    temp_file = "temp_data/temp{0}.txt".format(random.random())
    intersect_bed(exons_file, cds_file, overlap = 1, overlap_rec = True, output_file = temp_file, force_strand = True, write_both = True, no_dups = False, no_name_check = False)
    #filter out terminal exons
    #in theory, there shouldn't be any left after the previous step
    #in practice, there may be unannotated UTRs, so it looks like we have a fully coding terminal exon,
    #whereas in reality, the exon is only partially coding
    temp_file2 = "temp_data/temp{0}.txt".format(random.random())
    with open(temp_file2, "w") as outfile:
        #figure out the rank of the last exon for each transcript
        filt_exons = gen.read_many_fields(exons_file, "\t")
        filt_exons = [i for i in filt_exons if len(i) > 3]
        names = [i[3].split(".") for i in filt_exons]
        names = gen.list_to_dict(names, 0, 1, as_list = True)
        names = {i: max([int(j) for j in names[i]]) for i in names}
        coding_exons = gen.read_many_fields(temp_file, "\t")
        for exon in coding_exons:
            overlap_name = exon[10].split(".")
            if overlap_name[0] in names:
                name = exon[3].split(".")
                if name[-1] != "1":
                    last_exon = names[name[0]]
                    if int(name[-1]) != last_exon:
                        exon = [str(i) for i in exon[:7]]
                        outfile.write("{0}\n".format("\t".join(exon)))
    sort_bed(temp_file2, temp_file2)
    gen.run_process(["mergeBed", "-i", temp_file2, "-c", "4,5,6,7", "-o", "distinct,distinct,distinct,distinct"], file_for_output = output_file)
    gen.remove_file(temp_file)
    gen.remove_file(temp_file2)


def get_exon_coding(exons_bed, quality_filtered_bed, final_cds_bed, non_coding_bed_output, coding_bed_output):
    """
    Check whether exons are coding or not.

    Args:
        exons_bed (str): path to the full exon bed file
        cds_bed (str): path to the more limited bed file
        non_coding_bed_output (str): path to the non coding bed output file
    """

    # get the non coding exons
    remove_bed_intersects(exons_bed, final_cds_bed, non_coding_bed_output)
    # get coding exons
    get_coding_exons(exons_bed, final_cds_bed, coding_bed_output, remove_overlapping=True)
    # do a sanity check and ensure no coding exons are in non coding file and vice versa
    check_exon_files(coding_bed_output, non_coding_bed_output)


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


def extract_gtf_features_all(gtf, bed, exclude_XY=None):
    """
    Get all the features in a gtf file in bed format.

    Args:
        gtf (str): path to gtf file
        bed (str): path to output bed file
        exclude_XY (bool): if true, do not include X and Y chromosomes
    """

    # get all feature lines
    lines = [line for line in gen.read_many_fields(gtf, "\t") if "#" not in line[0]]
    # format as .bed and switch to base 0
    gtf_features = [["chr{0}".format(i[0]), int(i[3]) - 1, i[4], ".", i[2], i[6], i[8]] for i in lines if len(i) >= 3]
    trans_regex = re.compile("(?<=transcript_id \")ENST[0-9]*")
    exon_no_regex = re.compile("(?<=exon_number \")[0-9]*")
    gene_regex = re.compile("(?<=gene_id \")ENSG[0-9]*")
    biotype_regex = re.compile("(?<=gene_biotype \").*?(?=\";)")
    for pos, feature in enumerate(gtf_features):
        to_parse = feature[-1]
        trans = re.search(trans_regex, to_parse)
        exon_no = re.search(exon_no_regex, to_parse)
        gene_id = re.search(gene_regex, to_parse)
        biotype = re.search(biotype_regex, to_parse)
        if not trans:
            trans = "nan"
        else:
            trans = trans.group(0)
        if not exon_no:
            exon_no = "nan"
        else:
            exon_no = exon_no.group(0)
        if not gene_id:
            gene_id = "nan"
        else:
            gene_id = gene_id.group(0)
        if not biotype:
            biotype = "nan"
        else:
            biotype = biotype.group(0)
        gtf_features[pos][3] = "{0}.{1}.{2}".format(trans, exon_no, gene_id)
        gtf_features[pos] = gtf_features[pos][:-1]
        gtf_features[pos].append(biotype)

    # exclude X and Y chrs if required
    required_chrs = ["chr{0}".format(i) for i in range(1, 24)]
    if not exclude_XY:
        required_chrs.extend(["chrX", "chrY"])

    #write to bed
    with open(bed, "w") as outfile:
        for feature in gtf_features:
            if feature[0] in required_chrs:
                outfile.write("{0}\n".format("\t".join([str(i) for i in feature])))



def get_exon_flank_reading_frame(coding_exons_fasta, full_exons_fasta, output_file):
    """
    Get the reading frames that the flanks of an exon starts in. If the first position of the codon
    this is 0, second is 1 and last is 2.

    Args:
        coding_exons_fasta (str): path to file containing coding exon sequences
        full_exons_fasta (str): path to file containing all exon sequences
        output_file (str): path to the output file
    """

    # read in the coding exons and all exons
    full_names, full_seqs = gen.read_fasta(full_exons_fasta)
    full_coding_names = gen.read_fasta(coding_exons_fasta)[0]

    # create a dictionary that hold each of the full sequences
    full_exon_seqs = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict()))
    for i, name in enumerate(full_names):
        id = name.split(".")[0]
        exon_id = int(name.split(".")[1].split("(")[0])
        full_exon_seqs[id][exon_id] = full_seqs[i]

    # now get the reading frame that each flank starts in
    flanks_reading_frames = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: [])))
    for id in full_exon_seqs:
        seq_length = 0
        for exon_id, seq in sorted(full_exon_seqs[id].items()):
            # sequence needs to be longer than 138 nts to have two flanks
            if len(seq) > 138:
                five_prime_start = seq_length + 2
                three_prime_start = seq_length + len(seq) - 69
                flanks_reading_frames[id][exon_id] = [five_prime_start % 3, three_prime_start % 3]

            seq_length += len(seq)

    # now output the reading frames for just the coding exons
    with open(output_file, "w") as outfile:
        for name in full_coding_names:
            id = name.split(".")[0]
            exon_id = int(name.split(".")[1].split("(")[0])
            if exon_id in flanks_reading_frames[id]:
                outfile.write(">{0}\n{1},{2}\n".format(name, flanks_reading_frames[id][exon_id][0], flanks_reading_frames[id][exon_id][1]))


def get_genome_bed_from_fasta_index(features_bed, fasta_index, output_file):
    """
    Given a list of features, get the genome coordinates as a bed file.

    Args:
        features_bed (str): path to bed file containing features
        fasta_index (str): path to fasta index file
        output_file (str): path to output file
    """

    # get all the chromosomes required
    first_column = [i for i in list(set(gen.run_process(["awk", "{print $1}"], file_for_input=features_bed).split('\n'))) if len(i)]
    # get index lines
    index = gen.read_many_fields(fasta_index, "\t")
    with open(output_file, "w") as outfile:
        for i in index:
            if i[0] in first_column:
                start = 0
                length = int(i[1])
                out_info = [i[0], start, start+length, ".", "."]
                outfile.write("{0}\t+\n{0}\t-\n".format("\t".join(gen.stringify(out_info))))


def get_passed_NONCODE_codes(input_fasta, codes_file, mapping_file, output_fasta, code):
    """
    Only keep sequences that have particular NONCODE code

    Args:
        input_fasta (str): path to input fasta file
        codes_file (str): path to file containing code
        mapping_file (str): path to transcript-gene mapping file
        output_fasta (str): path to output fasta
        code (str): code to look for. As string because cant pass 0001 through
    """

    codes = {code[0]: code[1] for code in gen.read_many_fields(codes_file, "\t")}
    mappings = {name[0].split(".")[0]: name[1] for name in gen.read_many_fields(mapping_file, "\t")}

    names, seqs = gen.read_fasta(input_fasta)
    with open(output_fasta, "w") as outfile:
        for i,name in enumerate(names):
            gene = mappings[name]
            seq_code = codes[gene]
            if seq_code == code:
                outfile.write(">{0}\n{1}\n".format(name, seqs[i]))


def get_exon_reading_frame(coding_exons_fasta, full_exons_fasta, output_file):
    """
    Get the reading frame an exon starts in. If the first position of the codon
    this is 0, second is 1 and last is 2.

    Args:
        coding_exons_fasta (str): path to file containing coding exon sequences
        full_exons_fasta (str): path to file containing all exon sequences
        output_file (str): path to the output file
    """

    # read in the coding exons and all exons
    full_names, full_seqs = gen.read_fasta(full_exons_fasta)
    full_coding_names = gen.read_fasta(coding_exons_fasta)[0]

    # create a dictionary that hold each of the full sequences
    full_exon_seqs = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict()))
    for i, name in enumerate(full_names):
        id = name.split(".")[0]
        exon_id = int(name.split(".")[1].split("(")[0])
        full_exon_seqs[id][exon_id] = full_seqs[i]

    # now get the reading frame that they start in
    full_reading_frames = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict()))
    for id in full_exon_seqs:
        seq_length = 0
        for exon_id, seq in sorted(full_exon_seqs[id].items()):
            full_reading_frames[id][exon_id] = seq_length % 3
            seq_length += len(seq)

    # now output the reading frames for just the coding exons
    with open(output_file, "w") as outfile:
        for name in full_coding_names:
            id = name.split(".")[0]
            exon_id = int(name.split(".")[1].split("(")[0])
            # print(name, full_reading_frames[id][exon_id])
            outfile.write(">{0}\n{1}\n".format(name, full_reading_frames[id][exon_id]))


def get_introns_from_bed(input_bed, output_file):

    exons = gen.read_many_fields(input_bed, "\t")

    exon_list = collections.defaultdict(lambda: collections.defaultdict(lambda: []))
    strands = {}
    chrs = {}

    for exon in exons:
        transcript = exon[3].split(".")[0]
        exon_id = exon[3].split(".")[1]
        if "(" in exon_id:
            exon_id = exon_id.split("(")[0]
        exon_id = int(exon_id)
        start = int(exon[1])
        end = int(exon[2])
        chr = exon[0]
        strand = exon[5]
        strands[transcript] = strand
        chrs[transcript] = chr
        exon_list[transcript][exon_id] = [start, end]

    with open(output_file, "w") as outfile:
        for transcript in exon_list:
            for exon_id in sorted(exon_list[transcript]):
                # check whether the next exon exists
                if exon_id + 1 in exon_list[transcript]:
                    if strands[transcript] == "+":
                        intron_start = exon_list[transcript][exon_id][1]
                        intron_end = exon_list[transcript][exon_id+1][0]
                    else:
                        intron_start = exon_list[transcript][exon_id+1][1]
                        intron_end = exon_list[transcript][exon_id][0]
                    outdata = [chrs[transcript], intron_start, intron_end, "{0}.{1}-{2}".format(transcript, exon_id, exon_id+1), ".", strands[transcript]]
                    outfile.write("{0}\n".format("\t".join(gen.stringify(outdata))))



def get_region_stop_counts(input_list, window_start, window_end):
    """
    Get the number of stop coding in the defined regions for a dictionary
    containing sequences

    Args:
        input_list (dict): dictionary containing sequences dict[name] = sequence
        window_start (int): the nucleotide to start the flank region (base 1)
        window_end (int): the nucleotide to end the flank region (base 1)

    Returns:
        region_stop_counts (list): list containing [flank_count, core_count]
    """

    region_stop_counts = [0,0]
    for name, seq in input_list.items():
        five_prime_start = window_start - 1
        five_prime_end = window_end
        three_prime_start = len(seq) - window_end
        three_prime_end = len(seq) - window_start + 1

        # get the region sequences
        five_prime_flank = seq[five_prime_start:five_prime_end]
        three_prime_flank = seq[three_prime_start:three_prime_end]
        core = seq[five_prime_end:three_prime_start]

        # get the stop codon counts
        region_stop_counts[0] += seqo.get_stop_count(five_prime_flank)
        region_stop_counts[0] += seqo.get_stop_count(three_prime_flank)
        region_stop_counts[1] += seqo.get_stop_count(core)

    return region_stop_counts


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

    print("{0} sequences after filtering...".format(pass_count))


def filter_seq_lengths(input_file, output_file, length):
    """
    Given a fasta filter, filter and keep only sequences longer than the given length

    Args:
        input_file (str): path to input fasta
        output_file (str): path to output fasta
        length (int): length threshold to retain sequences longer than
    """

    names, seqs = gen.read_fasta(input_file)
    with open(output_file, "w") as outfile:
        for i, name in enumerate(names):
            if len(seqs[i]) >= length:
                outfile.write(">{0}\n{1}\n".format(name, seqs[i]))


def intersect_bed(bed_file1, bed_file2, overlap = False, overlap_rec = False, write_both = False, sort = False, output_file = None,
                             force_strand = False, no_name_check = False, no_dups = True, intersect = False, hit_count = False, bed_path = None, intersect_bam=None,
                  write_zero = False, write_bed = False, subtract=None, return_non_overlaps=None, write_none=False):
    """
    Use bedtools to intersect coordinates from two bed files.
    Return those lines in bed file 1 that overlap with intervals in bed file 2.
    Adapted from RS.

    Args:
        bed_file1 (str): path to first bed file (could be bam file if intersect_bam=True)
        bed_file2 (str): path to second bed file
        overlap (float): minimum overlap required as a fraction of the intervals in bed file 1 (EX: 0.8 means that the overlap has to be at least 80% of the intervals in bed file 1).
        overlap_rec (bool): if true, require that the overlap as a fraction of the intervals in file 2 be at least as high as the threshold indicated in -f.
        write_both (bool): if true, return not only the interval from bed file 1 but, tagged onto the end, also the interval from bed file 2 that it overlaps.
        sort (bool): if true, sort bed files before taking the intersection
        output_file (str): if exists, path to with which to write output file
        force_strand (bool): if true, check that the feature and the bed interval are on the same strand
        no_name_check (bool): if false, checks whether the chromosome names are the same in the too bed files
        no_dups (bool): if true, only returns each interval once. If false, intervals in bed file 1 that overlap several intervals in bed file 2 will be returned several times (as many times as there are overlaps with different elements in bed file 2)
        intersect (bool): if true, rather than returning the entire interval, only return the part of the interval that overlaps an interval in bed file 2.
        hit_count (bool): for each element in bed file 1, return the number of elements it overlaps in bed file 2
        intersect_bam (bool): if true, intersect a bam file with a bed file. Requires bam file to be called first
        write_zero (bool): like write_both but also write A intervals that don't overlap with any B intervals
        write_bed (bool): if true, when intersecting a bam file, write output as bed
        subtract (bool): if true, set argument to subtractBed
        return_non_overlaps (bool): if true, only return entries in bed file 1 that don't overlap bed file 2

    Returns:
        bedtools_output (list): list of bed lines from the output file
    """

    gen.create_directory("temp_data/")
    temp_file_name = "temp_data/temp_bed_file{0}.bed".format(random.random())
    #have it write the output to a temporary file
    bedtools_output = run_bedtools(bed_file1, bed_file2, force_strand, write_both, overlap, sort, no_name_check, no_dups, output_file = temp_file_name, intersect = intersect, hit_number = hit_count, bed_path = bed_path, intersect_bam = intersect_bam, write_zero = write_zero, overlap_rec = overlap_rec, write_bed = write_bed, subtract = subtract, return_non_overlaps = return_non_overlaps, write_none = write_none)
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
    """
    Keep entries from one bed file that strictly don't overlap a second bed file.

    Args:
        bed_file1 (str): path to the first bed file
        bed_file2 (str): path to the second bed file
        output_file (str): path to the output file
    """

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


def run_bedtools(A_file, B_file, force_strand = False, write_both = False, overlap = None, sort = False, no_name_check = False, no_dups = True, hit_number = False, output_file = None, intersect = False, bed_path = None, overlap_rec = None, intersect_bam = None, write_zero = None, write_bed = False, subtract=False, return_non_overlaps=None, write_none = None):
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
    bedtools_args = [arg, "-a", A_file,"-b", B_file]
    if not write_none:
        bedtools_args.append(write_option)
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
    gen.run_process(["sortBed", "-i", input_file], file_for_output = temp_file_name)
    gen.run_process(["mv", temp_file_name, output_file])
    gen.remove_file(temp_file_name)


def uniquify_lincRNA_transcripts(input_fasta, mapping_file, output_fasta):
    """
    Given the set of lincRNA and the transcript-gene mapping, filter to only
    leave one per gene, keeping the longest.

    Args:
        input_fasta (str): path to fasta file containing the transcripts
        mapping_file (str): file containing the transcript gene mappings
        output_fasta (str): path to output fasta
    """

    # get the mappings
    mappings = {line[0].split(".")[0]: line[1] for line in gen.read_many_fields(mapping_file, "\t")}
    # get the sequences
    names, seq_list = gen.read_fasta(input_fasta)
    seqs = {name: seq_list[i] for i, name in enumerate(names)}

    mapped_names = collections.defaultdict(lambda: [])
    mapped_seq_lengths = collections.defaultdict(lambda: [])

    # add the mappings to dicts
    for transcript in seqs:
        if transcript in mappings:
            mapped_names[mappings[transcript]].append(transcript)
            mapped_seq_lengths[mappings[transcript]].append(len(seqs[transcript]))

    # get the longest transcript for each case
    max_lengths = {loc: mapped_names[loc][mapped_seq_lengths[loc].index(max(mapped_seq_lengths[loc]))] for loc in mapped_seq_lengths}

    # now write sequence to file
    with open(output_fasta, "w") as outfile:
        for loc in max_lengths:
            outfile.write(">{0}\n{1}\n".format(max_lengths[loc], seqs[max_lengths[loc]]))


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
