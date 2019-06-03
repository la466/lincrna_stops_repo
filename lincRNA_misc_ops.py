import generic as gen
import conservation as cons
import file_ops as fo
import time
import os
import collections
import re

from Bio import SeqIO


def extract_bed_coordinates_block_format(input_bed, output_exons_bed, output_introns_bed):
    """
    Extract exons and introns from bed file containing both, with coordinates
    in block format

    Args:
        input_bed (str): path to input fasta file
        output_exons_bed (str): path to output exons file
        output_intons_bed (str): path to output introns file
    """

    # set up dictionary to hold coordinates
    exon_list = collections.defaultdict(lambda: collections.defaultdict())
    intron_list = collections.defaultdict(lambda: collections.defaultdict())
    # read in data
    data = gen.read_many_fields(input_bed, "\t")

    with open(output_exons_bed, "w") as output_exons:
        with open(output_introns_bed, "w") as output_introns:
            for line in data:
                start = int(line[1])
                id = line[3]
                strand = line[5]
                block_sizes = [int(i) for i in line[10].split(",") if len(i)]
                start_indices = [int(i) for i in line[11].split(",") if len(i)]
                # if on the reverse strand, need to reverse order
                if strand == "-":
                    block_sizes = block_sizes[::-1]
                    start_indices = start_indices[::-1]
                # now get a list of exon ids to use for intron calculations
                exon_ids = list(range(len(start_indices)))

                for i in range(len(start_indices)):
                    # now get the start and end of the exon coordinates
                    start_index = start + start_indices[i]
                    end_index = start_index + block_sizes[i]
                    # get the exon id
                    exon_id = i+1
                    # now write to the exons file
                    output_exons.write("{0}\t{1}\t{2}\t{3}.{4}\t.\t{5}\n".format(line[0], start_index, end_index, id, exon_id, strand))

                    if i+1 in exon_ids:
                        intron_id = "{0}-{1}".format(i+1, i+2)
                        if strand == "-":
                            intron_start = start + start_indices[i+1] + block_sizes[i+1]
                            intron_end = start_index
                        else:
                            intron_start = end_index
                            intron_end = start + start_indices[i+1]
                        output_introns.write("{0}\t{1}\t{2}\t{3}.{4}\t.\t{5}\n".format(line[0], intron_start, intron_end, id, intron_id, strand))

def build_transcripts(input_fasta, output_fasta):
    """
    Build transcripts from and exons file

    Args:
        input_fasta (str): path to input fasta file
        output_fasta (str): path to output fasta file
    """

    # read the input file
    names, seqs = gen.read_fasta(input_fasta)
    exon_list = collections.defaultdict(lambda: collections.defaultdict())
    for i, name in enumerate(names):
        exon_list[name.split(".")[0]][int(name.split(".")[-1].split("(")[0])] = seqs[i]
    with open(output_fasta, "w") as outfile:
        for id in exon_list:
            transcript_sequence = "".join([exon_list[id][i] for i in sorted(exon_list[id])])
            outfile.write(">{0}\n{1}\n".format(id, transcript_sequence))

def sort_by_exon_number(input_bed, single_exons_bed, multi_exons_bed):
    """
    Given a bed file containing exon coordinates, sort into those that
    only have one exon and those that have multiple.

    Args:
        input_bed (str): path to input bed file
        single_exons_bed (str): path to single exons output file
        multi_exons_bed (str): path to multi exons output file
    """
    # read in data
    data = gen.read_many_fields(input_bed, "\t")
    exon_list = collections.defaultdict(lambda: [])
    [exon_list[line[3].split(".")[0]].append(line) for line in data]
    # write to file
    with open(single_exons_bed, "w") as single_output:
        with open(multi_exons_bed, "w") as mutli_output:
            for id in exon_list:
                if len(exon_list[id]) > 1:
                    [mutli_output.write("{0}\n".format("\t".join(i))) for i in exon_list[id]]
                else:
                    [single_output.write("{0}\n".format("\t".join(i))) for i in exon_list[id]]

def sort_fasta_by_bed(input_bed, input_fasta, output_fasta):
    """
    Given a bed file, keep only the sequences that correspond to an entry in
    the bed file

    Args:
        input_bed (str): path to input bed file
        input_fasta (str): path to input fasta file
        output_fasta (str): path to output fasta file
    """

    # get the bed data
    bed_ids = list(set([i[3].split(".")[0] for i in gen.read_many_fields(input_bed, "\t")]))
    # get fasta entries
    names, seqs = gen.read_fasta(input_fasta)
    with open(output_fasta, "w") as outfile:
        for i, name in enumerate(names):
            if name.split(".")[0] in bed_ids:
                outfile.write(">{0}\n{1}\n".format(name, seqs[i]))


def extract_lncrna_only(input_file, output_file):
    """
    Given the file, take the id information and keep only those entries
    that correspond to lncRNA. Write to file
    """

    ids = []
    for entry in entries:
        type = re.findall("^ENSG\d+\.\d+:(.+)", entry[3])
        # if the type exists
        if len(type) != 0:
            splits = type[0].split(",")
            # and if there is only 1 entry
            if len(splits) == 1:
                # and that entry is lncRNA
                if splits[0] == "lncRNA":
                    ids.append(entry[1])
    with open(output_file, "w") as outfile:
        outfile.write("{0}\n".format("\t".join(sorted(ids))))


def extract_second_seqs(input_bed, input_file, genome_fasta, output_dir):
    """
    Extract the second set of sequences
    """
    # get a set of ids that correspond only to lincrna entries
    id_file = "{0}/lncrna_ids.txt".format(output_dir)
    extract_lncrna_only(input_file, id_file)

    # now keep only the bed entries that are in the id list
    filtered_bed = "{0}.filtered".format(input_bed)
    ids = gen.read_many_fields(id_file, "\t")
    bed_entries = gen.read_many_fields(input_bed, "\t")
    with open(filtered_bed, "w") as outfile:
        for entry in bed_entries:
            if entry[3] in ids:
                outfile.write("{0}\n".format("\t".join(entry)))

    # now write the bed to an exon bed
    exons_bed = "{0}.exons.bed".format(input_bed)
    fo.entries_to_bed(filtered_bed, exons_bed, hg38 = True)
    # now get the exon sequences
    exons_fasta = "{0}.exons.fasta".format(input_bed)
    fo.fasta_from_intervals(exons_bed, exons_fasta, genome_fasta, force_strand = True, names = True)

    # now generate the full transcript for multi exon transcripts
    transcripts_fasta = "{0}.multi_exon_transcripts.fasta".format(input_bed)
    names, seqs = gen.read_fasta(exons_fasta)
    seq_list = collections.defaultdict(lambda: collections.defaultdict())
    for i, name in enumerate(names):
        id = ".".join(name.split("(")[0].split(".")[:-1])
        exon = int(name.split("(")[0].split(".")[-1])
        seq_list[id][exon] = seqs[i]
    with open(transcripts_fasta, "w") as outfile:
        for id in sorted(seq_list):
            if len(seq_list[id]) > 1:
                exon_list = []
                for exon in sorted(seq_list[id]):
                    exon_list.append(seq_list[id][exon])
                seq = "".join(exon_list)
                if "N" not in seq and len(seq) >= 200:
                    # convert names to : here as otherwise it will run sorting later
                    id = ":".join(id.split("."))
                    outfile.write(">{0}\n{1}\n".format(id, seq))

    # blast to get paralogous families
    blast_db_path = "{0}/bast_db".format(output_directory)
    output_blast_file = "{0}/blast_output.csv".format(output_directory)
    families_file = "{0/families.txt".format(output_directory)
    gen.create_output_directories(blast_db_path)
    cons.filter_families(transcripts_fasta, output_blast_file, families_file, database_path = blast_db_path, clean_run = True)
