import generic as gen
import ops
import collections
import os
import path

def build_seqs_from_exons_fasta(input_fasta, output_fasta):
    """
    Build a list of sequences from a fasta file. Build using name.exon_no

    Args:
        input_fasta (str): fasta file path containing exon
        output_fasta (str): output fasta file path containing
    """

    # create a list sorted by name and exon number
    exon_dict = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict()))
    strand_dict = {}
    exon_names, exon_seqs = gen.read_fasta(input_fasta)
    for i, name in enumerate(exon_names):
        clean_name = name.split('(')[0].split('.')[0]
        exon_id = int(name.split('(')[0].split('.')[1])
        strand = name.split('(')[1].strip(')')

        exon_dict[clean_name][exon_id] = exon_seqs[i]
        strand_dict[clean_name] = strand

    # now build the sequences
    names = []
    seqs = []
    for id in exon_dict:
        names.append(id)
        strand = strand_dict[id]
        sequence_parts = []

        if strand == "-":
            for exon_no, exon_seq in reversed(sorted(exon_dict[id].items())):
                sequence_parts.append(exon_seq)
        else:
            for exon_no, exon_seq in sorted(exon_dict[id].items()):
                sequence_parts.append(exon_seq)
        sequence = "".join(sequence_parts)
        seqs.append(sequence)

    # write to file
    gen.write_to_fasta(names, seqs, output_fasta)


def entries_to_bed(source_path, output_file, exclude_XY=None, hg38=None, NONCODE=None):
    """
    Generate a file containing the exon info from a bed file

    Args:
        source_path (str): the source path for the origin .bed file
        output_file (str): output .bed file to contain the exon info
        exclude_XY (bool): if true, exclude cases on X and Y chr
        hg38 (bool): if true, using hg38
    """
    # read the file in
    lines = gen.read_many_fields(source_path, "\t")

    with open(output_file, "w") as outfile:
        for line in lines:
            features = gen.Bed_Feature(line)
            starts = [int(start) for start in features.featureStarts.split(',')[:-1]]
            sizes = [int(size) for size in features.featureSizes.split(',')[:-1]]
            expected = features.featureCount
            if len(sizes) == expected and len(starts) == expected:
                for i, start_pos in enumerate(starts):
                    # get the end position
                    start_pos = features.start + start_pos
                    end_pos = start_pos + sizes[i]
                    # create a list of the new bed line
                    write=True
                    if hg38:
                        features.chr = features.chr.strip("chr")
                    if NONCODE:
                        features.name = features.name.split(".")[0]
                        if "NONHSAT" not in features.name:
                            write=False
                    output_list = [features.chr, start_pos, end_pos, "{0}.{1}".format(features.name, i+1), ".", features.strand, features.thickStart, features.thickEnd, ".", features.featureCount, ".", "."]
                    # only add if a transcript, used for NONCODE sequences
                    if write:
                        # add the info if exists
                        if hasattr(features, "info"):
                            output_list.extend(features.info)
                        if exclude_XY:
                            if features.chr not in ["chrX", "chrY"]:
                                outfile.write("{0}\n".format("\t".join(gen.stringify(output_list))))
                        else:
                            outfile.write("{0}\n".format("\t".join(gen.stringify(output_list))))
            else:
                print('Error in the number of features')


def extract_seqs(source_path, genome_fasta, output_bed, output_fasta, output_seq_fasta, mapping_file, codes_file, exclude_XY=None, hg38=None, NONCODE=None):
    """
    Generate a file containing the exon sequences for a given .bed file

    Args:
        source_path (str): the source path for the origin .gtf file
        genome_fasta (str): the source path for the genome fasta
        output_bed (str): output .bed file to contain the exon info
        output_fasta (str): output fasta containing sequences
        output_seq_fasta (str):
        mapping_file (str):
        codes_file (str): used for NONCODE sequences to get the lincRNA
        exclude_XY (bool): if true, exclude cases on the X and Y chr
        hg38 (bool): if true, use hg38
        NONCODE (bool): if true, using NONCODE sequences
    """

    # create the exon bed file
    full_bed = "{0}/full_{1}".format("/".join(output_bed.split('/')[:-1]), output_bed.split("/")[-1])
    entries_to_bed(source_path, full_bed, exclude_XY, hg38=hg38, NONCODE=NONCODE)
    # generate the fasta from the file
    full_exon_fasta = "{0}/full_{1}".format("/".join(output_fasta.split('/')[:-1]), output_fasta.split("/")[-1])
    fasta_from_intervals(full_bed, full_exon_fasta, genome_fasta, names=True)
    # build the sequences from the exons
    full_seq_fasta = "{0}/full_{1}".format("/".join(output_seq_fasta.split('/')[:-1]), output_seq_fasta.split("/")[-1])
    build_seqs_from_exons_fasta(full_exon_fasta, full_seq_fasta)
    length_filter_fasta = "{0}/length_filtered_{1}".format("/".join(output_seq_fasta.split('/')[:-1]), output_seq_fasta.split("/")[-1])
    ops.filter_seq_lengths(full_seq_fasta, length_filter_fasta, 200)
    # filter to only keep one transcript per gene
    unique_transcripts_fasta = "{0}/unique_gene_filtered_{1}".format("/".join(output_seq_fasta.split('/')[:-1]), output_seq_fasta.split("/")[-1])
    ops.uniquify_lincRNA_transcripts(length_filter_fasta, mapping_file, unique_transcripts_fasta)

    if NONCODE:
        # get only those that are lincRNA
        ops.get_passed_NONCODE_codes(unique_transcripts_fasta, codes_file, mapping_file, output_seq_fasta, "0001")
    else:
        # otherwise dont need the step above so copy to file
        gen.run_process(["cp", unique_transcripts_fasta, output_seq_fasta])
    # filter bed file from fasta
    ops.filter_bed_from_fasta(full_bed, output_seq_fasta, output_bed)
    # now just get the exon seqs from these entries
    fasta_from_intervals(output_bed, output_fasta, genome_fasta, names=True)


def fasta_from_intervals(bed_file, fasta_file, genome_fasta, force_strand = True, names = False):
    """
    Takes a bed file and creates a fasta file with the corresponding sequences.
    Credit: Rosina Savisaar

    Args:
        bed_file (str): the bed file path to create fasta from
        fasta_file (str): the output fasta file path
        genome_fasta (str): the file path to the genome fasta file
        names (bool): if False, the fasta record names will be generated from the sequence coordinates.
        names (bool): if True, the fasta name will correspond to whatever is in the 'name' field of the bed file
    """

    #if the index file exists, check whether the expected features are present
    genome_fasta_index = genome_fasta + '.fai'
    if(os.path.exists(genome_fasta_index)):
        bed_chrs = sorted(list(set([entry[0] for entry in gen.read_many_fields(bed_file, "\t")])))
        index_chrs = sorted(list(set([entry[0] for entry in gen.read_many_fields(genome_fasta_index, "\t")])))
        if(not set(bed_chrs).issubset(set(index_chrs))):
            gen.remove_file(genome_fasta_index)

    bedtools_args = ["bedtools", "getfasta", "-s", "-fi", genome_fasta, "-bed", bed_file, "-fo", fasta_file]
    if not force_strand:
        del bedtools_args[2]
    if names:
        bedtools_args.append("-name")
    gen.run_process(bedtools_args)
    names, seqs = gen.read_fasta(fasta_file)
    seqs = [i.upper() for i in seqs]
    gen.write_to_fasta(names, seqs, fasta_file)


def order_temp_files(files):
    """
    Get a dictionary of filelists by simulation number

    Args:
        files (list): list of files with format paths_to_file/randomFloat.simNo.ext

    Returns:
        filelist (dict): dict[simNo] = file
    """

    filelist = {}
    for file in files:
        simulation_no = int(file.split('.')[-2])
        filelist[simulation_no] = file
    return filelist
