import generic as gen
import ops as ops
import numpy as np
import itertools as it
import os
import re
import collections

def non_coding_exons(genome_fasta, gtf, output_directory, clean_run=None):
    """
    Wrapper for looking at non coding exons

    Args:
        genome_fasta (str): path to the genome fasta file
        gtf (str): path to the genome gtf file
        output_directory (str): path to the output directory
        clean_run(bool): if true, remove the processing sequence files and start fresh
    """

    sequence_output_directory = "{0}/sequence_files".format(output_directory)
    if clean_run:
        gen.remove_directory(sequence_output_directory)

    # create the output directory if it doesnt exist
    gen.create_output_directories(sequence_output_directory)

    # get the coding and non coding exons
    coding_exons_bed = "{0}/coding_exons.bed".format(output_directory)
    non_coding_exons_bed = "{0}/non_coding_exons.bed".format(output_directory)
    coding_exons_fasta = "{0}/coding_exons.fasta".format(output_directory)
    non_coding_exons_fasta = "{0}/non_coding_exons.fasta".format(output_directory)
    if clean_run or not os.path.isfile(non_coding_exons_fasta):
        ops.get_non_coding_exons(genome_fasta, gtf, coding_exons_bed, non_coding_exons_bed, coding_exons_fasta, non_coding_exons_fasta, sequence_output_directory)
