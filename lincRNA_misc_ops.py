import generic as gen

import time
import os
import collections
import re

from Bio import SeqIO


def extract_exon_intron_fasta(input_fasta, output_exons_fasta, output_introns_fasta):
    """
    Extract exons and introns from file containing both.

    Args:
        input_fasta (str): path to input fasta file
        output_exons_fasta (str): path to output exons file
        output_intons_fasta (str): path to output introns file
    """

    try:
        # read in the file
        names, seqs = gen.read_many_fields(input_fasta)
    except:
        names, seqs = [], []
        for seq_record in SeqIO.parse(input_fasta, "fasta"):
            # this presumes the id is the last bit of the description
            names.append(seq_record.description.split(" ")[-1])
            seqs.append(str(seq_record.seq))


    # create a dictionary of entries
    seq_list = {name: seqs[i] for i, name in enumerate(names)}
    # now split each sequence into exons and introns
    exons_list = {}
    introns_list = {}
    for id in seq_list:
        exons = re.findall("([A-Z]+)", seq_list[id])
        if id == "XLOC_008301":
            print(seq_list[id])
            print(exons)
        introns = re.findall("([^A-Z]+)", seq_list[id])
        exons_list[id] = exons
        introns_list[id] = [i.upper() for i in introns]

    with open(output_exons_fasta, "w") as outfile:
        for id in exons_list:
            for i, exon in enumerate(exons_list[id]):
                outfile.write(">{0}.{1}\n{2}\n".format(id, i+1, exon))
    with open(output_introns_fasta, "w") as outfile:
        for id in introns_list:
            for i, intron in enumerate(introns_list[id]):
                outfile.write(">{0}.{1}\n{2}\n".format(id, i+1, intron))

def extract_sequences(input_fasta, output_fasta):
    """
    Extract sequences from exons file.

    Args:
        input_fasta (str): path to input fasta file
        output_fasta (str): path to output fasta file
    """

    # read the input file
    names, seqs = gen.read_fasta(input_fasta)
    # create a dictionary of entries
    exons_list = collections.defaultdict(lambda: [])
    [exons_list[name.split(".")[0]].append(seqs[i]) for i, name in enumerate(names)]
    # create dictionary containing the exons
    with open(output_fasta, "w") as outfile:
        for id in exons_list:
            seq = "".join(exons_list[id])

            if id == "XLOC_008301":
                print(id, len(seq), seq[:10], seq[-10:])
                [print(i) for i in exons_list[id]]
            outfile.write(">{0}\n{1}\n".format(id, "".join(exons_list[id])))
