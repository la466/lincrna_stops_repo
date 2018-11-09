import generic as gen
import file_ops as fo
import time
import os
import collections
import random
import re
from Bio.Seq import Seq
# from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Align.Applications import MuscleCommandline
# from Bio.SeqRecord import SeqRecord
# from Bio.Align import MultipleSeqAlignment
# from Bio import AlignIO


def get_conservation(transcript_list):
    """
    Get the conversation for a list of sequences and only keep those that pass

    Args:
        transcript_list (dict): dict containing transcript id, the cds and the ortholog seqs
    """

    for transcript_id in transcript_list:
        check_conservation(transcript_id, transcript_list[transcript_id][0], transcript_list[transcript_id][1])


def check_conservation(transcript_id, cds_seq, transcript_cds_orthologs):
    """
    Check the conservation of a sequence.
    """

    # setup the muscle alignment
    muscle_exe = "../tools/muscle3.8.31_i86darwin64"
    alignment_functions = Alignment_Functions(muscle_exe)

    # convert the sequences to protein sequence
    cds_iupac = Seq(cds_seq[0], IUPAC.unambiguous_dna)
    cds_protein_seq = cds_iupac.translate()
    ortholog_ids = transcript_cds_orthologs.keys()
    orthologs_iupac = [Seq(ortholog) for ortholog in transcript_cds_orthologs.values()]
    ortholog_protein_seqs = [seq.translate() for seq in orthologs_iupac]

    # for each of the ortholog sequences
    for i, ortholog_id in enumerate(ortholog_ids):
        # align the sequences
        alignment_functions.align_seqs(cds_protein_seq, ortholog_protein_seqs[i], transcript_id, ortholog_id)
        # extract the alignments
        alignment_functions.extract_alignments()
        # now we want to get the nucleotide sequences for the alignments
        alignment_functions.revert_alignment_to_nucleotides()
        # clean up the files
        alignment_functions.cleanup()




class Alignment_Functions(object):

    def __init__(self, muscle_exe):
        self.muscle_exe = muscle_exe

    def align_seqs(self, seq1, seq2, seq1_id = None, seq2_id = None, temp_input_file = None, temp_output_file = None):
        """
        Align two protein sequences using Muscle.

        Args:
            seq1 (str): protein sequence 1
            seq2 (str): protein sequence 2
            seq_id1 (str): if set, the id for sequence 1
            seq_id2 (str): if set, the id for sequence 2
            temp_input_file (str): if set, the path to the alignment input file
            temp_output_file (str): if set, the path to the alignment output file
        """

        # add the sequences to the object
        self.input_seqs = [seq1, seq2]
        alignment_input_file, alignment_output_file = align_sequences(self.muscle_exe, self.input_seqs[0], self.input_seqs[1], seq1_id, seq2_id, temp_input_file = temp_input_file, temp_output_file = temp_output_file)
        self.alignment_input_file = alignment_input_file
        self.alignment_output_file = alignment_output_file


    def cleanup(self):
        """
        Cleanup the output files
        """
        gen.remove_file(self.alignment_input_file)
        gen.remove_file(self.alignment_output_file)


    def extract_alignments(self, input_file = None):
        """
        Extract the alignments from a muscle output file

        Args:
            input_file (str): if set, the path to and alignment output file
        """
        # set the input file if not given
        if not input_file:
            input_file = self.alignment_input_file
        alignments = extract_alignments(input_file)
        self.protein_alignments = alignments


    def revert_alignment_to_nucleotides(self):
        pass



def align_sequences(muscle_exe, seq1, seq2, seq1_id = None, seq2_id = None, temp_input_file = None, temp_output_file = None):
    """
    Align two protein sequences using Muscle.

    Args:
        muscle_exe (str): path to muscle tool
        seq1 (str): protein sequence 1
        seq2 (str): protein sequence 2
        seq_id1 (str): if set, the id for sequence 1
        seq_id2 (str): if set, the id for sequence 2
        temp_input_file (str): if set, the path to the alignment input file
        temp_output_file (str): if set, the path to the alignment output file
    """

    # create temp files for running alignment
    temp_dir = "temp_alignment"
    gen.create_output_directories(temp_dir)
    # create the random alignment files
    random_alignment = random.random()
    if not temp_input_file:
        temp_input_file = "{0}/protein_alignment_input_{1}.fasta".format(temp_dir, random_alignment)
    if not temp_output_file:
        temp_output_file = "{0}/protein_alignment_output_{1}.fasta".format(temp_dir, random_alignment)
    # in case the sequence ids are not set
    if not seq1_id:
        seq1_id = "seq_id_{0}_1".format(random.random())
    if not seq2_id:
        seq2_id = "{0}_2".format(seq1_id[:-2])
    # write the temporary alignment file
    with open(temp_input_file, "w") as temp_file:
        temp_file.write(">{0}\n{1}\n>{2}\n{3}\n".format(seq1_id, seq1, seq2_id, seq2))
    # run muscle alignment
    muscle_output = MuscleCommandline(muscle_exe, input = temp_input_file, out = temp_output_file)
    # get object
    muscle_output()

    return temp_input_file, temp_output_file


def extract_alignments(input_file):
    """
    Extract the alignments from a muscle output file

    Args:
        input_file (str): the path to and alignment output file
    """

    with open(input_file) as input_content:
        content = input_content.read()
        # remove any gaps
        muscle_alignments = "".join(content)
        # now join the lines of the alignment if there are multiple lines
        muscle_alignments = re.sub("([A-Z\-])\n([A-Z\-])","\\1\\2", muscle_alignments)
        # now get the sequences, searching each line
        alignments = re.findall("^[A-Z\-]+(?=\n)", muscle_alignments, re.MULTILINE)

    return alignments
