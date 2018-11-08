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

    # convert the sequences to protein sequence
    cds_iupac = Seq(cds_seq[0], IUPAC.unambiguous_dna)
    cds_protein_seq = cds_iupac.translate()
    ortholog_ids = transcript_cds_orthologs.keys()
    orthologs_iupac = [Seq(ortholog) for ortholog in transcript_cds_orthologs.values()]
    ortholog_protein_seqs = [seq.translate() for seq in orthologs_iupac]

    print(cds_protein_seq)
    [print(i) for i in ortholog_protein_seqs]

    # create temp directory
    temp_dir = "temp_alignment"
    gen.remove_directory(temp_dir)
    gen.create_output_directories(temp_dir)

    # for each of the ortholog sequences
    for i, ortholog_id in enumerate(ortholog_ids):
        # write to temp .phylip file
        temp_alignment_file = "{0}/protein_alignment_{1}.fasta".format(temp_dir, random.random())
        with open(temp_alignment_file, "w") as temp_file:
            temp_file.write(">{0}\n{1}\n>{2}\n{3}\n".format(transcript_id, cds_protein_seq, ortholog_id, ortholog_protein_seqs[i]))

        muscle_exe = "../tools/muscle3.8.31_i86darwin64"
        temp_output_file = "{0}/alignment_output_{1}.fasta".format(temp_dir, random.random())
        muscle_object = MuscleCommandline(muscle_exe, input = temp_alignment_file, out = temp_output_file)
        stdout, stderr = muscle_object()
        print(stdout, stderr)

        file_from_muscle = open(temp_output_file)
        muscle_string = "".join(file_from_muscle)
        muscle_string = re.sub("([A-Z\-])\n([A-Z\-])","\\1\\2", muscle_string)
        print(muscle_string)
        aligned_prot_sequences = re.findall("^[A-Z\-]+(?=\n)", muscle_string, re.MULTILINE)
        print(aligned_prot_sequences)
