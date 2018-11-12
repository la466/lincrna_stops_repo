import generic as gen
import file_ops as fo
import sequence_ops as sequo
import time
import os
import collections
import random
import re
import sys
from useful_motif_sets import codon_map, stops
from Bio.Phylo.PAML import codeml
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


def get_conservation(transcript_list, output_file, max_dS_threshold = None, max_omega_threshold = None):
    """
    Get the conversation for a list of sequences and only keep those that pass

    Args:
        transcript_list (dict): dict containing transcript id, the cds and the ortholog seqs
        output_file (str): path to output file
        max_dS_threshold (float): if set, pass in the dS threshold you wish alignments to be below
        max_omega_threshold (float): if set, pass in the omega threshold you wish alignments to be below
    """

    print("Getting the most conserved ortholog for each transcript...")

    temp_dir = "temp_conservation_files"
    gen.create_output_directories(temp_dir)
    # get a list of the transcript ids
    transcript_ids = list(transcript_list.keys())
    # output_filelist = run_conservation_check(transcript_ids, transcript_list, max_dS_threshold, max_omega_threshold, temp_dir)
    outputs = gen.run_parallel_function(transcript_ids, [transcript_list, max_dS_threshold, max_omega_threshold, temp_dir], run_conservation_check)
    # remove the old output file if there is one
    gen.remove_file(output_file)
    # now concat the output files
    args = ["cat"]
    [args.append(i) for i in outputs]
    gen.run_process(args, file_for_output = output_file)
    gen.remove_directory(temp_dir)


def run_conservation_check(input_list, transcript_list, max_dS_threshold, max_omega_threshold, temp_dir):
    """
    Wrapper to run the conservation check in parallel

    Args:
        input_list (list): list of transcript ids to iterate over
        transcript_list (dict): dict containing transcript id, the cds and the ortholog seqs
        output_file (str): path to output file
        max_dS_threshold (float): if set, pass in the dS threshold you wish alignments to be below
        max_omega_threshold (float): if set, pass in the omega threshold you wish alignments to be below
    """
    # create a list to keep temporary outputs
    temp_filelist = []
    temp_instance_dir = "temp_codeml_dir.{0}".format(random.random())
    gen.create_output_directories(temp_instance_dir)

    if input_list:
        temp_file = "{0}/best_ortholog_match.{1}.bed".format(temp_dir, random.random())
        temp_filelist.append(temp_file)
        with open(temp_file, "w") as outfile:
            # get best ortholog for each transcript
            for i, transcript_id in enumerate(input_list):
                print("{0}/{1}".format(i+1, len(input_list)))
                ortholog_id = check_conservation(transcript_id, transcript_list[transcript_id][0], transcript_list[transcript_id][1], temp_instance_dir,  max_dS_threshold = max_dS_threshold, max_omega_threshold = max_omega_threshold)
                if ortholog_id:
                    outfile.write("{0}\t{1}\n".format(transcript_id, ortholog_id))
    gen.remove_directory(temp_instance_dir)
    return temp_filelist


def check_conservation(transcript_id, cds_seq, transcript_cds_orthologs, temp_dir, max_dS_threshold = None, max_omega_threshold = None):
    """
    Check the conservation of a sequence.

    Args:
        transcript_id (str): id of the transcript in question
        cds_seq (str): nucleotide sequence of the transcript
        transcript_cds_orthologs (dict): dict containing the sequences of orthologous sequences
        max_dS_threshold (float): if set, the dS threshold you wish alignments to be below
        max_omega_threshold (float): if set, the omega threshold you wish alignments to be below
    Returns:
        best_ortholog_id (str): id of the ortholog that gives the most conserved alignment
    """

    muscle_exe = "../tools/muscle3.8.31_i86{0}64".format(sys.platform)
    if not os.path.isfile(muscle_exe):
        print("Could not find the MUSCLE exe {0}...".format(muscle_exe))
        raise Exception

    # setup the muscle alignment
    alignment_functions = Alignment_Functions(muscle_exe)

    # convert the sequences to protein sequence
    cds_iupac = Seq(cds_seq[0], IUPAC.unambiguous_dna)
    cds_protein_seq = cds_iupac.translate()
    ortholog_ids = list(transcript_cds_orthologs.keys())
    orthologs_iupac = [Seq(ortholog) for ortholog in transcript_cds_orthologs.values()]
    ortholog_protein_seqs = [seq.translate() for seq in orthologs_iupac]
    # set up lists to hold the transcript_ds scores and omega scores
    dS_scores = []
    omega_scores = []
    # for each of the ortholog sequences
    for i, ortholog_id in enumerate(ortholog_ids):
        # align the sequences
        alignment_functions.align_seqs(cds_protein_seq, ortholog_protein_seqs[i], transcript_id, ortholog_id)
        # extract the alignments
        alignment_functions.extract_alignments()
        # now we want to get the nucleotide sequences for the alignments
        aligned_sequences = alignment_functions.revert_alignment_to_nucleotides(input_seqs = [cds_seq[0], transcript_cds_orthologs[ortholog_id]])
        # clean up the files
        alignment_functions.cleanup()
        # we have the aligned sequences
        #write the nucleotide alignment to a phylip file
        aligned_sequences_iupac = [Seq("".join(i),IUPAC.unambiguous_dna) for i in aligned_sequences]
        alignment = MultipleSeqAlignment([SeqRecord(aligned_sequences_iupac[0], id = "seq"), SeqRecord(aligned_sequences_iupac[1], id = "orth_seq")])
        temp_phylip_file = "{0}/{1}_{2}.phy".format(temp_dir, transcript_id, ortholog_id)
        temp_output_file = "{0}/phy_{1}_{2}.out".format(temp_dir, transcript_id, ortholog_id)
        fo.write_to_phylip(alignment, temp_phylip_file)
        # # run paml on sequences
        paml = sequo.PAML_Functions(input_file = temp_phylip_file, output_file = temp_output_file, working_dir = temp_dir)
        # run codeml
        codeml_output = paml.run_codeml()
        dS_scores.append(codeml_output["NSsites"][0]["parameters"]["dS"])
        omega_scores.append(codeml_output["NSsites"][0]["parameters"]["omega"])

    # get the minimum dS and omega scores
    min_dS = min(dS_scores)
    min_omega = min(omega_scores)

    # get the alignment with the min dS
    # do it this way so that if you dont have the filters, you just keep the most conserved
    best_ortholog_id = ortholog_ids[dS_scores.index(min_dS)]
    # check that there is an ortholog with omega < 0.5
    if min_omega >= max_omega_threshold:
        best_ortholog_id = None
    if min_dS >= max_dS_threshold:
        best_ortholog_id = None

    return best_ortholog_id



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


    def extract_alignments(self, alignment_output_file = None):
        """
        Extract the alignments from a muscle output file

        Args:
            input_file (str): if set, the path to and alignment output file
        """

        # set the input file if not given
        if not alignment_output_file:
            alignment_output_file = self.alignment_output_file
        alignments = extract_alignments(alignment_output_file)
        self.protein_alignments = alignments


    def revert_alignment_to_nucleotides(self, input_seqs, input_alignments = None):
        """
        Revert the aligned protein sequences back to nucleotides

        Args:
            input_alignments (list): list containing the two protein alignments
        """

        if not input_alignments:
            input_alignments = self.protein_alignments
        aligned_seqs = [revert_alignment_to_nucleotides(seq, input_alignments[i]) for i, seq in enumerate(input_seqs)]
        self.aligned_seqs = aligned_seqs
        return aligned_seqs




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

    Returns:
        alignments (list): list containing [original_cds_alignment, ortholog_cds_alignment]
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


def revert_alignment_to_nucleotides(input_seq, input_alignment):
    """
    Given a protein alignment, convert the alignment back to the nucleotide
    sequence that it was originally aligned from.

    Args:
        input_seq (str): input nucleotide sequence
        input_alignment (str): input protein alignment
    """

    # break each down into list of items
    input_nts = list(input_seq)
    input_alignment_items = list(input_alignment)
    # create a list to hold the aligned sequence
    aligned_sequence = []
    sequence_position = 0

    for i, amino_acid in enumerate(input_alignment_items):
        if amino_acid == "-":
            aligned_sequence.append("---")
        else:
            codon = "".join([input_nts[sequence_position], input_nts[sequence_position+1], input_nts[sequence_position+2]])
            if codon in stops or amino_acid != codon_map[codon]:
                print("There is an error in the alignment...")
                print(input_seq)
                print(input_alignment)
                print(i, amino_acid, codon, codon_map[codon])
                raise Exception
            else:
                aligned_sequence.append(codon)
                sequence_position += 3
    aligned_sequence = "".join(aligned_sequence)
    return aligned_sequence
