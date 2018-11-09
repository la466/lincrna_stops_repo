from conservation import *
import numpy as np
import unittest
import collections

muscle_exe = "../tools/muscle3.8.31_i86darwin64"

class Tests(unittest.TestCase):

    def test_align_sequences(self):
        input_fasta1 = "test_data/conservation/test_align_sequences/input1.fasta"
        input_fasta2 = "test_data/conservation/test_align_sequences/input2.fasta"
        input_fasta3 = "test_data/conservation/test_align_sequences/input3.fasta"
        inputs = [input_fasta1, input_fasta2, input_fasta3]
        expected_file = "test_data/conservation/test_align_sequences/expected.fasta"
        observed_dir = "test_data/conservation/test_align_sequences"
        expected_seqs = gen.read_fasta(expected_file)[1]
        expected = [[expected_seqs[i], expected_seqs[i+1]] for i in range(0, len(expected_seqs), 2)]
        observed = []
        input_names1, input_seqs1 = gen.read_fasta(input_fasta1)
        input_names2, input_seqs2 = gen.read_fasta(input_fasta2)
        for i, input in enumerate(inputs):
            temp_input_file = "{0}/temp_input{1}.fasta".format(observed_dir, i)
            observed_file = "{0}/observed{1}.fasta".format(observed_dir, i)
            gen.remove_file(observed_file)
            input_names, input_seqs = gen.read_fasta(input)
            align_sequences(muscle_exe, input_seqs[0], input_seqs[1], input_names[0], input_names[1], temp_input_file = temp_input_file, temp_output_file = observed_file)
            observed.append(gen.read_fasta(observed_file)[1])
            gen.remove_file(temp_input_file)
        self.assertEqual(expected, observed)

    def test_extract_alignments(self):
        input_fasta1 = "test_data/conservation/test_extract_alignments/input1.fasta"
        input_fasta2 = "test_data/conservation/test_extract_alignments/input2.fasta"
        input_fasta3 = "test_data/conservation/test_extract_alignments/input3.fasta"
        inputs = [input_fasta1, input_fasta2, input_fasta3]
        expected_file = "test_data/conservation/test_extract_alignments/expected.fasta"
        expected_seqs = gen.read_fasta(expected_file)[1]
        expected = [[expected_seqs[i], expected_seqs[i+1]] for i in range(0, len(expected_seqs), 2)]
        observed = []
        for input in inputs:
            observed.append(extract_alignments(input))
        self.assertEqual(expected, observed)
