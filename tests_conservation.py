from conservation import *
import numpy as np
import unittest
import collections
import sys

muscle_exe = "../tools/muscle3.8.31_i86{0}64".format(sys.platform)

class Tests(unittest.TestCase):

    def test_add_family_hits(self):
        input_file = "test_data/conservation/test_add_family_hits/input.bed"
        inputs = gen.read_many_fields(input_file, "\t")
        families = [[]]
        expected_results = [["ENST0003", "ENST0013"], ["ENST0004", "ENST0005"], ["ENST0004", "ENST0008"], ["ENST0005", "ENST0004"], ["ENST0006", "ENST0009"], ["ENST0008", "ENST0004"], ["ENST0008", "ENST0011"], ["ENST0009", "ENST0006"], ["ENST0010", "ENST0012"], ["ENST0011", "ENST0008"], ["ENST0012", "ENST0010"], ["ENST0013", "ENST0003"]]
        expected_families = [["ENST0001", "ENST0002", "ENST0003"]]
        id = "ENST0001"
        observed_results, observed_families = add_family_hits(inputs, families, id)
        observed_families = [sorted(i) for i in observed_families]
        self.assertEqual(expected_results, observed_results)
        self.assertEqual(expected_families, observed_families)
        inputs = [["ENST0003", "ENST0013"], ["ENST0004", "ENST0005"], ["ENST0004", "ENST0008"], ["ENST0005", "ENST0004"], ["ENST0006", "ENST0009"], ["ENST0008", "ENST0004"], ["ENST0008", "ENST0011"], ["ENST0009", "ENST0006"], ["ENST0010", "ENST0012"], ["ENST0011", "ENST0008"], ["ENST0012", "ENST0010"], ["ENST0013", "ENST0003"]]
        families = [["ENST0001", "ENST0002", "ENST0003"], []]
        expected_results = [["ENST0003", "ENST0013"], ["ENST0006", "ENST0009"], ["ENST0008", "ENST0011"], ["ENST0009", "ENST0006"], ["ENST0010", "ENST0012"], ["ENST0011", "ENST0008"], ["ENST0012", "ENST0010"], ["ENST0013", "ENST0003"]]
        expected_families = [["ENST0001", "ENST0002", "ENST0003"], ["ENST0004", "ENST0005", "ENST0008"]]
        id = "ENST0004"
        observed_results, observed_families = add_family_hits(inputs, families, id)
        observed_families = [sorted(i) for i in observed_families]
        self.assertEqual(expected_results, observed_results)
        self.assertEqual(expected_families, observed_families)
        inputs = [["ENST0003", "ENST0013"], ["ENST0006", "ENST0009"], ["ENST0008", "ENST0011"], ["ENST0009", "ENST0006"], ["ENST0010", "ENST0012"], ["ENST0011", "ENST0008"], ["ENST0012", "ENST0010"], ["ENST0013", "ENST0003"]]
        families = [["ENST0001", "ENST0002", "ENST0003"], ["ENST0004", "ENST0005", "ENST0008"], ["ENST0006", "ENST0009"]]
        expected_results = [["ENST0003", "ENST0013"], ["ENST0008", "ENST0011"], ["ENST0010", "ENST0012"], ["ENST0011", "ENST0008"], ["ENST0012", "ENST0010"], ["ENST0013", "ENST0003"]]
        expected_families = [["ENST0001", "ENST0002", "ENST0003"], ["ENST0004", "ENST0005", "ENST0008"], ["ENST0006", "ENST0009"]]
        id = "ENST0009"
        observed_results, observed_families = add_family_hits(inputs, families, id)
        observed_families = [sorted(i) for i in observed_families]
        self.assertEqual(expected_results, observed_results)
        self.assertEqual(expected_families, observed_families)

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

    def test_group_ids_into_families(self):
        input_file1 = "test_data/conservation/test_group_ids_into_families/input1.bed"
        input_file2 = "test_data/conservation/test_group_ids_into_families/input2.bed"
        expected_file1 = "test_data/conservation/test_group_ids_into_families/expected1.bed"
        expected_file2 = "test_data/conservation/test_group_ids_into_families/expected2.bed"
        inputs = gen.read_many_fields(input_file1, "\t")
        expected = [sorted(i) for i in gen.read_many_fields(expected_file1, "\t")].sort(key=lambda x: x[0])
        observed = [sorted(i) for i in group_ids_into_families(inputs)].sort(key=lambda x: x[0])
        self.assertEqual(expected, observed)
        inputs = gen.read_many_fields(input_file2, "\t")
        expected = [sorted(i) for i in gen.read_many_fields(expected_file2, "\t")].sort(key=lambda x: x[0])
        observed = [sorted(i) for i in group_ids_into_families(inputs)].sort(key=lambda x: x[0])
        self.assertEqual(expected, observed)


    def test_remove_self_blast_hits(self):
        input_file = "test_data/conservation/test_remove_self_blast_hits/input.csv"
        expected_file = "test_data/conservation/test_remove_self_blast_hits/expected.csv"
        expected = gen.read_many_fields(expected_file, ",")
        entries = gen.read_many_fields(input_file, ",")
        observed = remove_self_blast_hits(entries)
        self.assertEqual(expected, observed)

    def test_revert_alignment_to_nucleotides(self):
        input_seqs_file = "test_data/conservation/test_revert_alignment_to_nucleotides/input_seqs.fasta"
        input_alignments_file = "test_data/conservation/test_revert_alignment_to_nucleotides/input_protein_alignment.fasta"
        expected_file = "test_data/conservation/test_revert_alignment_to_nucleotides/expected.fasta"
        expected = gen.read_fasta(expected_file)[1]
        input_seqs = gen.read_fasta(input_seqs_file)[1]
        input_alignments = gen.read_fasta(input_alignments_file)[1]
        observed = [revert_alignment_to_nucleotides(seq, input_alignments[i]) for i, seq in enumerate(input_seqs)]
        self.assertEqual(expected, observed)
