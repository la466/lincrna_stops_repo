from file_ops import *
import unittest

class Test_File_Ops(unittest.TestCase):

    def test_build_seqs_from_exons_fasta(self):
        input_file = "test_data/file_ops/test_build_seqs_from_exons_fasta/input.fasta"
        expected = "test_data/file_ops/test_build_seqs_from_exons_fasta/expected.fasta"
        observed = "test_data/file_ops/test_build_seqs_from_exons_fasta/observed.fasta"
        gen.remove_file(observed)
        build_seqs_from_exons_fasta(input_file, observed)
        expected = gen.read_fasta(expected)
        observed = gen.read_fasta(observed)
        self.assertEqual(observed,expected)

    def test_entries_to_bed(self):
        input_file = "test_data/file_ops/test_entries_to_bed/input.bed"
        expected = "test_data/file_ops/test_entries_to_bed/expected.bed"
        observed = "test_data/file_ops/test_entries_to_bed/observed.bed"
        gen.remove_file(observed)
        entries_to_bed(input_file, observed)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)

    def test_fasta_from_intervals(self):
        observed = "test_data/file_ops/test_fasta_from_intervals/observed_converted_fasta.fasta"
        bed_file = "test_data/file_ops/test_fasta_from_intervals/test_bed_for_fasta_conversion.bed"
        expected = "test_data/file_ops/test_fasta_from_intervals/expected_converted_fasta.fasta"
        gen.remove_file(observed)
        fasta_from_intervals(bed_file, observed, "test_data/file_ops/test_fasta_from_intervals/test_genome.fa")
        expected = gen.read_fasta(expected)
        observed = gen.read_fasta(observed)
        self.assertEqual(observed,expected)

    def test_order_temp_files(self):
        input = ["./path/0.34242.3.txt", "./path/0.321111.1.txt", "./path/0.983243.4.txt", "./path/0.2131313.2.txt"]
        expected = {3: "./path/0.34242.3.txt", 1: "./path/0.321111.1.txt", 4: "./path/0.983243.4.txt", 2: "./path/0.2131313.2.txt"}
        observed = order_temp_files(input)
        self.assertEqual(expected, observed)
