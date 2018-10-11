from seq_ops import *
import numpy as np
import unittest

class Test_Seq_Ops(unittest.TestCase):

    def test_calc_seq_gc(self):
        input_file = "test_data/seq_ops/test_calc_seq_gc/input.fasta"
        names, seqs = gen.read_fasta(input_file)
        expected = [np.divide(4,8), np.divide(3,9), np.divide(0,6), np.divide(7,13)]
        observed = [calc_seq_gc(seq) for seq in seqs]
        self.assertEqual(expected, observed)

    def test_fasta_to_dict(self):
        input_file = "test_data/seq_ops/test_fasta_to_dict/input.fasta"
        expected = {"seq1": "ACGGACGAAGCGAGATTTAAAA", "seq2": "ACGCGAGGAAGCGTACGTCGATAACGAGA", "seq3": "GAGTTGTGACGAG"}
        observed = fasta_to_dict(input_file)
        self.assertEqual(expected, observed)

    def test_generate_nt_matched_seq(self):
        input_file = "test_data/seq_ops/test_generate_nt_matched_seq/input.fasta"
        expected_file = "test_data/seq_ops/test_generate_nt_matched_seq/expected.fasta"
        input_names, input_seqs = gen.read_fasta(input_file)
        expected_names, expected = gen.read_fasta(expected_file)
        dinucleotides = ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]
        dinucleotide_probabilities = [0.01, 0.08, 0.06, 0.1, 0.02, 0.15, 0.03, 0.05, 0.2, 0.02, 0.015, 0.015, 0.0625, 0.0625, 0.0825, 0.0425]
        nucleotides = ["A", "C", "G", "T"]
        nucleotide_probabilities = [0.1,0.4,0.25,0.25]
        seed = 5
        observed = []
        for seq in input_seqs:
            observed.append(generate_nt_matched_seq(seq, dinucleotides, dinucleotide_probabilities, nucleotides, nucleotide_probabilities, seed=seed))
        self.assertEqual(expected, observed)

    def test_get_dinucleotide_content(self):
        input_file = "test_data/seq_ops/test_get_dinucleotide_content/input.fasta"
        names, seqs = gen.read_fasta(input_file)
        expected = {
            "AA": np.divide(3,43), "AC": np.divide(2,43), "AG": np.divide(5,43), "AT": np.divide(2,43),
            "CA": np.divide(2,43), "CC": np.divide(2,43), "CG": np.divide(3,43), "CT": np.divide(2,43),
            "GA": np.divide(4,43), "GC": np.divide(3,43), "GG": np.divide(3,43), "GT": np.divide(3,43),
            "TA": np.divide(3,43), "TC": np.divide(2,43), "TG": np.divide(2,43), "TT": np.divide(2,43)
        }
        observed = get_dinucleotide_content(seqs)
        self.assertEqual(expected, observed)

    def test_get_dinucleotide_content_as_string(self):
        input_file = "test_data/seq_ops/test_get_dinucleotide_content/input.fasta"
        names, seqs = gen.read_fasta(input_file)
        expected = {
            "AA": np.divide(3,48), "AC": np.divide(3,48), "AG": np.divide(5,48), "AT": np.divide(2,48),
            "CA": np.divide(2,48), "CC": np.divide(2,48), "CG": np.divide(4,48), "CT": np.divide(2,48),
            "GA": np.divide(4,48), "GC": np.divide(3,48), "GG": np.divide(4,48), "GT": np.divide(4,48),
            "TA": np.divide(4,48), "TC": np.divide(2,48), "TG": np.divide(2,48), "TT": np.divide(2,48)
        }
        observed = get_dinucleotide_content(seqs, as_string=True)
        self.assertEqual(expected, observed)

    def test_get_longest_orf(self):
        input_file = "test_data/seq_ops/test_get_longest_orf/input.fasta"
        names, seqs = gen.read_fasta(input_file)
        expected = [21,30,10,np.nan]
        observed = []
        for seq in seqs:
            observed.append(get_longest_orf(seq))
        self.assertEqual(expected, observed)

    def test_get_longest_orf_strict_stop(self):
        input_file = "test_data/seq_ops/test_get_longest_orf_strict_stop/input.fasta"
        names, seqs = gen.read_fasta(input_file)
        expected = [21,30,np.nan,np.nan]
        observed = []
        for seq in seqs:
            observed.append(get_longest_orf(seq, strict_stop=True))
        self.assertEqual(expected, observed)

    def test_get_longest_orfs_dict(self):
        input_file = "test_data/seq_ops/test_get_longest_orfs_dict/input.fasta"
        seqs = fasta_to_dict(input_file)
        expected = {"seq1|21": 21, "seq2|30": 30, "seq3|10": 10, "seq4|na": np.nan}
        observed = get_longest_orfs(seqs)
        self.assertEqual(expected, observed)

    def test_get_longest_orfs_list(self):
        input_file = "test_data/seq_ops/test_get_longest_orfs_list/input.fasta"
        names, seqs = gen.read_fasta(input_file)
        expected = [21,30,10,np.nan]
        observed = get_longest_orfs(seqs)
        self.assertEqual(expected, observed)

    def test_get_nucleotide_content(self):
        input_file = "test_data/seq_ops/test_get_nucleotide_content/input.fasta"
        names, seqs = gen.read_fasta(input_file)
        expected = {"A": np.divide(7,18), "C": np.divide(3,18), "G": np.divide(5,18), "T": np.divide(3,18)}
        observed = get_nucleotide_content(seqs)
        self.assertEqual(expected, observed)

    def test_get_stop_count(self):
        input_file = "test_data/seq_ops/test_get_stop_count/input.fasta"
        names, seqs = gen.read_fasta(input_file)
        expected = [1,5,2,0,3]
        observed = []
        for seq in seqs:
            observed.append(get_stop_count(seq))
        self.assertEqual(expected, observed)

    def test_get_stop_counts_list(self):
        input_file = "test_data/seq_ops/test_get_stop_counts_list/input.fasta"
        names, seqs = gen.read_fasta(input_file)
        expected = [1,6,7,0,10]
        observed = get_stop_counts(seqs)
        self.assertEqual(expected, observed)

    def test_get_stop_counts_dict(self):
        input_file = "test_data/seq_ops/test_get_stop_counts_dict/input.fasta"
        seqs = fasta_to_dict(input_file)
        expected = {"seq1": 1, "seq2": 6, "seq3": 7, "seq4": 0, "seq5": 10}
        observed = get_stop_counts(seqs)
        self.assertEqual(expected, observed)

    def test_get_unique_seqs(self):
        input_file = "test_data/seq_ops/test_get_unique_seqs/input.fasta"
        expected_file = "test_data/seq_ops/test_get_unique_seqs/expected.fasta"
        names, seqs = gen.read_fasta(input_file)
        expected = fasta_to_dict(expected_file)
        observed = get_unique_seqs(names, seqs)
        self.assertEqual(expected, observed)
