from seq_ops import *
import unittest

class Test_Seq_Ops(unittest.TestCase):

    def test_get_shortest_orf(self):
        input_file = "test_data/seq_ops/test_get_shortest_orf/input.fasta"
        names, seqs = gen.read_fasta(input_file)
        expected = [9,30,10,"na"]
        observed = []
        for seq in seqs:
            observed.append(get_shortest_orf(seq))
        self.assertEqual(expected, observed)

    def test_get_shortest_orf_strict_stop(self):
        input_file = "test_data/seq_ops/test_get_shortest_orf_strict_stop/input.fasta"
        names, seqs = gen.read_fasta(input_file)
        expected = [9,30,"na","na"]
        observed = []
        for seq in seqs:
            observed.append(get_shortest_orf(seq, strict_stop=True))
        self.assertEqual(expected, observed)

    def test_fasta_to_dict(self):
        input_file = "test_data/seq_ops/test_fasta_to_dict/input.fasta"
        expected = {"seq1": "ACGGACGAAGCGAGATTTAAAA", "seq2": "ACGCGAGGAAGCGTACGTCGATAACGAGA", "seq3": "GAGTTGTGACGAG"}
        observed = fasta_to_dict(input_file)
        self.assertEqual(expected, observed)
