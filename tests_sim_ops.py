import sim_ops_containers as simopc
import sim_ops as so
import generic as gen
import numpy as np
import unittest
import re

class Test_Sim_Ops(unittest.TestCase):

    def test_sim_orf_length(self):
        input_file = "test_data/sim_ops/test_sim_orf_length/input.fasta"
        expected_file = "test_data/sim_ops/test_sim_orf_length/expected.txt"
        observed_file = "test_data/sim_ops/test_sim_orf_length/observed.txt"
        gen.remove_file(observed_file)
        names, seqs = gen.read_fasta(input_file)

        simulations = 5
        seeds = list(range(simulations))
        seq_seeds = [[10,20,30],[11,21,31],[12,22,32],[13,23,33],[14,24,34]]
        simopc.sim_orf_length(input_file, simulations, observed_file, parallel = False, seeds=seeds, seq_seeds=seq_seeds)
        expected = gen.read_many_fields(expected_file, ",")
        observed = gen.read_many_fields(observed_file, ",")
        self.assertEqual(expected, observed)

    def test_sim_stop_count(self):
        input_file = "test_data/sim_ops/test_sim_stop_count/input.fasta"
        expected_file = "test_data/sim_ops/test_sim_stop_count/expected.txt"
        observed_file = "test_data/sim_ops/test_sim_stop_count/observed.txt"
        gen.remove_file(observed_file)
        names, seqs = gen.read_fasta(input_file)

        simulations = 5
        seeds = list(range(simulations))
        seq_seeds = [[10,20,30],[11,21,31],[12,22,32],[13,23,33],[14,24,34]]
        simopc.sim_stop_count(input_file, simulations, observed_file, parallel = False, seeds=seeds, seq_seeds=seq_seeds)
        expected = gen.read_many_fields(expected_file, ",")
        observed = gen.read_many_fields(observed_file, ",")
        self.assertEqual(expected, observed)

    def test_sim_cds_seqs(self):
        input_file = "test_data/sim_ops/test_sim_cds_seqs/input.fasta"
        expected_file = "test_data/sim_ops/test_sim_cds_seqs/expected.fasta"
        seq_list = gen.read_fasta(input_file)[1]
        expected = gen.read_fasta(expected_file)[1]
        codons = [re.findall(".{3}", seq) for seq in seq_list]
        codon_list = [codon_set[1:-1] for codon_set in codons]
        starts = [i[0] for i in codons]
        stops = [i[-1] for i in codons]
        observed = so.sim_cds_seqs(codon_list, starts, stops, seed=1)
        self.assertEqual(expected, observed)
