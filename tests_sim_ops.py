from sim_ops_containers import *
import numpy as np
import unittest

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
        sim_orf_length(input_file, simulations, observed_file, parallel = False, seeds=seeds, seq_seeds=seq_seeds)
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
        sim_stop_count(input_file, simulations, observed_file, parallel = False, seeds=seeds, seq_seeds=seq_seeds)
        expected = gen.read_many_fields(expected_file, ",")
        observed = gen.read_many_fields(observed_file, ",")
        self.assertEqual(expected, observed)
