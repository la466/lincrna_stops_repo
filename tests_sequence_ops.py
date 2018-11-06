from sequence_ops import *
import numpy as np
import unittest
import collections

class Tests(unittest.TestCase):

    def test_filter_cds(self):
        input_file = "test_data/sequence_ops/test_filter_cds/input.fasta"
        expected_file = "test_data/sequence_ops/test_filter_cds/expected.fasta"
        observed_file = "test_data/sequence_ops/test_filter_cds/observed.fasta"
        gen.remove_file(observed_file)
        filter_cds(input_file, observed_file)
        # use an arbritary delimiter
        expected = gen.read_many_fields(expected_file, "#")
        observed = gen.read_many_fields(observed_file, "#")
        self.assertEqual(expected, observed)

    def test_filter_gtf_features(self):
        input_file = "test_data/sequence_ops/test_filter_gtf_features/input.gtf"
        expected_file = "test_data/sequence_ops/test_filter_gtf_features/expected.bed"
        input_list = ["ENST00000456328", "ENST0003246"]
        expected = gen.read_many_fields(expected_file, "\t")
        observed = filter_gtf_features(input_list, input_file, filter_transcripts = True)
        self.assertEqual(expected, observed)

    def test_filter_one_transcript_per_gene(self):
        input_file = "test_data/sequence_ops/test_filter_one_transcript_per_gene/input.bed"
        expected_file = "test_data/sequence_ops/test_filter_one_transcript_per_gene/expected.bed"
        inputs = collections.defaultdict(lambda: [])
        [inputs[i[3]].append(i) for i in gen.read_many_fields(input_file, "\t")]
        expected = {i[3]: i for i in gen.read_many_fields(expected_file, "\t")}
        observed = filter_one_transcript_per_gene(inputs)
        self.assertEqual(expected, observed)

    def test_list_transcript_ids_from_features(self):
        input_file = "test_data/sequence_ops/test_list_transcript_ids_from_features/input.gtf"
        expected = ["ENST00000456328", "ENST00032323", "ENST0003246"]
        observed = list_transcript_ids_from_features(input_file)
        self.assertEqual(expected, observed)

    # for the Genome_Features class

    def test_build_cds(self):
        input_bed_file = "test_data/sequence_ops/test_Genome_Functions/test_build_cds/input.bed"
        input_fasta_file = "test_data/sequence_ops/test_Genome_Functions/test_build_cds/input.fasta"
        expected_file = "test_data/sequence_ops/test_Genome_Functions/test_build_cds/expected.fasta"
        observed_file = "test_data/sequence_ops/test_Genome_Functions/test_build_cds/observed.fasta"
        test = Genome_Functions(None, input_fasta_file)
        test.load_genome_dataset("test", input_bed_file)
        observed_seqs = test.build_cds(input_bed_file, observed_file)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_get_cds_features(self):
        input_file = "test_data/sequence_ops/test_Genome_Functions/test_get_cds_features/input.bed"
        expected_file = "test_data/sequence_ops/test_Genome_Functions/test_get_cds_features/expected.bed"
        inputs = gen.read_many_fields(input_file, "\t")
        test = Genome_Functions(None, None)
        test.load_genome_dataset("test", input_file)
        entries = gen.read_many_fields(expected_file, "\t")
        expected = collections.defaultdict(lambda: [])
        for i in entries:
            i[1], i[2] = int(i[1]), int(i[2])
            expected[i[3]].append(i)
        observed = test.get_cds_features()
        self.assertEqual(expected, observed)

    def test_get_stop_codon_features(self):
        input_file = "test_data/sequence_ops/test_Genome_Functions/test_get_stop_codon_features/input.bed"
        expected_file = "test_data/sequence_ops/test_Genome_Functions/test_get_stop_codon_features/expected.bed"
        inputs = gen.read_many_fields(input_file, "\t")
        test = Genome_Functions(None, None)
        test.load_genome_dataset("test", input_file)
        expected = {i[3]: [[i[0], int(i[1]), int(i[2]), i[3], i[4], i[5], i[6], i[7]]] for i in gen.read_many_fields(expected_file, "\t")}
        observed = test.get_stop_codon_features()
        self.assertEqual(expected, observed)

    def test_get_transcript_features(self):
        input_file = "test_data/sequence_ops/test_Genome_Functions/test_get_transcript_features/input.bed"
        expected_file = "test_data/sequence_ops/test_Genome_Functions/test_get_transcript_features/expected.bed"
        test = Genome_Functions(None, None)
        test.load_genome_dataset("test", input_file)
        expected = {i[3]: [[i[0], int(i[1]), int(i[2]), i[3], i[4], i[5], i[6], i[7]]] for i in gen.read_many_fields(expected_file, "\t")}
        observed = test.get_transcript_features()
        self.assertEqual(expected, observed)
