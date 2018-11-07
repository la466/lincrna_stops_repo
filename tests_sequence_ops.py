from sequence_ops import *
import numpy as np
import unittest
import collections

class Tests(unittest.TestCase):

    def test_build_coding_sequences(self):
        input_bed_file = "test_data/sequence_ops/test_build_coding_sequences/input.bed"
        input_fasta_file = "test_data/sequence_ops/test_build_coding_sequences/input.fasta"
        expected_file = "test_data/sequence_ops/test_build_coding_sequences/expected.fasta"
        observed_file = "test_data/sequence_ops/test_build_coding_sequences/observed.fasta"
        gen.remove_file(observed_file)
        build_coding_sequences(input_bed_file, input_fasta_file, observed_file)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_extract_cds_features(self):
        input_file = "test_data/sequence_ops/test_extract_cds_features/input.bed"
        expected_file = "test_data/sequence_ops/test_extract_cds_features/expected.bed"
        input_list = ["ENST00000620200", "ENST00000433179"]
        entries = gen.read_many_fields(expected_file, "\t")
        expected = collections.defaultdict(lambda: [])
        [expected[i[3]].append(i) for i in entries]
        observed = extract_cds_features(input_file, input_list)
        self.assertEqual(expected, observed)

    def test_extract_gtf_features(self):
        input_file = "test_data/sequence_ops/test_extract_gtf_features/input.gtf"
        expected_file = "test_data/sequence_ops/test_extract_gtf_features/expected.bed"
        input_list = ["ENST00000456328", "ENST0003246"]
        expected = gen.read_many_fields(expected_file, "\t")
        observed = extract_gtf_features(input_list, input_file, filter_transcripts = True)
        self.assertEqual(expected, observed)

    def test_extract_stop_codon_features(self):
        input_file = "test_data/sequence_ops/test_extract_stop_codon_features/input.bed"
        expected_file = "test_data/sequence_ops/test_extract_stop_codon_features/expected.bed"
        input_list = ["ENST00000620200", "ENST00000433179"]
        inputs = gen.read_many_fields(input_file, "\t")
        expected = {i[3]: [i] for i in gen.read_many_fields(expected_file, "\t")}
        observed = extract_stop_codon_features(inputs, input_list)
        self.assertEqual(expected, observed)

    def test_extract_transcript_features(self):
        input_file = "test_data/sequence_ops/test_extract_transcript_features/input.bed"
        expected_file = "test_data/sequence_ops/test_extract_transcript_features/expected.bed"
        input_list = ["ENST00000620200", "ENST00000433179"]
        inputs = gen.read_many_fields(input_file, "\t")
        expected = {i[3]: [i] for i in gen.read_many_fields(expected_file, "\t")}
        observed = extract_transcript_features(inputs, input_list)
        self.assertEqual(expected, observed)

    def test_filter_one_transcript_per_gene(self):
        input_file = "test_data/sequence_ops/test_filter_one_transcript_per_gene/input.bed"
        expected_file = "test_data/sequence_ops/test_filter_one_transcript_per_gene/expected.bed"
        inputs = collections.defaultdict(lambda: [])
        [inputs[i[3]].append(i) for i in gen.read_many_fields(input_file, "\t")]
        expected = {i[3]: i for i in gen.read_many_fields(expected_file, "\t")}
        observed = filter_one_transcript_per_gene(inputs)
        self.assertEqual(expected, observed)

    def test_get_clean_genes(self):
        input_bed_file = "test_data/sequence_ops/test_get_clean_genes/input.bed"
        input_fasta_file = "test_data/sequence_ops/test_get_clean_genes/input.fasta"
        expected_file = "test_data/sequence_ops/test_get_clean_genes/expected.bed"
        observed_file = "test_data/sequence_ops/test_get_clean_genes/observed.bed"
        gen.remove_file(observed_file)
        get_clean_genes(input_fasta_file, input_bed_file, observed_file)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_get_orthologous_pairs(self):
        input_bed_file = "test_data/sequence_ops/test_get_orthologous_pairs/input.bed"
        input_pairs_file = "test_data/sequence_ops/test_get_orthologous_pairs/input.csv"
        expected_file = "test_data/sequence_ops/test_get_orthologous_pairs/expected.bed"
        observed_file = "test_data/sequence_ops/test_get_orthologous_pairs/observed.bed"
        gen.remove_file(observed_file)
        get_orthologous_pairs(input_bed_file, input_pairs_file, observed_file)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_list_transcript_ids_from_features(self):
        input_file = "test_data/sequence_ops/test_list_transcript_ids_from_features/input.gtf"
        expected = ["ENST00000456328", "ENST00032323", "ENST0003246"]
        observed = list_transcript_ids_from_features(input_file)
        self.assertEqual(expected, observed)

    def test_quality_filter_cds_sequences(self):
        input_file = "test_data/sequence_ops/test_quality_filter_cds_sequences/input.fasta"
        expected_file = "test_data/sequence_ops/test_quality_filter_cds_sequences/expected.fasta"
        observed_file = "test_data/sequence_ops/test_quality_filter_cds_sequences/observed.fasta"
        gen.remove_file(observed_file)
        quality_filter_cds_sequences(input_file, observed_file)
        # use an arbritary delimiter
        expected = gen.read_many_fields(expected_file, "#")
        observed = gen.read_many_fields(observed_file, "#")
        self.assertEqual(expected, observed)
