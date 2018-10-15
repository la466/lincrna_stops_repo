from ops import *
import numpy as np
import unittest

class Test_Ops(unittest.TestCase):

    def test_build_sequences(self):
        input_file = "test_data/ops/test_build_sequences/input.fasta"
        input_stops_file = "test_data/ops/test_build_sequences/input_stops.fasta"
        expected_file = "test_data/ops/test_build_sequences/expected.fasta"
        observed_file = "test_data/ops/test_build_sequences/observed.fasta"
        gen.remove_file(observed_file)
        build_sequences(input_file, input_stops_file, observed_file)
        #use an arbritary seperator
        expected = gen.read_many_fields(expected_file, "#")
        observed = gen.read_many_fields(observed_file, "#")
        self.assertEqual(expected, observed)

    def test_extract_features_cds(self):
        input_file = "test_data/ops/test_extract_features_cds/input.gtf"
        expected_file = "test_data/ops/test_extract_features_cds/expected.bed"
        observed_file = "test_data/ops/test_extract_features_cds/observed.bed"
        gen.remove_file(observed_file)
        extract_features(input_file, ['CDS'], observed_file)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_extract_features_cds_stop_codon(self):
        input_file = "test_data/ops/test_extract_features_cds_stop_codon/input.gtf"
        expected_file = "test_data/ops/test_extract_features_cds_stop_codon/expected.bed"
        observed_file = "test_data/ops/test_extract_features_cds_stop_codon/observed.bed"
        gen.remove_file(observed_file)
        extract_features(input_file, ['CDS', 'stop_codon'], observed_file)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_extract_features_stop_codon(self):
        input_file = "test_data/ops/test_extract_features_stop_codon/input.gtf"
        expected_file = "test_data/ops/test_extract_features_stop_codon/expected.bed"
        observed_file = "test_data/ops/test_extract_features_stop_codon/observed.bed"
        gen.remove_file(observed_file)
        extract_features(input_file, ['stop_codon'], observed_file)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_extract_gtf_features_exon(self):
        gtf_file = "test_data/ops/test_extract_gtf_features_exon/input.gtf"
        expected_file = "test_data/ops/test_extract_gtf_features_exon/expected.bed"
        observed_file = "test_data/ops/test_extract_gtf_features_exon/observed.bed"
        gen.remove_file(observed_file)
        extract_gtf_features(gtf_file, ['exon'], observed_file)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_extract_gtf_features_exon_intron(self):
        gtf_file = "test_data/ops/test_extract_gtf_features_exon_intron/input.gtf"
        expected_file = "test_data/ops/test_extract_gtf_features_exon_intron/expected.bed"
        observed_file = "test_data/ops/test_extract_gtf_features_exon_intron/observed.bed"
        gen.remove_file(observed_file)
        extract_gtf_features(gtf_file, ['exon', 'intron'], observed_file)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_extract_gtf_features_stop_codon(self):
        gtf_file = "test_data/ops/test_extract_gtf_features_stop_codon/input.gtf"
        expected_file = "test_data/ops/test_extract_gtf_features_stop_codon/expected.bed"
        observed_file = "test_data/ops/test_extract_gtf_features_stop_codon/observed.bed"
        gen.remove_file(observed_file)
        extract_gtf_features(gtf_file, ['stop_codon'], observed_file)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_filter_bed_from_fasta(self):
        input_bed = "test_data/ops/test_filter_bed_from_fasta/input.bed"
        input_fasta = "test_data/ops/test_filter_bed_from_fasta/input.fasta"
        expected_file = "test_data/ops/test_filter_bed_from_fasta/expected.bed"
        observed_file = "test_data/ops/test_filter_bed_from_fasta/observed.bed"
        gen.remove_file(observed_file)
        filter_bed_from_filtered_cds(input_bed, input_fasta, observed_file)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_filter_coding_sequences(self):
        input_file = "test_data/ops/test_filter_coding_sequences/input.fasta"
        expected_file = "test_data/ops/test_filter_coding_sequences/expected.fasta"
        observed_file = "test_data/ops/test_filter_coding_sequences/observed.fasta"
        gen.remove_file(observed_file)
        filter_coding_sequences(input_file, observed_file)
        # use an arbritary delimiter
        expected = gen.read_many_fields(expected_file, "#")
        observed = gen.read_many_fields(observed_file, "#")
        self.assertEqual(expected, observed)

    def test_link_transcripts_to_genes(self):
        input_file = "test_data/ops/test_link_transcripts_to_genes/input.bed"
        expected = {"ENSG1": ["ENST100", "ENST102"], "ENSG2": ["ENST2"], "ENSG3": ["ENST106"]}
        observed = link_transcripts_to_genes(input_file)
        self.assertEqual(expected, observed)

    def test_uniquify_transcripts(self):
        input_file = "test_data/ops/test_uniquify_transcripts/input.fasta"
        expected_file = "test_data/ops/test_uniquify_transcripts/expected.fasta"
        observed_file = "test_data/ops/test_uniquify_transcripts/observed.fasta"
        gen.remove_file(observed_file)
        input_links = {"ENSG1": ["ENST0001", "ENST0004", "ENST0005"], "ENSG2": ["ENST0003"], "ENSG3": ["ENST0007", "ENST0008"], "ENSG4": ["ENST0002", "ENST0006", "ENST0009", "ENST00010"]}
        uniquify_transcripts(input_file, input_links, observed_file)
        # use an arbritary delimiter
        expected = gen.read_many_fields(expected_file, "#")
        observed = gen.read_many_fields(observed_file, "#")
        self.assertEqual(expected, observed)
