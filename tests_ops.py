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

    def test_check_exon_files_fail(self):
        input_file1 = "test_data/ops/test_check_exon_files_fail/input1.bed"
        input_file2 = "test_data/ops/test_check_exon_files_fail/input2.bed"
        self.assertRaises(Exception, check_exon_files, [input_file1, input_file2])

    def test_check_exon_files_pass(self):
        input_file1 = "test_data/ops/test_check_exon_files_pass/input1.bed"
        input_file2 = "test_data/ops/test_check_exon_files_pass/input2.bed"
        observed = check_exon_files(input_file1, input_file2)
        self.assertTrue(observed)

    def test_count_stop_codons(self):
        input_file = "test_data/ops/test_count_stop_codons/input.fasta"
        seqs = gen.read_fasta(input_file)[1]
        expected = [0,1,4]
        observed = [count_stop_codons(seq) for seq in seqs]
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
        filter_bed_from_fasta(input_bed, input_fasta, observed_file)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_filter_bed_transcript_id(self):
        input_bed1 = "test_data/ops/test_filter_bed_transcript_id/input1.bed"
        input_bed2 = "test_data/ops/test_filter_bed_transcript_id/input2.bed"
        expected_file = "test_data/ops/test_filter_bed_transcript_id/expected.bed"
        observed_file = "test_data/ops/test_filter_bed_transcript_id/observed.bed"
        gen.remove_file(observed_file)
        filter_bed_transcript_id(input_bed1, input_bed2, observed_file)
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

    def test_get_coding_exons(self):
        exon_file = "test_data/ops/test_get_coding_exons/exons.bed"
        CDS_file = "test_data/ops/test_get_coding_exons/CDSs.bed"
        expected_file = "test_data/ops/test_get_coding_exons/expected.bed"
        observed_file = "test_data/ops/test_get_coding_exons/observed.bed"
        gen.remove_file(observed_file)
        get_coding_exons(exon_file, CDS_file, observed_file)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_get_exon_reading_frame(self):
        input_full_exons_file = "test_data/ops/test_get_exon_reading_frame/input_full_exons.fasta"
        input_coding_exons_file = "test_data/ops/test_get_exon_reading_frame/input_coding_exons.fasta"
        expected_file = "test_data/ops/test_get_exon_reading_frame/expected.fasta"
        observed_file = "test_data/ops/test_get_exon_reading_frame/observed.fasta"
        gen.remove_file(observed_file)
        get_exon_reading_frame(input_coding_exons_file, input_full_exons_file, observed_file)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_get_region_stop_counts(self):
        input_file = "test_data/ops/test_get_region_stop_counts/input.fasta"
        input_names, input_seqs = gen.read_fasta(input_file)
        region_start, region_end = 3, 69
        input_list = {name: input_seqs[i] for i, name in enumerate(input_names) if len(input_seqs[i]) > region_end}
        expected = [33,7]
        observed = get_region_stop_counts(input_list, region_start, region_end)
        self.assertEqual(expected, observed)

    def test_link_transcripts_to_genes(self):
        input_file = "test_data/ops/test_link_transcripts_to_genes/input.bed"
        expected = {"ENSG1": ["ENST100", "ENST102"], "ENSG2": ["ENST2"], "ENSG3": ["ENST106"]}
        observed = link_transcripts_to_genes(input_file)
        self.assertEqual(expected, observed)

    def test_remove_bed_intersects(self):
        input1 = "test_data/ops/test_remove_bed_intersects/input1.bed"
        input2 = "test_data/ops/test_remove_bed_intersects/input2.bed"
        expected_file = "test_data/ops/test_remove_bed_intersects/expected.bed"
        observed_file = "test_data/ops/test_remove_bed_intersects/observed.bed"
        gen.remove_file(observed_file)
        remove_bed_intersects(input1, input2, observed_file)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_remove_bed_overlaps(self):
        input_file = "test_data/ops/test_remove_bed_overlaps/input.bed"
        expected_file = "test_data/ops/test_remove_bed_overlaps/expected.bed"
        observed_file = "test_data/ops/test_remove_bed_overlaps/observed.bed"
        gen.remove_file(observed_file)
        remove_bed_overlaps(input_file, observed_file)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_remove_bed_overlaps_2(self):
        input_file = "test_data/ops/test_remove_bed_overlaps_2/input.bed"
        expected_file = "test_data/ops/test_remove_bed_overlaps_2/expected.bed"
        observed_file = "test_data/ops/test_remove_bed_overlaps_2/observed.bed"
        gen.remove_file(observed_file)
        remove_bed_overlaps(input_file, observed_file)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_sort_bed(self):
        input_file = "test_data/ops/test_sort_bed/input.bed"
        expected_file = "test_data/ops/test_sort_bed/expected.bed"
        observed_file = "test_data/ops/test_sort_bed/observed.bed"
        gen.remove_file(observed_file)
        sort_bed(input_file, observed_file)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_uniquify_lincRNA_transcripts(self):
        input_file = "test_data/ops/test_uniquify_lincRNA_transcripts/input.fasta"
        map_file = "test_data/ops/test_uniquify_lincRNA_transcripts/map.txt"
        expected_file = "test_data/ops/test_uniquify_lincRNA_transcripts/expected.fasta"
        observed_file = "test_data/ops/test_uniquify_lincRNA_transcripts/observed.fasta"
        gen.remove_file(observed_file)
        uniquify_lincRNA_transcripts(input_file, map_file, observed_file)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
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
