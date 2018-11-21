import generic as gen
import sequence_ops as sequo
import conservation as cons
import file_ops as fo
import time
import os
import collections

def check_conservation(human_cds_fasta, ortholog_cds_fasta, ortholog_transcripts_links_file, output_file, output_fasta, max_dS_threshold = None, max_omega_threshold = None, clean_run = None):
    """
    Wrapper to check conservation between cds and orthologous cds. Only keep transcripts
    that have dS less than max_dS_threshold and omega less than max_omega_threshold.

    Args:
        human_cds_fasta (str): path to fasta containing human cds
        ortholog_cds_fasta (str): path to fasta containing ortholog cds
        ortholog_transcripts_links_file (str): path to bed file containing the links between transcript ids and orthologs
        output_file (str): path to file containing the retained transcript ids
        output_fasta (str): path to fasta file containing the retained sequences
        max_dS_threshold (float): if set, the max that at least one ortholog must have dS below
        max_omega_threshold (float): if set, the max that at least one ortholog must have omega below
    """
    # first generate a dictionary containing each of the transcripts with the orthologous transcripts
    transcripts_with_orthologs = sequo.get_transcript_and_orthologs(human_cds_fasta, ortholog_cds_fasta, ortholog_transcripts_links_file)
    # remove lambda for pickling
    transcripts_with_orthologs = {i: transcripts_with_orthologs[i] for i in transcripts_with_orthologs}
    if not os.path.isfile(output_file) or clean_run:
        cons.get_conservation(transcripts_with_orthologs, output_file, max_dS_threshold = max_dS_threshold, max_omega_threshold = max_omega_threshold)
    fo.filter_fasta_from_bed(output_file, human_cds_fasta, output_fasta, filter_column = 0)


def extract_clean_sequences(gtf_file, genome_file, ortholog_gtf_file, ortholog_genome_file, orthologs_file, output_directory, clean_run = None):

    # # get a list of transcript features we want to keep
    transcript_ids = sequo.list_transcript_ids_from_features(gtf_file)

    human_dataset_name = "human"
    human_directory = "{0}/genome_sequences/{1}".format(output_directory, human_dataset_name)

    ortholog_dataset_name = "macaque"
    ortholog_directory = "{0}/genome_sequences/{1}".format(output_directory, ortholog_dataset_name)

    # generate the human dataset
    human, human_filelist = sequo.generate_genome_dataset(gtf_file, genome_file, human_dataset_name, human_directory, id_list = transcript_ids, filter_by_transcript = True, filter_one_per_gene = True, clean_run = clean_run)

    # now create a file that contains all the gene orthologs for the orthologous species
    orthologous_pairs_file = "{0}/{1}.{2}.orthologs.bed".format(human_directory, human_dataset_name, ortholog_dataset_name)
    if not os.path.isfile(orthologous_pairs_file) or clean_run:
        sequo.get_orthologous_pairs(human_filelist["unique_transcript_gene_list_file"], orthologs_file, orthologous_pairs_file)

    # get a list of unique ortholog gene ids to keep
    ortholog_gene_ids = list(set([i[1] for i in gen.read_many_fields(orthologous_pairs_file, "\t")]))

    # now get the ortholog sequences
    ortholog, ortholog_filelist = sequo.generate_genome_dataset(ortholog_gtf_file, ortholog_genome_file, ortholog_dataset_name, ortholog_directory, id_list = ortholog_gene_ids, filter_by_gene = True, clean_run = clean_run)

    # retrieve a list of transcripts and their ortholog transcript links,
    # leaving only those orthologs that appear in the cds fasta (some may have been filtered)
    orthologs_transcript_links_file = "{0}/genome_sequences/{1}.{2}.transcripts.ortholog_links.bed".format(output_directory, human_dataset_name, ortholog_dataset_name)
    if not os.path.isfile(orthologs_transcript_links_file) or clean_run:
        sequo.get_ortholog_transcript_pairings(human_filelist["unique_transcript_gene_list_file"], ortholog_filelist["unique_transcript_gene_list_file"], orthologous_pairs_file, ortholog.cds_fasta, orthologs_transcript_links_file)

    # check orthologs conservation
    human_ids_after_conservation_file = "{0}/genome_sequences/{1}.{2}.conservation_filtering.step3.bed".format(output_directory, human_dataset_name, ortholog_dataset_name)
    human_cds_after_ortholog_filter_fasta = "{0}/genome_sequences/{1}/{1}.cds.conservation_filtered.step3.fasta".format(output_directory, human_dataset_name)
    if not os.path.isfile(human_ids_after_conservation_file) or not os.path.isfile(human_cds_after_ortholog_filter_fasta) or clean_run:
        check_conservation(human.cds_fasta, ortholog.cds_fasta, orthologs_transcript_links_file, human_ids_after_conservation_file, human_cds_after_ortholog_filter_fasta, max_dS_threshold = 0.2, max_omega_threshold = 0.5, clean_run = clean_run)

    # filter the sequences into families
    human_blast_file = "{0}/genome_sequences/{1}/{1}.cds.blast_all_against_all.csv".format(output_directory, human_dataset_name)
    human_families_file = "{0}/genome_sequences/{1}/{1}.cds.families.bed".format(output_directory, human_dataset_name)
    human_blast_database_path = "{0}/genome_sequences/{1}/blast_all_against_all".format(output_directory, human_dataset_name)
    if not os.path.isfile(human_blast_file) or not os.path.isfile(human_families_file) or clean_run:
        cons.filter_families(human_cds_after_ortholog_filter_fasta, human_blast_file, human_families_file, database_path = human_blast_database_path, clean_run = clean_run)
