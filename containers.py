import generic as gen
import sequence_ops as sequo
import conservation as cons
import file_ops as fo
import time
import os
import collections

def extract_sequences(gtf_file, genome_file, ortholog_gtf_file, ortholog_genome_file, orthologs_file, output_directory, clean_run = None):

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

    # # check orthologs conservation
    # check_conservation(human.cds_fasta, ortholog.cds_fasta, orthologs_transcript_links_file)



def check_conservation(human_cds_fasta, ortholog_cds_fasta, ortholog_transcripts_links, output_directory, clean_run = None):

    # first generate a dictionary containing each of the transcripts with the orthologous transcripts
    transcripts_with_orthologs = sequo.get_transcript_and_orthologs(human_cds_fasta, ortholog_cds_fasta, ortholog_transcripts_links)

    cons.get_conservation(transcripts_with_orthologs)
