import generic as gen
import sequence_ops as sequo
import file_ops as fo
import time
import os

def extract_sequences(gtf_file, genome_file, ortholog_gtf_file, ortholog_genome_file, orthologs_file, output_directory, clean_run = None):

    # get a list of transcript features we want to keep
    transcript_ids = sequo.list_transcript_ids_from_features(gtf_file)

    human_dataset_name = "test"
    human_directory = "{0}/genome_sequences/{1}".format(output_directory, human_dataset_name)

    ortholog_dataset_name = "macaque"
    ortholog_directory = "{0}/genome_sequences/{1}".format(output_directory, ortholog_dataset_name)

    # generate the human dataset
    human_filelist = sequo.generate_genome_dataset(gtf_file, genome_file, transcript_ids, human_dataset_name, human_directory, clean_run = clean_run)

    # now create a file that contains all the gene orthologs for the orthologous species
    orthologous_pairs_file = "{0}/{1}.{2}.orthologs.bed".format(human_directory, human_dataset_name, ortholog_dataset_name)
    if not os.path.isfile(orthologous_pairs_file) or clean_run:
        sequo.get_orthologous_pairs(human_filelist["unique_transcript_gene_list_file"], orthologs_file, orthologous_pairs_file)
    #
    # # get a list of ortholog gene ids to keep
    # ortholog_ids = [i[1] for i in gen.read_many_fields(orthologous_pairs_file, "\t")]
