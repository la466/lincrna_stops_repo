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
    """
    Wrapper to extract clean CDS sequences.

    Args:
        gtf_file (str): path to gtf file containing genome features
        genome_file (str): path to fasta file containing genome sequence
        ortholog_gtf_file (str): path to ortholog gtf file containing ortholog genome features
        ortholog_genome_file (str): path to ortholog fasta file containing ortholog genome sequence
        orthologs_file (str): path to ensembl biomart file that contains focal-ortholog pairs
        output_directory (str): path to output directory
        clean_run (bool): if set, run
    """

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

    human_single_exon_bed = "{0}/genome_sequences/{1}/{1}.cds.single_exons.bed".format(output_directory, human_dataset_name)
    human_multi_exon_bed = "{0}/genome_sequences/{1}/{1}.cds.multi_exons.bed".format(output_directory, human_dataset_name)
    human_single_exon_fasta = "{0}/genome_sequences/{1}/{1}.cds.single_exons.fasta".format(output_directory, human_dataset_name)
    human_multi_exon_fasta = "{0}/genome_sequences/{1}/{1}.cds.multi_exons.fasta".format(output_directory, human_dataset_name)
    if not os.path.isfile(human_single_exon_bed) or not os.path.isfile(human_multi_exon_bed) or not os.path.isfile(human_single_exon_fasta) or not os.path.isfile(human_multi_exon_fasta) or clean_run:
        sequo.filter_by_exon_number(human_filelist["full_cds_features_bed"], human_cds_after_ortholog_filter_fasta, human_single_exon_bed, human_multi_exon_bed, human_single_exon_fasta, human_multi_exon_fasta)


def extract_exons(gtf_file, genome_file, output_directory, output_file, clean_run = None):

    print("Extracting exons...")

    human_dataset_name = "human"
    human_directory = "{0}/genome_sequences/{1}".format(output_directory, human_dataset_name)
    seqs_fasta = "{0}/genome_sequences/{1}/{1}.cds.conservation_filtered.step3.fasta".format(output_directory, human_dataset_name)
    if not os.path.isfile(output_file) or clean_run:
        # load the human dataset
        human, human_filelist = sequo.generate_genome_dataset(gtf_file, genome_file, human_dataset_name, human_directory, clean_run = False)
        # get the exons from the dataset
        ids_to_extract = gen.read_fasta(seqs_fasta)[0]
        exons = sequo.extract_gtf_feature(human_filelist["dataset_features_bed"], "exon", ids_to_keep = ids_to_extract)
        with open(output_file, "w") as outfile:
            for id in exons:
                [outfile.write("{0}\n".format("\t".join(i))) for i in exons[id]]
    


def extract_lincRNA_sequences(input_bed, genome_fasta, single_exon_bed, multi_exon_bed, single_exon_fasta, multi_exon_fasta, single_exon_families, multi_exon_families, clean_run = None):
    """
    Wrapper to extract lincRNA sequences

    Args:
        input_bed (str): path to bed file containing lincRNA coordinates
        genome_fasta (str): path to genome fasta file
        single_exon_bed (str): path to bed file containing coordinates of single exon cases
        multi_exon_bed (str): path to bed file containing coordinates of multi exon cases
        single_exon_fasta (str): path to fasta file containing sequences of single exon cases
        multi_exon_fasta (str): path to fasta file containing sequences of multi exon cases
        single_exon_families (str): path to bed file containing groupings of single exon paralogs
        multi_exon_families (str): path to bed file containing groupings of multi exon paralogs
        clean_run (bool): if set, run
    """

    print("Extracting lincRNA sequences...")

    if not os.path.isfile(single_exon_bed) or not os.path.isfile(multi_exon_bed) or not os.path.isfile(single_exon_fasta) or not os.path.isfile(multi_exon_fasta) or clean_run:
        # filter bed file to only containing those with strand information and then by number of exons
        sequo.filter_bed_file(input_bed, filter_columns = [5,9], filter_values = [["+", "-"], ["1"]], inclusive = [True, True], output_file = single_exon_bed)
        sequo.filter_bed_file(input_bed, filter_columns = [5,9], filter_values = [["+", "-"], ["1"]], inclusive = [True, False], output_file = multi_exon_bed)
        # have to convert to hg38
        fo.convert_bed(single_exon_bed)
        fo.convert_bed(multi_exon_bed)
        # now get sequence files
        # single exon cases is easier
        fo.fasta_from_intervals(single_exon_bed, single_exon_fasta, genome_fasta, names = True)
        # multi exon not so much, have to extract each piece by itself
        multi_exon_exons_bed = "{0}.exons.bed".format(".".join(multi_exon_bed.split(".")[:-1]))
        multi_exon_exons_fasta = "{0}.exons.fasta".format(".".join(multi_exon_fasta.split(".")[:-1]))
        sequo.extract_multi_exons_entries_to_bed(multi_exon_bed, multi_exon_exons_bed)
        fo.fasta_from_intervals(multi_exon_exons_bed, multi_exon_exons_fasta, genome_fasta, names = True)
        # build sequences from exons
        sequo.build_sequences_from_exon_fasta(multi_exon_exons_fasta, multi_exon_fasta)

    single_exon_blast_file = "{0}.blast_all_against_all.csv".format(".".join(single_exon_families.split(".")[:-1]))
    single_exon_blast_database = "{0}/single_exon_blast_all_against_all".format("/".join(single_exon_families.split("/")[:-1]))
    multi_exon_blast_file = "{0}.blast_all_against_all.csv".format(".".join(multi_exon_families.split(".")[:-1]))
    multi_exon_blast_database = "{0}/multi_exon_blast_all_against_all".format("/".join(multi_exon_families.split("/")[:-1]))
    # now group into paralagous families
    if not os.path.isfile(single_exon_families) or os.path.isfile(multi_exon_families) or clean_run:
        cons.filter_families(single_exon_fasta, single_exon_blast_file, single_exon_families, database_path = single_exon_blast_database, clean_run = clean_run)
        cons.filter_families(multi_exon_fasta, multi_exon_blast_file, multi_exon_families, database_path = multi_exon_blast_database, clean_run = clean_run)
