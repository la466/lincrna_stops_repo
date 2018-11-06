import generic as gen
import sequence_ops as sequo
import file_ops as fo
import time

def extract_sequences(gtf_file, genome_file, output_directory, clean_run = None):

    # get a list of transcript features we want to keep
    transript_ids = sequo.list_transcript_ids_from_features(gtf_file)
    # create the genome features directory if it doesnt already exist
    genome_features_directory = "{0}/genome_features".format(output_directory)
    if clean_run:
        gen.remove_directory(genome_features_directory)

    gen.create_output_directories(genome_features_directory)
    # get the genome features
    human_genome = sequo.Genome_Functions(gtf_file, genome_file)
    # create the features dataset
    # we only want to do this once if possible, or if we want a clean run
    human_dataset_name = "human"
    human_dataset_features_bed = "{0}/{1}.{2}_features_dataset.bed".format(genome_features_directory, gtf_file.split("/")[-1], human_dataset_name)
    human_genome.generate_genome_dataset(human_dataset_name, human_dataset_features_bed, input_list = transript_ids, clean_run = clean_run)
    # load the dataset
    human_genome.load_genome_dataset(human_dataset_name, human_dataset_features_bed)

    # get the CDS features
    cds_features = human_genome.get_cds_features()
    # write the cds_features to an output_file
    full_cds_features_bed = "{0}/{1}.cds.full_features.bed".format(genome_features_directory, human_dataset_name)
    sequo.write_features_to_bed(cds_features, full_cds_features_bed, clean_run = clean_run)

    # build the cds sequences
    full_cds_fasta = "{0}/{1}.cds.no_filtering_full.fasta".format(genome_features_directory, human_dataset_name)
    human_genome.build_cds(full_cds_features_bed, full_cds_fasta)

    # filter the sequences
    quality_filtered_cds_fasta = "{0}/{1}.cds.quality_filtered.step1.fasta".format(genome_features_directory, human_dataset_name)
    sequo.filter_cds(full_cds_fasta, quality_filtered_cds_fasta)

    clean_transcript_ids, clean_transcript_seqs = gen.read_fasta(quality_filtered_cds_fasta)
    # get a list of transcripts from the gtf file and keep only those that we want
    transcript_list = human_genome.get_transcript_features()
    transcript_list = {i: transcript_list[i] for i in transcript_list if i in clean_transcript_ids}
    # keep only the longest transcript if more than one per gene
    unique_gene_transcripts = sequo.filter_one_transcript_per_gene(transcript_list)
    transcript_list = unique_gene_transcripts
    unique_gene_cds = {name: clean_transcript_seqs[i] for i, name in enumerate(clean_transcript_ids) if name in transcript_list}
    # write the unique transcripts to file
    unique_cds_fasta = "{0}/{1}.cds.unique_gene_transcripts.step2.fasta".format(genome_features_directory, human_dataset_name)
    fo.write_fasta(unique_gene_cds, unique_cds_fasta)
