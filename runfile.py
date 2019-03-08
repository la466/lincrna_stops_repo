import generic as gen

arg_list = {
    "lincrna_dints": ["python3", "lincRNA.py", "None", "source_data/cabili_clean_alignments.fasta", "clean_run/tests/lincrna/substitution_rates", "--substitution_rate", "-motif_file", "source_data/motif_sets/int3.txt", "-families_file", "clean_run/genome_sequences/lincrna/cabili/multi_exon_families.txt", "-required_simulations", 1000, "-output_prefix", "lincrna"],
    "pc_dints": ["python3", "lincRNA.py", "None", "clean_run/genome_sequences/human.macaque.alignments.fasta", "clean_run/tests/lincrna/substitution_rates", "--substitution_rate", "-motif_file", "source_data/motif_sets/int3.txt", "-families_file", "clean_run/genome_sequences/human/human.cds.families.bed", "-required_simulations", 1000, "-output_prefix", "pc"],
}

for arg in arg_list:
    print("Running {0}".format(arg))
    gen.run_process(arg_list[arg])
