import generic as gen
import copy

args = ["python3", "main_tests.py", "clean_run", "clean_run", "-input_fasta", "clean_run/genome_sequences/human/human.cds.clean_coding_exons.fasta", "-input_fasta2", "clean_run/genome_sequences/human/human.clean_introns.fasta", "-families_file", "clean_run/genome_sequences/human/human.cds.families.bed", "-output_prefix", "pc", "--intron_length_test", "-ese_file"]
ese_files = [
    # "source_data/motif_sets/int3.txt",
    # "source_data/motif_sets/RESCUE.txt",
    "source_data/motif_sets/PESE.txt",
    "source_data/motif_sets/ESR.txt",
    "source_data/motif_sets/combined_eses.txt",
]


for file in ese_files:
    print("Running {0}".format(file))
    file_args = copy.deepcopy(args)
    file_args.append(file)

    print(" ".join(file_args))
    gen.run_process(file_args)
