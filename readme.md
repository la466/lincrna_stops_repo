DISCLAIMER The code in this repository has not been packaged to be run out of the box. For instance, there are software dependencies that haven't been explicitly documented, and the documentation regarding how to run the scripts is also insufficiently documented for external use. The purpose of the repository is to serve as an explicit record of the analyses performed in the paper. It is thus primarily a supplement to the methods. Some functions in some files may not be used. Credit to Rosina Savisaar for some functions that have been transferred from other projects.

----------

Source code for the Abrahams and Hurst manuscript in preparation May 2019.

test_data/: contains test cases for functions.

conservation.py: contains functions for comparing between sequences, e.g. blast, ds

containers.py: container functions for genome sequence operations

extract_sequences.py: contains functions to extract sequences

file_ops.py: contains functions for working with files

generic.py: contains generic utility unctions

lincRNA_misc_ops.py: contains functions that were generally used to test lincRNA

lincRNA_misc.py: contains functions mainly to further extract lincRNA sequences

lincRNA_tests_ops.py: contains functions for lincRNA tests

lincRNA.py: main wrapper for lincRNA tests

main_tests_ops.py: contains functions used for other tests

main_tests.py: main wrapper for other tests

motif_tests_ops.py: contains functions used on motif sets

motif_tests.py: main wrapper for tests on motifs

ops_containers.py: contains functions mainly wrapper other ops

ops.py: contains functions used to extract sequences

R_ese_density_intron_length.R: intron length and ESE density tests

R_intron_hexamers.R: intron hexamer tests

R_lincrna_orf_lengths.R: intron orf length tests

R_motif_codon_nd.R: tests for codon sets in motifs

R_motif_overlap_capability.R: tests looking at ability for motifs to overlap

R_motif_sets_stop_densities.R: stop densities in motif sets

R_purine_content.R: purine content tests

regex_patterns.py: utility regex pattens file

seq_ops.py: contains functions for sequence operations

sequence_ops.py: contains more functions on sequences operations

sim_ops_containers.py: contains functions normally used to wrap simulation functions

sim_ops: contains functions used for simulations

tests_*.py: unittests for functions in respective files

useful_motif_sets.py: contains some useful motifs and lists of things
