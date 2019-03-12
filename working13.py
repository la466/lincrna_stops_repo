import generic as gen
import seq_ops as seqo
import sequence_ops as sequo
import sim_ops_containers as soc
import SNP_ops as so
import main_tests_ops as mto
import ops
import numpy as np
import collections
import re
import os
import pandas as pd
# import scipy.stats as stats
# import matplotlib.pyplot as plt
import conservation
import os
import zipfile
from useful_motif_sets import stops, codon_map
from regex_patterns import codon_pattern
import containers as cont
import itertools as it
from Bio import SeqIO
import scipy.stats


output_directory = "clean_run/tests/clinvar"
gen.create_output_directories(output_directory)

clinvar_file = "source_data/clinvar_20190305.vcf"
exons_bed = "clean_run/genome_sequences/lincrna/cabili/multi_exons.bed"

disease_snp_intersect_file_vcf = "{0}/disease_snp_intersect.vcf".format(output_directory)
disease_snp_intersect_file_bed = "{0}/disease_snp_intersect.bed".format(output_directory)

print("Intersecting snps with exons")
so.intersect_snps_parallel(exons_bed, clinvar_file, disease_snp_intersect_file_vcf)
so.intersect_vcf_to_bed(exons_bed, disease_snp_intersect_file_vcf, disease_snp_intersect_file_bed, change_names = True)
