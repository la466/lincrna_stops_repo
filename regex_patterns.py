import re

"""
Set of regex patterns
"""

# get codons, that is any set of 3 characters (assumes pre filtering)
codon_pattern = re.compile(".{3}")
# get the exon number
exon_number_pattern = re.compile("(?<=exon_number)\s\"(\d+)\"")
# get the gene id
gene_id_pattern = re.compile("ENS\w*G\d+")
# get the transcript id
transcript_id_pattern = re.compile("ENS\w*T\d*")
# get stop codons
stop_codon_pattern = re.compile("(TAA|TAG|TGA)")
