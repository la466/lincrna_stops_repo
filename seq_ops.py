import generic as gen
import file_ops as fo
import ops
import numpy as np
import itertools as it
import re
import collections

def calc_seq_gc(seq):
    """
    Calculate GC content of a sequence

    Args:
        seq (str): input string

    Return:
        gc (float): gc content of string
    """
    gc = np.divide(sum([1 for i in list(seq) if i in ["G", "C"]]), len(seq))
    return gc


def get_non_transcribed_regions(input_gtf, input_fasta, output_features_bed, output_bed, output_fasta, output_directory):
    """
    Get the sequences that are not featured in a gtf file.

    Args:
        input_gtf (str): path to gtf file containing features
        input_fasta (str): path to genome fasta sequence
        output_features_bed (str): path to bed file containing the features
        output_bed (str): path to bed file containing regions not containing a feature
        output_fasta (str): path to fasta file containing sequences that dont contain a feature
        output_directory (str): path to the given output directory
    """

    # get all features in a bed file
    if not os.path.isfile(output_featured_bed):
        ops.extract_gtf_features_all(input_gtf, output_features_bed, exclude_XY=True)
    # create genome bed from fasta
    genome_bed = "{0}/genome.bed".format(output_directory)
    genome_index = "{0}.fai".format(input_fasta)
    if not os.path.isfile(genome_bed):
        ops.get_genome_bed_from_fasta_index(output_features_bed, genome_index, genome_bed)
    # # subtract the features from a bed file that simpy contains the whole genome coordinates
    # ops.intersect_bed(genome_bed, output_features_bed, output_file=output_bed, force_strand=True, subtract=True, no_dups=False, write_none=True)
    # # now get the fasta sequences of the remaining regions
    # fo.fasta_from_intervals(output_bed, output_fasta, input_fasta)
    #

def generate_nt_matched_seq(seq, dinucleotide_choices, dicnucleotide_probabilities, nucleotide_choices, nucleotide_probabilities, seed=None):
    """
    Generate a sequence based on dinucleotide and nucleotide frequencies.
    Args:
        seq (str): the sequence to simulate
        dinucleotide_choices (list): a list of sorted dinucleotides
        dicnucleotide_probabilities (list): a list of associated dinucleotide probablities
        nucleotide_choices (list): a list of sorted nucleotides
        nucleotide_probabilities (list): a list of associated nucleotide probabilities
        seed (bool): optionally set the randomisation seed
    Returns:
        sim_sequence (str): the simulated sequence
    """

    # optionally set seed
    if seed:
        np.random.seed(seed)

    required_dints = int(len(seq)/2)
    sim_content = list(np.random.choice(dinucleotide_choices, required_dints, p=dicnucleotide_probabilities))
    if len(seq) % 2 != 0:
        sim_content.append(np.random.choice(nucleotide_choices, p=nucleotide_probabilities))
    sim_sequence = "".join(sim_content)

    return sim_sequence


def generate_nt_matched_seq_mm(seq, nts, dinucleotide_probs, seed=None):
    """
    Generate a sequence based on dinucleotide and nucleotide frequencies.
    Args:
        seq (str): the sequence to simulate
        dinucleotide_choices (list): a list of sorted dinucleotides
        dicnucleotide_probabilities (list): a list of associated dinucleotide probablities
        nucleotide_choices (list): a list of sorted nucleotides
        nucleotide_probabilities (list): a list of associated nucleotide probabilities
        seed (bool): optionally set the randomisation seed
    Returns:
        sim_sequence (str): the simulated sequence
    """

    # optionally set seed
    if seed:
        np.random.seed(seed)

    sequence = seq[:2]
    while len(sequence) < len(seq):
        choice = np.random.choice(nts, p=dinucleotide_probs[seq[-2:]])
        sequence += choice
    return sequence


def generate_nt_matched_seq_frame(seq, dinucleotide_choices, dicnucleotide_probabilities, nucleotide_choices, nucleotide_probabilities, seq_frame=None, seed=None):
    """
    Generate a sequence based on dinucleotide and nucleotide frequencies.

    Args:
        seq (str): the sequence to simulate
        dinucleotide_choices (list): a list of sorted dinucleotides
        dicnucleotide_probabilities (list): a list of associated dinucleotide probablities
        nucleotide_choices (list): a list of sorted nucleotides
        nucleotide_probabilities (list): a list of associated nucleotide probabilities
        seq_frame (int): frame which the exon starts
        seed (bool): optionally set the randomisation seed

    Returns:
        sim_sequence (str): the simulated sequence
    """

    codon_regex = re.compile(".{3}")

    # optionally set seed
    if seed:
        np.random.seed(seed)

    generated = False
    failed = False
    attempts = 0
    while not generated:
        required_dints = int(len(seq)/2)
        sim_content = list(np.random.choice(dinucleotide_choices, required_dints, p=dicnucleotide_probabilities))
        if len(seq) % 2 != 0:
            sim_content.append(np.random.choice(nucleotide_choices, p=nucleotide_probabilities))
        sim_sequence = "".join(sim_content)

        codons = []
        codons.append(sim_sequence[:3-seq_frame])
        for i in range(0, len(sim_sequence), 3):
            codons.append(sim_sequence[i+3-seq_frame:i+6-seq_frame])

        stops_present = list(set(codons) & set(["TAA", "TAG", "TGA"]))
        if not len(stops_present):
            generated = True
        else:
            attempts += 1
            if attempts >= 200:
                failed = True
                generated = True



    return sim_sequence, failed


def get_dinucleotide_content(seqs, as_string=None):
    """
    Get the dinucleotide content of a list of seqs

    Args:
        seqs (list): list of sequences

    Returns:
        dinucleotide_content (dict): dict[dinucleotide] = dinucleotide proportion
    """

    dinucleotide_list = ["".join(i) for i in it.product("ACGT", repeat=2)]
    dinucleotide_count = collections.defaultdict(lambda: 0)

    if as_string:
        seqs = ["".join(seqs)]

    total_dints = 0
    for seq in seqs:
        dinucleotide_regex = re.compile('.{2}')
        dinucleotides1 = re.findall(dinucleotide_regex, seq)
        dinucleotides2 = re.findall(dinucleotide_regex, seq[1:])

        total_dints += len(dinucleotides1) + len(dinucleotides2)

        # get the counts
        for dint in dinucleotide_list:
            dinucleotide_count[dint] += dinucleotides1.count(dint)
            dinucleotide_count[dint] += dinucleotides2.count(dint)

    dinucleotide_content = {}
    # get proportions
    for i, count in dinucleotide_count.items():
        dinucleotide_content[i] = np.divide(count, total_dints)

    return dinucleotide_content


def get_nucleotide_content(seqs):
    """
    Get the nucleotide content of a list of seqs

    Args:
        seqs (list): list of sequences

    Returns:
        nucleotide_content (dict): dict[nucleotide] = nucleotide proportion
    """

    nts = ["A", "C", "G", "T"]
    seqs_string = "".join(seqs)
    seqs_nts = list(seqs_string)
    nucleotide_content = {}
    for nt in nts:
        nucleotide_content[nt] = np.divide(seqs_nts.count(nt), len(seqs_nts))
    return nucleotide_content


def get_dinucleotide_probabilities_markov(seqs):

    counts = collections.defaultdict(lambda: collections.defaultdict(lambda: 0))
    probabilities = collections.defaultdict(lambda: [])
    for seq in seqs:
        for i in range(len(seq)-2):
            counts[seq[i:i+2]][seq[i+2]] += 1
    for dint in counts:
        dint_counts = sum(counts[dint].values())
        for nt, count in sorted(counts[dint].items()):
            probabilities[dint].append(np.divide(count, dint_counts))
    return probabilities


def get_nucleotide_probabilities_markov(seqs):
    counts = collections.defaultdict(lambda: collections.defaultdict(lambda: 0))
    probabilities = collections.defaultdict(lambda: [])
    for seq in seqs:
        for i in range(len(seq)-1):
            counts[seq[i]][seq[i+1]] += 1
    for nt in counts:
        nt_counts = sum(counts[nt].values())
        for second_nt, count in sorted(counts[nt].items()):
            probabilities[nt].append(np.divide(count, nt_counts))
    return probabilities


def get_longest_orf(seq, strict_stop=None):
    """
    Get the longest ORF in a sequence

    Args:
        seq (str): the sequence

    Returns:
        longest (int): the longest ORF length in nucleotides
    """

    orfs = []
    stops = ["TAA", "TAG", "TGA"]
    reg = re.compile('^ATG(.{3})+?(?=(TAA|TAG|TGA))')

    start_re = re.compile("ATG")
    starts = start_re.finditer(seq)
    for start in starts:
        start_index = start.start()
        query_seq = seq[start_index:]
        orf_search = re.finditer(reg, query_seq)
        lengths = [len(match.group(0)) for match in orf_search]
        if len(lengths):
            orfs.append(lengths[0])
        else:
            # if we haven't specific that there needs to be a stop
            if not strict_stop:
                orfs.append(len(query_seq))

    # return the longest orf
    if len(orfs):
        return max(orfs)
    else:
        return np.nan


def get_longest_orfs(seqs):
    """
    Get the longest open reading frames in a list of sequences

    Args:
        seqs (list/dict): a list/dictionary of sequences

    Returns:
        longest_orfs (list/dict): A list/dictionary contaning then longest length reading frame for each seq
    """

    if isinstance(seqs, list):
        longest_orfs = []
        for seq in seqs:
            longest_orfs.append(get_longest_orf(seq))
        return longest_orfs

    if isinstance(seqs, dict):
        longest_orfs = {}
        for seq in seqs:
            longest_orfs[seq] = get_longest_orf(seqs[seq])
        return longest_orfs


def get_stop_count(seq):
    """
    Get the number of stop codons in any reading frame for a sequence.

    Args:
        seq (str): given sequence

    Returns:
        stop_count (int): number of stop codons found in any reading frame
    """

    stop_regex = re.compile('(TAA|TAG|TGA)')
    return len(re.findall(stop_regex, seq))


def get_stop_counts(seqs):
    """
    Count the number of stop codons in any reading frame.

    Args:
        seqs (list/dict): list or dictionary holding the sequences

    Returns:
        stop_counts (list/dict): list or dictionary with the stop counts.
    """

    if isinstance(seqs, list):
        return [get_stop_count(seq) for seq in seqs]

    if isinstance(seqs, dict):
        stop_counts = {}
        for id, seq in seqs.items():
            stop_counts[id] = get_stop_count(seq)
        return stop_counts


def get_unique_seqs(names, seqs_list):
    """
    From a list of fasta entries, return only unique sequences.

    Args:
        names (list): list of sequence identifiers
        seqs_lists (list): list of corresponding seqs

    Returns:
        seqs (dict): dict[name] = seq of unique sequences
    """

    unique_seqs = []
    seqs = {}
    for i, name in enumerate(names):
        if seqs_list[i] not in unique_seqs:
            unique_seqs.append(seqs_list[i])
            seqs[name] = seqs_list[i]
    return seqs


def fasta_to_dict(fasta_path):
    """
    Take a fasta file and return a dictionary of seqs with the
    fasta info as a key

    Args:
        fasta_path (str): path to fasta file

    Returns:
        out_dict (dict): dict[fasta_name] = fasta_seq
    """

    names, seqs = gen.read_fasta(fasta_path)
    out_dict = {}
    for i, name in enumerate(names):
        out_dict[name] = seqs[i]
    return out_dict
