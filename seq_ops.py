import generic as gen
import re

def get_shortest_orf(seq, strict_stop=None):
    """
    Get the shortest ORF in a sequence

    Args:
        seq (str): the sequence

    Returns:
        shortest_orf (int): the shortest ORF length in nucleotides
    """
    orfs = []
    stops = ["TAA", "TAG", "TGA"]
    start_re = re.compile("ATG")
    codon_re = re.compile(".{3}")
    starts = start_re.finditer(seq)
    for start in starts:
        start_index = start.start()
        query_seq = seq[start_index:]
        query_codons = codon_re.findall(query_seq)
        stops_exist = [stop for stop in stops if stop in query_codons]
        # if a stop codon exists
        if len(stops_exist):
            for i, codon in enumerate(query_codons):
                if codon in stops:
                    # only cound coding codons, i.e not the stop codon
                    nts = i*3
                    #only want cases where there is something more than a START:STOP
                    if nts > 3:
                        orfs.append(nts)
        # or if the sequence doesnt have a stop, use the full length of the rest of the sequence
        else:
            # if we haven't specified a strict stop
            if not strict_stop:
                orfs.append(len(query_seq))
    # return the shortest orf
    if len(orfs):
        return min(orfs)
    else:
        return "na"


def get_shortest_orfs(seqs):
    """
    Get the shortest open reading frames in a list of sequences

    Args:
        seqs (list): a list of sequences

    Returns:
        shortest_orfs (list): A list contaning then shortest length reading frame for each seq
    """

    shortest_orfs = []
    for seq in seqs:
        shortest_orfs.append(get_shortest_orf(seq))
    return shortest_orfs


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
