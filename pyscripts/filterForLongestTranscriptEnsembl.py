#!/bin/env python3

import argparse
import os
import textwrap
import gzip

"""
filter a protein fasta file for the longest isoform
"""

class SEQ:
    def __init__(self, seq: str, header: str):
        self.seq    = seq
        self.header = header
        self.length = len(seq)

def get_arguments() -> str:
    """get the arguments"""

    d = "Some python code to filter protein sequences to the longest isoform from an ensembl peptide fasta"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-f", "--fasta", help="A *.pep.all.fa file", required=True)

    args  = parser.parse_args()
    fasta = args.fasta

    assert os.path.isfile(fasta), f"Could not locate file: {fasta}"
    
    return fasta

def get_gene_id(header: str) -> str:
    """return the gene id from an ensembl header"""

    sects = header.split(' ')
    gene  = sects[3]
    gene  = gene.split(':')[1]
    gene  = gene.split('.')[0]
    return gene

def write_fasta(seqs: dict, fasta: str) -> int:
    """write the Seqs dictionary to a fasta file"""

    path  = os.path.dirname(fasta)
    fname = os.path.basename(fasta)
    fname = f"{path}/Filtered.{fname}"

    if (fname.endswith(".gz")):
        fh = gzip.open(fname, "wt")
    else:
        fh = open(fname, 'w')

    for seq in seqs.values():
        fh.write(seq.header + '\n')
        fh.write(textwrap.fill(seq.seq, width=60))
        fh.write('\n')

    fh.close()

    return fasta


def filter_fasta(fasta: str) -> int:
    """filter a fasta file for a single representative isoform per gene"""

    seqs = dict()
    fh   = gzip.open(fasta, "rt") if fasta.endswith(".gz") else open(fasta, 'r')
    cur  = ''
    trn  = ''
    seq  = ''

    for line in fh:
        if (len(line) == 0) or (line[0] == '#'): 
            continue
        if (line[0] == '>'):
            if (cur == ''):
                cur = get_gene_id(line.strip())
                trn = line.strip()
                continue
            elif (cur not in seqs) or (seqs[cur].length < len(seq)):
                seqs[cur] = SEQ(seq, trn)
            cur = get_gene_id(line.strip())
            trn = line.strip()
            seq = ''
        else:
            seq += line.strip()
    
    if (cur not in seqs) or (seqs[cur].length < len(seq)):
        seqs[cur] = SEQ(seq, trn)

    fh.close()

    write_fasta(seqs, fasta)

    return 0

def main() -> int:
    """Filter a peptides file obtained from ensmbl for the longest isoform"""

    # get the file path
    fasta = get_arguments()

    # filter
    filter_fasta(fasta)

    return 0

if __name__ == "__main__":
    main()
