#!/bin/env python3

import argparse
import os
import glob
import textwrap
import gzip

"""
filter a protein fasta file for the longest transcript
"""

class rec:
    def __init__(self, seq: str, header: str):
        self.seq    = seq
        self.header = header
        self.len    = len(seq)

        
def get_arguments() -> str:
    """get the arguments"""

    d = "Some python code to filter protein sequences to the longest transcript"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-d", "--dir", help="directory containing .aa files", required=True)

    args = parser.parse_args()
    fdir = args.dir

    assert os.path.isdir(fdir), f"Could not locate directory: {fdir}"
    
    return fdir

def get_aa_files(fdir: str) -> list:
    """get the aa files from the provided directory"""

    aaFiles = list()

    if (fdir[-1] == '/'):
        aaFiles = glob.glob(f"{fdir}*.aa")
        if (len(aaFiles) == 0):
            aaFiles = glob.glob(f"{fdir}*.aa.gz")
    else:
        aaFiles = glob.glob(f"{fdir}/*.aa")
        if (len(aaFiles) == 0):
            aaFiles = glob.glob(f"{fdir}/*.aa.gz")

    assert (len(aaFiles) > 0), f"Could not locate any .aa files in {fdir}"

    return aaFiles

def get_geneName(seqId: str) -> str:
    """remove the transcript extension from a sequence id"""

    c = seqId.split(' ')

    if (len(c) == 1):
        return seqId
    
    geneId = c[3].split(':')[1]

    return geneId

def write_fasta(seqs: dict, fasta: str) -> None:
    """write the Seqs dictionary to a fasta file"""

    fname = os.path.basename(fasta)
    fname = f"Filtered.{fname}"

    fh = open(fname, 'w')

    for seqId, seq_rec in seqs.items():
        fh.write(seq_rec.header + '\n')
        fh.write(textwrap.fill(seq_rec.seq, width=60))
        fh.write('\n')

    fh.close()


def filter_fasta(fasta: str) -> None:
    """filter a fasta file for a single representative transcript per gene"""

    seqs = dict()
    fh   = gzip.open(fasta, "rt") if fasta.endswith(".gz") else open(fasta, 'r')
    head = ''
    cur  = ''
    seq  = ''

    for line in fh:
        if (len(line) == 0) or (line[0] == '#'): 
            continue
        line = line.strip()
        if (line[0] == '>'):
            if (cur == ''):
                cur  = get_geneName(line)
                head = line
                continue
            elif (cur not in seqs) or (seqs[cur].len < len(seq)):
                seqs[cur] = rec(seq, head)
            cur  = get_geneName(line)
            head = line
            seq  = ''
        else:
            seq += line
    
    if (cur not in seqs) or (seqs[cur].len < len(seq)):
        seqs[cur] = rec(seq, head)

    fh.close()

    write_fasta(seqs, fasta)

def process_files(aaFiles: list) -> None:
    """func() that will iterate over the amino acid sequences"""

    for fasta in aaFiles:
        filter_fasta(fasta)


def main() -> int:
    """Filter a AA file for a single transcript (longest) per gene"""

    # get arguments
    fdir = get_arguments()

    # get the fasta files
    aaFiles = get_aa_files(fdir)

    # begin the process
    process_files(aaFiles)

    return 0

if __name__ == "__main__":
    main()
