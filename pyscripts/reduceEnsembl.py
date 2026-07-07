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

def get_arguments() -> tuple[str, str]:
    """get the arguments"""

    d = "Some python code to filter protein sequences to the longest isoform from an ensembl peptide " \
        "fasta and found in the assocaited annotation"
    
    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-f", "--fasta", help="A *.pep.all.fa file", required=True)
    parser.add_argument("-g", "--gtf", help="A *.gtf.gz file", required=True)

    args  = parser.parse_args()
    fasta = args.fasta
    gtf   = args.gtf

    assert os.path.isfile(fasta), f"Could not locate file: {fasta}"
    assert os.path.isfile(gtf), f"Could not locate file: {gtf}"
    
    return (fasta, gtf)

def get_attributes(field: str) -> dict[str, str]:
    """return a dict of the transcript's attributes"""

    attr      = dict()
    subfields = field.split(';')

    for subfield in subfields:
        subfield  = subfield.strip(' "')
        space     = subfield.find(' ') # find first space
        key       = subfield[:space]
        attrib    = subfield[space + 1:].strip('"')
        attr[key] = attrib

    return attr

def load_isoforms(gtf: str) -> set[str]:
    """return a set of isoform IDs to retain"""

    isforms = set()
    fh      = gzip.open(gtf, "rt") if gtf.endswith(".gz") else open(gtf, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.split('\t')
        if (fields[2] != "transcript"):
            continue
        attrib = fields[-1].strip()
        attrib = get_attributes(attrib)
        if ("transcript_id" not in attrib):
            continue
        tID = attrib["transcript_id"]
        if ("transcript_version" in attrib):
            tID = tID + '.' + attrib["transcript_version"]
        isforms.add(tID)

    fh.close()

    return isforms

def get_seq_ids(header: str) -> tuple[str, str]:
    """return the transcript and gene id from an ensembl header"""

    sects = header.split(' ')
    gene  = sects[3]
    gene  = gene.split(':')[1]
    gene  = gene.split('.')[0]
    tID   = sects[4]
    tID   = tID.split(':')[1]

    return (gene, tID)

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

def filter_fasta(fasta: str, isoforms: set[str]) -> int:
    """filter a fasta file for a single representative isoform per gene"""

    seqs = dict()
    fh   = gzip.open(fasta, "rt") if fasta.endswith(".gz") else open(fasta, 'r')
    gID  = ''
    tID  = ''
    seq  = ''

    for line in fh:
        if (len(line) == 0) or (line[0] == '#'): 
            continue
        line = line.strip()
        if (line[0] == '>'):
            if (gID == ''):
                gID, tID = get_seq_ids(line)
                header   = line
                continue
            elif (gID not in seqs) or (seqs[gID].length < len(seq)):
                if (tID in isoforms):
                    seqs[gID] = SEQ(seq, header)
            gID, tID = get_seq_ids(line.strip())
            header   = line
            seq = ''
        else:
            seq += line
    
    if (gID not in seqs) or (seqs[gID].length < len(seq)):
        if (tID in isoforms):
            seqs[gID] = SEQ(seq, header)

    fh.close()

    write_fasta(seqs, fasta)

    return 0

def main() -> int:
    """Filter a peptides file obtained from ensmbl for the longest isoform & annotation"""

    # get the file paths
    fasta, gtf = get_arguments()

    # load all the protein IDs
    isoforms = load_isoforms(gtf)

    # filter
    filter_fasta(fasta, isoforms)

    return 0

if __name__ == "__main__":
    main()
