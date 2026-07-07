#!/bin/env python3

import argparse
import os
import gzip

"""
filter a protein fasta file for only transcripts in the gene-protein map
"""

def get_arguments() -> tuple[str, str]:
    """get the arguments"""

    d = "Some python code to filter protein sequences for transcripts present in the gene protein map"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-f", "--fasta", help="fasta file", required=True)
    parser.add_argument("-g", "--gpmap", help="gene-protein map", required=True)

    args  = parser.parse_args()
    fasta = args.fasta
    gpMap = args.gpmap

    assert os.path.isfile(fasta), f"Could not locate file: {fasta}"
    assert os.path.isfile(gpMap), f"Could not locate file: {gpMap}"
    
    return (fasta, gpMap)

def read_gpMap(gpMap: str) -> set[str]:
    """read in the transcript ids into a set"""

    transcripts = set()
    fh          = open(gpMap, 'r')

    for line in fh:
        fields = line.strip().split('\t')
        transcripts.add(fields[1])
    fh.close() 

    return transcripts

def get_gene_name(line: str) -> str:
    """get the gene name from the header"""

    fields = line.split(' ')
    gene   = fields[3].split(':')[1].split('.')[0]

    return gene

def filter_fasta(fasta: str, transcripts: set[str]) -> int:
    """filter a fasta file for transcripts in the set"""

    name = os.path.basename(fasta)
    name = "filtered." + name

    fh   = gzip.open(fasta, "rt") if fasta.endswith(".gz") else open(fasta, 'r')
    ofh  = open(name, 'w')
    keep = False

    for line in fh:
        if (len(line) == 0) or (line[0] == '#'): 
            continue
        if (line[0] == '>'):
            # gene = get_gene_name(line)
            transcript = line[1:].split(' ')[0].split('.')[0]
            keep = transcript in transcripts
            if (keep):
                ofh.write(line)
        elif (keep):
            ofh.write(line)
    
    fh.close()
    ofh.close()

    return 0

def main() -> int:
    """Filter fasta file for transcripts in gene protein map"""

    # get arguments
    fasta, gpMap = get_arguments()

    # get the transcripts
    transcripts = read_gpMap(gpMap)

    # filter
    filter_fasta(fasta, transcripts)

    return 0

if __name__ == "__main__":
    main()
