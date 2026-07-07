#!/bin/env python3

import argparse
import os
import gzip
from collections import defaultdict

def get_arguments() -> tuple:
    """get the arguments"""

    d = "find the missing orthologs from two orthologs.tsv files (i.e., find what is missing in a that are present in b"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-a", "--a_orthos", help="First orthologs.tsv file", required=True)
    parser.add_argument("-b", "--b_orthos", help="Second orthologs.tsv file", required=True)
    args    = parser.parse_args()
    a = args.a_orthos
    b = args.b_orthos

    assert os.path.isfile(a), f"Could not locate file {a}"
    assert os.path.isfile(b), f"Could not locate file {b}"

    return (a, b)


def get_b_orthos(b: str) -> dict:
    """return a Dict[str: set()] for each spp"""
    Bgenes = defaultdict(set)
    fh     = gzip.open(b, "rt") if b.endswith(".gz") else open(b, 'r')

    for line in fh:
        if (line[0] == '#'):
            continue
        fields = line.split('\t')
        spp    = fields[4]
        gen    = fields[5]
        Bgenes[spp].add(gen)

    fh.close()

    return Bgenes

def compare_orthos(a: str, Bgenes: dict) -> None:
    """write out genes in orthos1 that are not found in orthos2"""
    
    fh   = gzip.open(a, "rt") if a.endswith(".gz") else open(a, 'r')
    ofh  = open("MissingOrthologs.txt", 'w')
    seen = defaultdict(set)
    ofh2 = open("multiples.txt", 'w')
    mcnt = 0

    for line in fh:
        if (line[0] == '#'):
            continue
        
        fields = line.split('\t')
        spp    = fields[4]
        gen    = fields[5]

        if (gen not in Bgenes[spp]):
            ofh.write(f"{spp} {gen}\n")

        if (gen in seen[spp]):
            mcnt += 1
            ofh2.write(f"{spp} {gen}\n")
        seen[spp].add(gen)

    fh.close()
    ofh.close()
    ofh2.close()

    print("Number of multiples:", mcnt)

def main() -> int:
    """entry point to this little subprogram"""

    # get arguments
    a, b = get_arguments()

    # get the genes from the previous run
    Bgenes = get_b_orthos(b)

    # do the comparison
    compare_orthos(a, Bgenes)

    return 0

if __name__ == "__main__":
    main()
