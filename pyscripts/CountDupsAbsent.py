#!/bin/env python3

import argparse
import os
import gzip
from   collections import defaultdict

def get_arguments() -> tuple:
    """get the arguments"""

    d = "Count the number of genes in duplicates not found in orthologs,tsv"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-o", "--orthos", help="orthologs.tsv file", required=True)
    parser.add_argument("-d", "--dups", help="duplicates.tsv file", required=True)
    args  = parser.parse_args()
    ortho = args.orthos
    dups  = args.dups

    assert os.path.isfile(ortho), f"Could not locate file {ortho}"
    assert os.path.isfile(dups), f"Could not locate file {dups}"

    return (ortho, dups)


def get_orthos(ortho: str) -> dict:
    """return a Dict[str: set()] for each spp"""
    genDict = defaultdict(set)
    fh      = gzip.open(ortho, "rt") if ortho.endswith(".gz") else open(ortho, 'r')

    for line in fh:
        if (line[0] == '#'):
            continue
        fields = line.split('\t')
        spp    = fields[4]
        gen    = fields[5]
        genDict[spp].add(gen)

    fh.close()

    return genDict

def compare_orthos(dups: str, genDict: dict) -> None:
    """write out genes in orthos1 that are not found in orthos2"""
    
    fh  = gzip.open(dups, "rt") if dups.endswith(".gz") else open(dups, 'r')
    ofh = open("MissingDuplicates.txt", 'w')
    cnt = 0

    for line in fh:
        if (line[0] == '#'):
            continue
        
        fields = line.split('\t')
        spp    = fields[1].strip()
        gen    = fields[5].strip()

        if (gen not in genDict[spp]):
            cnt += 1
            ofh.write(f"{spp} {gen}\n")

    fh.close()
    ofh.close()

    print("Missing:", cnt, "genes")

def main() -> int:
    """entry point to this little subprogram"""

    # get arguments
    ortho, dups = get_arguments()

    # get the genes from orthologs.tsv
    genDict = get_orthos(ortho)

    # do the comparison
    compare_orthos(dups, genDict)

    return 0

if __name__ == "__main__":
    main()
