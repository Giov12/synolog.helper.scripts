#!/bin/env python3

import argparse
import os
import gzip
import glob
from   collections import defaultdict

"""
    count the number of orthogorups from a reformated tsv file with >=2 species
"""

        
def get_arguments() -> str:
    """get the arguments"""

    d = "Some python code to compare either orthomcl or orthofinder to synolog"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-r", "--reformated",  help="reformated tsv of either orthofinder or orthomcl", required=True)

    args = parser.parse_args()
    omtd = args.reformated

    assert os.path.isfile(omtd), f"Could not locate file: {omtd}"

    return omtd

def parse_reformat_tsv(tsv: str) -> None:
    """count the number of ortholog groups with at least 2 spp"""

    fh      = gzip.open(tsv, "rt") if tsv.endswith(".gz") else open(tsv, 'r')
    tot     = 0
    ngenes  = 0

    for line in fh:
        fields = line.strip().split('\t')
        if (line[0] == '#'):
            continue
        spcnt = 0
        gcnt  = 0
        for i in range(1, len(fields)):
            if (fields[i] != ''):
                spcnt += 1
                subfield = fields[i].split(", ")
                for gene in subfield:
                    gcnt += 1
        if (spcnt > 1):
            tot += 1
            ngenes += gcnt
        
    fh.close()

    print(f"Total number of multi-species orthogroups:", tot)
    print(f"Total number of genes in multi-species orthogroups:", ngenes)

def main() -> int:
    """Count the number of multi-species orthogroups from either the orthofinder or orthomcl software"""

    # get arguments
    omtd = get_arguments()

    # count
    parse_reformat_tsv(omtd)

    return 0

if __name__ == "__main__":
    main()
