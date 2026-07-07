#!/bin/env python3

import argparse
import os
import gzip
from   collections import defaultdict

def get_arguments() -> tuple[str, str]:
    """get the arguments"""

    d = "Report the genes missing in the orthologs.tsv file found in a pairwise synteny cluster"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-o", "--orthos",         help="orthologs.tsv file", required=True)
    parser.add_argument("-m", "--membership", help="<prefix>__consyn_cluster_membership.tsv file", required=True)
    args  = parser.parse_args()
    ofile = args.orthos
    mfile = args.membership

    assert os.path.isfile(ofile), f"Could not locate file {ofile}"
    assert os.path.isfile(mfile), f"Could not locate file {mfile}"

    return (ofile, mfile)


def get_org_names(mfile: str) -> tuple[str, str]:
    """return the names of org_a and org_b"""

    bname = os.path.basename(mfile)
    bname = bname.replace("_consyn_cluster_membership.tsv", '')

    orgs  = bname.split('-')

    return (orgs[0], orgs[1])

def make_outname(mfile: str) -> str:
    """helper function to create an output file name"""

    org_a, org_b = get_org_names(mfile)

    return f"{org_a}-{org_b}_missingMembers.txt"

def parse_membership(mfile: str, genDict: dict[str: set[str]]) -> None:
    """return a dict of the org and a set of the gene_ids"""

    out  = make_outname(mfile)
    ofh  = open(out, 'w')
    fh   = open(mfile, 'r')
    cntA = 0
    cntB = 0

    for line in fh:
        if (len(line) == 0 or (line[0] == '#')): 
            continue
        fields = line.split('\t')
        org_a  = fields[2]
        org_b  = fields[3]
        gen_a  = fields[4]
        gen_b  = fields[10]
        if (gen_a not in genDict[org_a]):
            ofh.write(f"{org_a} {gen_a}\n")
            cntA += 1
        if (gen_b not in genDict[org_b]):
            ofh.write(f"{org_b} {gen_b}\n")
            cntB += 1

    fh.close()
    ofh.close()

    print("Found", cntA, "missing members for", org_a)
    print("Found", cntB, "missing members for", org_b)


def get_orthos(ofile: str, org_a: str, org_b: str) -> dict[str: set[str]]:
    """return a Dict[str: set()] for org a and org b"""

    genDict = defaultdict(set)
    fh      = gzip.open(ofile, "rt") if ofile.endswith(".gz") else open(ofile, 'r')

    for line in fh:
        if (line[0] == '#'):
            continue
        fields = line.split('\t')
        org    = fields[4]
        if (org == org_a or org == org_b):
            gen = fields[5]
            genDict[org].add(gen)

    fh.close()

    return genDict

def main() -> int:
    """entry point to this little helper script"""

    # get arguments
    ofile, mfile = get_arguments()

    # get the names of org a and b
    org_a, org_b = get_org_names(mfile)

    # get the genes from orthologs.tsv
    genDict = get_orthos(ofile, org_a, org_b)

    # do the comparison
    parse_membership(mfile, genDict)

    return 0

if __name__ == "__main__":
    main()
