#!/bin/env python3

import argparse
import os
import gzip
from   collections import defaultdict

def get_arguments() -> tuple:
    """get the arguments"""

    d = "a companion script for compareOrthoMethods.py"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-r", "--revisit", help="RevisitList.txt file", required=True)
    parser.add_argument("-o", "--outdir", help="outdirectory", default="./")
    parser.add_argument("-d", "--diamond", help="diamond results in outformat tsv [6]", required=True)
    parser.add_argument("-f", "--fasta", help="directory containing the fasta files for the sequences of interest", default='')
    args = parser.parse_args()
    revf = args.revisit
    odir = args.outdir
    diam = args.diamond

    assert os.path.isfile(revf), f"Could not locate file {revf}"
    assert os.path.isfile(diam), f"Could not locate file {diam}"
    assert os.path.isdir(odir), f"Could not locate directory {odir}"

    return (revf, diam, odir)

def read_diamond(diam: str) -> dict:
    """read the hits into a dict for later printing"""

    dhits = defaultdict(list)
    fh    = gzip.open(diam, "rt") if diam.endswith(".gz") else open(diam, 'r')

    for line in fh:
        fields = line.split('\t')
        dhits[fields[0]].append(line)
    fh.close()

    return dhits

def write_revisit(revf: str, dhits: dict, odir: str) -> None:
    """return a dictionary of the gene ids that need to be written"""

    i    = 0
    fh   = gzip.open(revf, "rt") if revf.endswith(".gz") else open(revf, 'r')
    odir = odir[:-1] if odir.endswith('/') else odir

    for line in fh:
        ofh = open(f"{odir}/Ortho{i}.hits.tsv", 'w')
        genes = line.strip().split(", ")
        for entry in genes:
            entry_hits = dhits[entry]
            for hit in entry_hits:
                ofh.write(hit)
        ofh.close()
        i += 1

    fh.close()


def main() -> int:
    """entry point to this little subprogram"""

    # get arguments
    revf, diam, odir = get_arguments()

    # read in the diamond hits
    dhits = read_diamond(diam)

    # write out the hits
    write_revisit(revf, dhits, odir)

    return 0

if __name__ == "__main__":
    main()
