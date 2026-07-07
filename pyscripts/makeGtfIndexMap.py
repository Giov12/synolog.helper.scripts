#!/bin/env python3

import argparse
import os
import gzip


def get_arguments() -> str:
    """get the arguments"""

    d = "create a mapping of the gene indices based on the chromosome from a gtf"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-g", "--gtf", help="gtf file", required=True)
    args = parser.parse_args()
    gtf  = args.gtf

    assert os.path.isfile(gtf), f"Could not locate file: {gtf}"

    return gtf

def get_gi(field: str) -> str:
    """get the gene_id if available"""

    if (field.count(';') == 0):
        return field.strip()
    
    subfields = field.split(';')
    gi        = "NA"
    
    for subfield in subfields:
        subfield = subfield.strip().split(' ')
        if (subfield[0] == "gene_id"):
            gi = subfield[1].replace('"', '')
            break

    return gi

def make_outname(gtf: str) -> str:
    """create the output name for the gtf file"""

    bname = os.path.basename(gtf)
    r     = ".gtf.gz" if gtf.endswith(".gz") else ".gtf"
    return bname.replace(r, ".GeneIndicesMap.tsv")

def make_map(gtf: str) -> None:
    """the work horse of this program"""

    fh  = gzip.open(gtf, "rt") if gtf.endswith(".gz") else open(gtf, 'r')
    cur = None
    idx = 0
    out = make_outname(gtf)
    ofh = open(out, 'w')

    for line in fh:
        if (len(line) == 0) or (line[0] == '#'):
            continue
        fields = line.split('\t')
        if (fields[2] != "gene"):
            continue
        Chr       = fields[0]
        subfields = fields[-1]
        gene_id   = get_gi(subfields)
        if (cur == None):
            idx += 1
            cur = Chr
            ofh.write(f"{Chr}\t{idx}\t{gene_id}\n")
        elif (cur == Chr):
            idx += 1
            ofh.write(f"{Chr}\t{idx}\t{gene_id}\n")
        else:
            idx = 1
            cur = Chr
            ofh.write(f"{Chr}\t{idx}\t{gene_id}\n")

    fh.close()
    ofh.close()

def main() -> int:
    """get the gtf file and create a tsv file to serve as a map"""

    # get arguments
    gtf = get_arguments()

    # get going
    make_map(gtf)

    return 0

if __name__ == "__main__":
    main()
