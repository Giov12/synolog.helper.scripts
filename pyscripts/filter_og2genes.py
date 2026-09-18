#!/bin/env python3

import argparse
import os
import gzip
import sys

og_ids    = set()
spp_ids   = set()
og2genes  = '' # input files
og_ids_f  = ''
spp_ids_f = '' 

def get_arguments() -> int:
    """get the arguments"""

    global og2genes, og_ids_f, spp_ids_f

    d = "Collect the gene records from odb12v2_OG2genes.tab.gz using OG and species IDs"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-o", "--ogs", help="odb12v2_OG2genes.tab.gz file", required=True, type=str)
    parser.add_argument("-g", "--groups", help="single columns list of OG IDs", required=True, type=str)
    parser.add_argument("-s", "--spp", help="single columns list of species IDs", required=True, type=str)
    args      = parser.parse_args()
    og2genes  = args.ogs
    og_ids_f  = args.groups
    spp_ids_f = args.spp

    assert os.path.isfile(og2genes), f"Could not locate file {og2genes}"
    assert os.path.isfile(og_ids_f), f"Could not locate file {og_ids_f}"
    assert os.path.isfile(spp_ids_f), f"Could not locate file {spp_ids_f}"

    return 0

def get_og_ids() -> int:
    """load the OG ids into memory"""

    global og_ids, og_ids_f

    fh = gzip.open(og_ids_f, "rt") if og_ids_f.endswith(".gz") else open(og_ids_f, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        og_id = line.strip()
        og_ids.add(og_id)

    fh.close()

    print(f"Loaded {len(og_ids)} target groups")

    return 0

def get_spp_ids() -> int:
    """load the species ids into memory"""

    global spp_ids, spp_ids_f

    fh = gzip.open(spp_ids_f, "rt") if spp_ids_f.endswith(".gz") else open(spp_ids_f, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        spp_id = line.strip()
        spp_ids.add(spp_id)

    fh.close()

    print(f"Loaded {len(spp_ids)} target species")

    return 0

def make_output_name() -> str:
    """helper function to create the output name"""

    global og_ids, spp_ids

    out = f"filtered.{len(og_ids)}.OGs.{len(spp_ids)}.spp.odb12v2_OG2genes.tab.gz"

    return out

def get_og_genes() -> int:
    """the work horse of this application"""

    global og_ids, spp_ids, og2genes

    fh  = gzip.open(og2genes, "rt") if og2genes.endswith(".gz") else open(og2genes, 'r')
    ofh = gzip.open(make_output_name(), "wt")
    tot = 0

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue

        cnt = 0
        for i in range(len(line)):
            if (line[i] == '\t' or line[i] == ':'):
                cnt += 1

        # safety check
        if (cnt != 2):
            msg = f"Invalid line found:\n{line}"
            sys.exit(msg)

        fields = line.split('\t')
        og_id  = fields[0]

        # ogroup checkpoint
        if (og_id not in og_ids):
            continue

        # species checkpoint
        idx = fields[1].find(':')
        spp_id = fields[1][:idx]
        if (spp_id not in spp_ids):
            continue

        ofh.write(line)
        tot += 1

    fh.close()
    ofh.close()

    print(f"Retained {tot} records")

    return 0

def main() -> int:
    """entry point to this little subprogram"""

    # get arguments
    get_arguments()

    # load the target OGs into memory
    get_og_ids()

    # load the species IDs into memory
    get_spp_ids()

    # start processing
    get_og_genes()

    return 0

if __name__ == "__main__":
    main()
