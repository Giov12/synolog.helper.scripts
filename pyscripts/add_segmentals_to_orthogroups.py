#!/bin/env python3

import argparse
import os
import glob
import gzip
import matplotlib.pyplot as plt
from  collections import defaultdict

orthologs   = {}               # str -> str   [str will be an string int | group ID]
orthogroups = defaultdict(set) # str -> set[str]
candidates  = defaultdict(list) # str -> list[tuple[str, str]]
path        = ''
retroset    = set()

def get_arguments() -> int:
    """get the arguments"""

    global path

    desc  = "Add unplaced segmental duplications to an orthogroup"
    phelp = "Path to the directory containing the orthogroups.tsv and *segmentalduplications.tsv files from synolog"

    parser = argparse.ArgumentParser(description = desc)
    parser.add_argument("-p", "--path", help=phelp, required=True)
    args = parser.parse_args()
    path = args.path

    assert os.path.isdir(path), f"Could not locate directory {path}"

    if (path[-1] == '/'):
        path = path[:-1] # remove path char

    return 0

def get_segmental_duplications() -> list[str]:
    """get the segmental duplication files in the provided path"""

    global path

    suffix   = "_segmentalduplications.tsv"
    segfiles = glob.glob(f"{path}/*{suffix}")

    assert len(segfiles) > 0, f"Could not locate any *{suffix} files in {path}/"

    return segfiles

def load_orthologs() -> int:
    """load the current members into memory"""

    global path, orthologs, orthogroups, retroset

    infile = f"{path}/orthologs.tsv"

    assert os.path.isfile(infile), f"Could not locate orthologs.tsv in {path}/"

    fh = open(infile, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#' or line[0] == '/n'):
            continue
        fields = line.split('\t')
        assert len(fields) == 12, f"Error: Expected 12 columns in orthologs.tsv. Offending line:\n{line}"
        ogrp     = fields[0]
        spp      = fields[4]
        gene_id  = fields[5]
        elemType = fields[6]
        member   = f"{spp}_{gene_id}"
        if (elemType != 'R'):
            orthogroups[ogrp].add(member)
            orthologs[member] = ogrp
        else:
            retroset.add(member)

    fh.close()

    print(f"Loaded {len(orthologs)} orthologs across {len(orthogroups)}")
    print(f"Loaded {len(retroset)} genes flagged as retrogenes")
    
    return 0

def load_candidates() -> int:
    """find segmental duplications not in orthogroups.tsv"""

    global orthologs, retroset, candidates

    segfiles = get_segmental_duplications()

    for segfile in segfiles:
        fh = open(segfile, 'r')

        for line in fh:
            if (len(line) == 0 or line[0] == '#'):
                continue
            fields = line.split('\t')
            assert len(fields) == 22, \
                f"Error: Expected 22 columns in {os.path.basename(segfile)}. Offending line:\n{line}"
            spp    = fields[2]
            dup_id = fields[3]
            cand   = f"{spp}_{dup_id}"
            if (cand in orthologs or cand in retroset):
                continue # was placed into orthologs.tsv
            gene_id = f"{spp}_{fields[9]}"         # duplicate of this gene
            other   = f"{fields[15]}_{fields[16]}" # species B
            candidates[cand].append((gene_id, other))

        fh.close()

    print(f"Loaded {len(candidates)} candidate segmental duplications")
    return 0

def main() -> int:
    """entry point to this helper script"""

    # get the input path
    get_arguments()

    # load the core members
    load_orthologs()

    # filter the segmental duplications for candidate genes
    load_candidates()

    # now to 


    return 0

if __name__ == "__main__":
    main()
