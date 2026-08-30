#!/bin/env python3

import argparse
import os

gpmap_f  = ''  # file to build prot -> gene mapping
odir     = '.'
assign_f = '' # file with RefOG <tab> Prots not found in annotations

def get_arguments() -> int:
    """get the arguments"""

    global gpmap_f, odir, assign_f

    d = "create RefOG files for genes that were not found in annotations (used for low confident assignments)"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-g", "--gpMap", help="A gpMap file that makes genes to proteins for relabeling orthogroup elements [optional]", default='', type=str)
    parser.add_argument("-a", "--assignments", help="A tsv file that contains the original RefOG file name and the protein ID", required=True, type=str)
    parser.add_argument("-o", "--outdir", help="directory to place the low confident orthogroups", type=str, default='./')

    args     = parser.parse_args()
    gpmap_f  = args.gpMap
    assign_f = args.assignments
    odir    = args.outdir

    assert os.path.isfile(assign_f), f"Could not locate file: {assign_f}"
    assert os.path.isdir(odir), f"Could not locate directory: {odir}"
    if (gpmap_f != ''):
        assert os.path.isfile(gpmap_f), f"Could not locate file: {gpmap_f}"

    if (odir[-1] == '/'):
        odir = odir[:-1]

    return 0


def load_mapping(emap: str) -> dict[str, str]:
    """read the gene to transcript map into memory"""

    egMap = dict()
    fh    = open(emap, 'r')

    for line in fh:
        fields = line.split('\t')
        gene   = fields[0]
        prot   = fields[1].strip()
        egMap[prot] = gene

    fh.close()

    return egMap

def build_low_assignments(pgMap: dict[str, str]) -> int:
    """this function will handle all the reformating"""

    global odir, assign_f

    file_handles = list()
    fh           = open(assign_f, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.split('\t')
        refgrp = fields[0]
        protID = fields[1].strip()
        geneID = pgMap.get(protID, protID) # use self if not available
        add    = True
        for i in range(len(file_handles)):
            if (file_handles[i][0] == refgrp):
                file_handles[i][1].write(geneID + '\n')
                add = False
                break
        if (add):
            outfile = f"{odir}/Recoded.{refgrp}"
            file_handles.append([refgrp, open(outfile, 'w')])
            file_handles[-1][1].write(geneID + '\n')

    fh.close()

    # shut down all the IO streams
    for i in range(len(file_handles)):
        file_handles[i][1].close()

    return 0

def main() -> int:
    """get the paths and perform functionality"""

    # get arguments
    get_arguments()

    # load an optional prot -> gene mapping for re-IDing
    pgMap = load_mapping()

    # build the orthogroups
    build_low_assignments(pgMap)

    return 0

if __name__ == "__main__":
    main()
