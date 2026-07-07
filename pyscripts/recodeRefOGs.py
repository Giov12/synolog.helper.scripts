#!/bin/env python3

import argparse
import os
import glob
import sys


genes = set()


def get_arguments() -> tuple[str, str, str, str]:
    """get the arguments"""

    d = "create the reference orthogroups swapping transcript IDs for gene IDs"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-g", "--geneMaps", help="directory containing gpMap.tsv files", required=True)
    parser.add_argument("-r", "--RefOGs", help="directory containing the RefOG txt files", required=True)
    parser.add_argument("-m", "--map", help="gene/transcript map for RefOF txt files", required=True)
    parser.add_argument("-o", "--outdir", help="directory to place the recoded orthogroups", type=str, default='./')

    args = parser.parse_args()
    gdir = args.geneMaps
    rdir = args.RefOGs
    tmap = args.map
    odir = args.outdir

    assert os.path.isdir(gdir), f"Could not locate directory: {gdir}"
    assert os.path.isdir(rdir), f"Could not locate directory: {rdir}"
    assert os.path.isdir(odir), f"Could not locate directory: {odir}"
    assert os.path.isfile(tmap), f"Could not locate file: {tmap}"

    if (odir[-1] == '/'):
        odir = odir[:-1]

    return (gdir, rdir, odir, tmap)

def load_genes(gdir: str) -> dict[str, str]:
    """construct a transcript to gene map"""

    global genes

    pgMap = dict()
    
    if (gdir[-1] == '/'):
        maps = glob.glob(f"{gdir}*.gpMap.tsv")
    else:
        maps = glob.glob(f"{gdir}/*.gpMap.tsv")

    if (len(maps) == 0):
        msg = f"Could not find any *gpMap.tsv files in {gdir}"
        sys.exit(msg)
    
    for map_ in maps:
        fh = open(map_, 'r')
        for line in fh:
            line    = line.strip()
            fields  = line.split('\t')
            isoform = fields[1].split('.')[0]
            gene    = fields[0]
            pgMap[isoform] = gene
            genes.add(gene)
        fh.close()

    return pgMap

def read_map(tmap: str) -> dict[str, str]:
    """read the gene to transcript map into memory"""

    tgMap = dict()
    fh    = open(tmap, 'r')

    for line in fh:
        fields = line.split('\t')
        gene   = fields[0]
        prot   = fields[1].strip()
        tgMap[prot] = gene

    fh.close()

    return tgMap

def prot_to_genes(pgMap: dict[str, str], tgMap: dict[str, str], rdir: str, odir: str) -> int:
    """this function will handle all the reformating"""

    global genes

    if (rdir[-1] == '/'):
        ogrps = glob.glob(f"{rdir}RefOG*.txt")
    else:
        ogrps = glob.glob(f"{rdir}/RefOG*.txt")

    if (len(ogrps) == 0):
        msg = f"Could not find any RefOG*.txt files in {rdir}"
        sys.exit(msg)

    for ogrp in ogrps:
        fname = os.path.basename(ogrp)
        out   = f"{odir}/Recoded.{fname}"
        inFH  = open(ogrp, 'r')
        outFH = open(out, 'w')

        for line in inFH:
            isoform = line.strip()
            gene    = ''
            if (isoform in pgMap):
                gene = pgMap[isoform]
            elif (isoform in tgMap):
                gene = tgMap[isoform]
                if (gene not in genes):
                    gene = '#' + tgMap[isoform]
            if (gene == ''):
                print("Couldn't find anything for", isoform)
                continue
            outFH.write(gene + '\n')
        inFH.close()
        outFH.close()

    return 0

def main() -> int:
    """get the paths and perform functionality"""

    # get arguments
    gdir, rdir, odir, tmap = get_arguments()

    # construct the protein mapping
    pgMap = load_genes(gdir)

    # get the original mapping
    tgMap = read_map(tmap)

    # process all the groups
    prot_to_genes(pgMap, tgMap, rdir, odir)

    return 0

if __name__ == "__main__":
    main()
