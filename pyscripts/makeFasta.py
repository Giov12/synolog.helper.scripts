#!/bin/env python3

import argparse
import os
import glob
import gzip
import textwrap

def get_arguments() -> tuple:
    """get the arguments"""

    d = "a companion script for compareOrthoMethods.py"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-r", "--revisit", help="RevisitList.txt file", required=True)
    parser.add_argument("-o", "--outdir", help="outdirectory", default="./")
    parser.add_argument("-m", "--maps", help="directory containing any map gene->protein id mappping", default='')
    parser.add_argument("-f", "--fasta", help="directory containing the fasta files for the sequences of interest", default='')
    args = parser.parse_args()
    revf = args.revisit
    fdir = args.fasta
    mdir = args.maps
    odir = args.outdir

    assert os.path.isfile(revf), f"Could not locate file {revf}"
    assert os.path.isdir(fdir), f"Could not locate directory {fdir}"
    assert os.path.isdir(mdir), f"Could not locate directory {mdir}"
    assert os.path.isdir(odir), f"Could not locate directory {odir}"

    return (revf, fdir, mdir, odir)

def get_maps(mdir: str) -> dict:
    """get any gene->protein mapping for easier comparison"""

    if (mdir == ''): return {}

    mfiles = []

    if (mdir.endswith('/')):
        mfiles = glob.glob(f"{mdir}*gpMap.tsv")
    else:
        mfiles = glob.glob(f"{mdir}/*gpMap.tsv")

    assert len(mfiles) > 0, f"Could not locate any *gpMap.tsv files at {mdir}"

    gpMap = {}

    for f in mfiles:
        fh         = open(f, 'r')
        spp        = os.path.basename(f).split('.')[0]
        gpMap[spp] = {}

        for line in fh:
            fields           = line.strip().split('\t')
            gene             = fields[0]
            prot             = fields[1]
            gpMap[spp][prot] = gene

        fh.close()

    return gpMap

def get_fasta(fdir: str, gpMap: dict) -> dict:
    """get the fasta files""" 

    fastas  = []
    seqDict = {}

    if (fdir[-1] != '/'):
        fastas = glob.glob(f"{fdir}/F*aa.gz")
    else:
        fastas = glob.glob(f"{fdir}F*aa.gz")

    assert (len(fastas) > 0), f"Could not locate any *.aa.gz files at {fdir}"

    for fasta in fastas:
        spp = os.path.basename(fasta).split('.')[2]
        seqDict[spp] = read_fasta(fasta, gpMap, spp)

    return seqDict

def strip_transcript(seqId: str) -> str:
    """remove the transcript extension from a sequence id"""

    c = seqId.split(' ')[0]
    c = c.split('.')

    if (len(c) == 1):
        return seqId
    
    if (c[0][0] == '>'):
        c[0] = c[0][1:]

    # if (c[-1][0].lower() in ['t', '1']):
    #     return '.'.join(c[:-1])
    
    return c[0]

def read_fasta(fasta: str, gpMap: dict, spp: str) -> dict:
    """return a dict with sequence header as key and seq as value"""

    seqDict = {}

    fh  = gzip.open(fasta, "rt") if fasta.endswith(".gz") else open(fasta, 'r')
    seq = ''
    hdr = ''

    for line in fh:
        line = line.strip()
        if (line[0] == '>'):
            if (seq != ''):
                seqDict[gene] = seq
                seq = ''
            hdr = strip_transcript(line)
            if (hdr == "ENSLOCP000000000011"):
                print(line)
            if (spp in gpMap):
                gene = gpMap[spp][hdr]
            else:
                gene = hdr
        else:
            seq += line

    fh.close()

    seqDict[gene] = seq

    return seqDict

def get_ids(seqDict: dict) -> dict:
    """return a dict of sets with the sequence headers"""
    transKey = {spp : set() for spp in seqDict}

    for spp, seqs in seqDict.items():
        for hdr in seqs:
            transKey[spp].add(hdr)

    return transKey

def read_revisit(revf: str, seqDict: dict, odir: str) -> None:
    """return a dictionary of the gene ids that need to be written"""

    i    = 0
    fh   = gzip.open(revf, "rt") if revf.endswith(".gz") else open(revf, 'r')
    odir = odir[:-1] if odir.endswith('/') else odir

    for line in fh:
        ofh = open(f"{odir}/Ortho{i}.aa", 'w')
        genes = line.strip().split(", ")
        for entry in genes:
            spp, gene = entry.split(':')
            seq       = seqDict[spp][gene]
            ofh.write(f">{entry}\n")
            ofh.write(textwrap.fill(seq, width=60))
            ofh.write('\n')

        ofh.close()
        i += 1

    fh.close()


def main() -> int:
    """entry point to this little subprogram"""

    # get arguments
    revf, fdir, mdir, odir = get_arguments()

    # read in the maps
    gpMap = get_maps(mdir)

    # get the seqs
    seqDict = get_fasta(fdir, gpMap)

    # write out the individual grps
    read_revisit(revf, seqDict, odir)

    return 0

if __name__ == "__main__":
    main()
