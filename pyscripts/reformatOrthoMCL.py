#!/bin/env python3

import argparse
import os
import glob
import gzip
from collections import defaultdict


def get_arguments() -> tuple:
    """get the arguments"""

    d = "reformat the OrthoMCL final output to a tsv for easier comparisons"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-o", "--orthomcl", help="orthomcl.final.out file", required=True)
    parser.add_argument("-s", "--spp", help="text file containing names of species for columns (1 spp per line)", required=True)
    parser.add_argument("-m", "--maps", help="directory containing any map gene->protein id mappping", default='')
    args = parser.parse_args()
    omcl = args.orthomcl
    spp  = args.spp
    mdir = args.maps

    assert os.path.isfile(omcl), f"Could not locate directory {omcl}"
    assert os.path.isfile(spp), f"Could not locate directory {spp}"

    if (mdir != ''):
        assert os.path.isdir(mdir), f"Could not locate directory {mdir}"

    return (omcl, spp, mdir)

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
            fields = line.strip().split('\t')
            gene   = fields[0]
            prot   = fields[1]
            gpMap[spp][prot] = gene

        fh.close()

    return gpMap

def get_spp(spp: str) -> list:
    """return a list of all the focal spp"""

    sppList = []

    fh = gzip.open(spp, "rt") if spp.endswith(".gz") else open(spp, 'r')

    for line in fh:
        sppList.append(line.strip())

    fh.close()

    return sppList


def reformat_omcl(omcl: str, gpMap: dict, sppList: list) -> None:
    """reformat each line & write out afterwards"""

    outfile = "Reformat." + os.path.basename(omcl)
    ofh     = open(outfile, 'w')
    fh      = gzip.open(omcl, "rt") if omcl.endswith(".gz") else open(omcl, 'r')

    header  = ["#GroupID"] + sppList
    header  = '\t'.join(header) + '\n'
    ofh.write(header)

    for line in fh:
        fields = line.strip().split(' ')
        gi     = fields[0][:-1]
        row    = [gi]
        spp_dic = {s : [] for s in sppList}
        for i in range(1, len(fields)):
            subfield = fields[i].split('|')
            spp      = subfield[0]
            # gene     = subfield[1].split('.')[0]
            gene     = subfield[1]
            if (spp in gpMap):
                if (gene not in gpMap[spp]):
                    gene = ''.join(gene.split('.')[:-1])
                gene = gpMap[spp][gene]
            spp_dic[spp].append(gene)
        for spp in sppList:
            row.append(', '.join(spp_dic[spp]) if spp_dic[spp] else '')
        outline = '\t'.join(row) + '\n'
        ofh.write(outline)

    fh.close()
    ofh.close()


def main() -> int:
    """get the gtf file and create a tsv file to serve as a map"""

    # get arguments
    omcl, spp, mdir = get_arguments()

    # get the files
    gpMap   = get_maps(mdir)
    sppList = get_spp(spp)

    reformat_omcl(omcl, gpMap, sppList)

    return 0

if __name__ == "__main__":
    main()
