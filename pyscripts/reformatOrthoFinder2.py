#!/bin/env python3

import argparse
import os
import glob
import gzip

def get_arguments() -> tuple:
    """get the arguments"""

    d = "reformat the OrtoFinder output to a tsv for easier comparisons"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-o", "--orthofinder", help="N0.tsv file", required=True)
    parser.add_argument("-s", "--spp", help="text file containing names of species for columns (1 spp per line)", required=True)
    parser.add_argument("-m", "--maps", help="directory containing any map gene->protein id mappping", default='')
    args = parser.parse_args()
    ofdr = args.orthofinder
    spp  = args.spp
    mdir = args.maps

    assert os.path.isfile(ofdr), f"Could not locate file {ofdr}"
    assert os.path.isfile(spp),  f"Could not locate file {spp}"

    if (mdir != ''):
        assert os.path.isdir(mdir), f"Could not locate directory {mdir}"

    return (ofdr, spp, mdir)

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


def reformat_finder(ofdr: str, gpMap: dict, sppList: list) -> None:
    """reformat each line & write out afterwards"""

    outfile = "Reformat." + os.path.basename(ofdr)
    ofh     = open(outfile, 'w')
    fh      = gzip.open(ofdr, "rt") if ofdr.endswith(".gz") else open(ofdr, 'r')

    header  = ["#GroupID"] + sppList
    header  = '\t'.join(header) + '\n'
    ofh.write(header)

    split_ID = False
    COLUMN   = 3

    for linenum, line in enumerate(fh):
        if (linenum == 0): 
            continue
        fields  = line.strip().split('\t')
        row     = [fields[0]]
        spp_dic = {s : [] for s in sppList}
        count   = 0
        for i in range(COLUMN, len(fields)): # use 3 if using N0.tsv
            subfield = fields[i].split(", ")
            spp      = sppList[i - COLUMN]
            for isoform in subfield:
                if (isoform == ''): continue
                gene = ''
                if (split_ID):
                    isoform = isoform.split('.')[0]
                if (isoform.startswith("mRNA.")):
                    isoform = isoform.replace("mRNA.", '')
                if (spp in gpMap and isoform in gpMap[spp]):
                    gene = gpMap[spp][isoform]
                else:
                    tmp = isoform.split('.')[0]
                    if (spp in gpMap and tmp in gpMap[spp]):
                        gene = gpMap[spp][tmp]
                if (gene != ''):
                    spp_dic[spp].append(gene)
                    count += 1
        if (count == 0):
            continue
        for spp in sppList:
            row.append(', '.join(spp_dic[spp]) if spp_dic[spp] else '')
        outline = '\t'.join(row) + '\n'
        ofh.write(outline)

    fh.close()
    ofh.close()

def main() -> int:
    """entry point to this little subprogram"""

    # get arguments
    ofdr, spp, mdir = get_arguments()

    # get the files
    gpMap   = get_maps(mdir)
    sppList = get_spp(spp)

    reformat_finder(ofdr, gpMap, sppList)

    return 0

if __name__ == "__main__":
    main()
