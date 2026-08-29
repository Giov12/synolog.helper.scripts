#!/bin/env python3

import argparse
import os
import glob
import gzip
import sys

ofdr     = ''    # orthogroups.tsv
spp      = ''    # text file with new species column IDs
mdir     = ''    # directory with gpMap files
split_ID = False # split gene IDs to get first field in `.` separated string

def get_arguments() -> int:
    """get the arguments"""

    global ofdr, spp, mdir, split_ID

    d = "reformat the OrtoFinder output to a tsv for easier comparisons"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-o", "--orthofinder", help="orthogroups.tsv file", required=True)
    parser.add_argument("-s", "--spp", help="text file containing names of species for columns (1 spp per line)", required=True)
    parser.add_argument("-m", "--maps", help="directory containing any map gene->protein id mappping", default='')
    parser.add_argument("--split", help="split gene ID attributes if mRNA are present", action="store_true")
    args     = parser.parse_args()
    ofdr     = args.orthofinder
    spp      = args.spp
    mdir     = args.maps
    split_ID = args.split

    assert os.path.isfile(ofdr), f"Could not locate file {ofdr}"
    assert os.path.isfile(spp),  f"Could not locate file {spp}"

    if (mdir != ''):
        assert os.path.isdir(mdir), f"Could not locate directory {mdir}"

    return 0

def get_maps() -> dict[dict[str: str]]:
    """get any gene->protein mapping for easier comparison"""

    global mdir

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

def get_spp() -> list[str]:
    """return a list of all the focal spp"""

    global spp

    sppList = []

    fh = gzip.open(spp, "rt") if spp.endswith(".gz") else open(spp, 'r')

    for line in fh:
        sppList.append(line.strip())

    fh.close()

    return sppList

def reformat_finder(gpMap: dict[dict[str, str]], sppList: list[str]) -> int:
    """reformat each line & write out afterwards"""

    global ofdr, split_ID

    outfile = "Reformat." + os.path.basename(ofdr)
    ofh     = open(outfile, 'w')
    fh      = gzip.open(ofdr, "rt") if ofdr.endswith(".gz") else open(ofdr, 'r')
    COLUMN  = 1                       # orthofinder2 value was 3
    header  = ["#GroupID"] + sppList  # create the header
    header  = '\t'.join(header) + '\n'
    ofh.write(header)

    for linenum, line in enumerate(fh):
        if (linenum == 0): 
            continue
        fields  = line.strip().split('\t')
        if (fields[0] == "Orthogroup"):
            continue          # this is the header
        row     = [fields[0]] # row begins with orthogroup ID
        spp_dic = {s : [] for s in sppList}
        count   = 0
        for i in range(COLUMN, len(fields)):
            subfield = fields[i].split(", ")
            spp      = sppList[i - COLUMN]
            for isoform in subfield:
                if (isoform == ''): 
                    continue
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
                else:
                    msg = f"Unable to locate gene ID for {isoform}"
                    sys.exit(msg)
        if (count == 0): # skip empty rows
            continue
        for spp in sppList:
            row.append(', '.join(spp_dic[spp]) if spp_dic[spp] else '')
        outline = '\t'.join(row) + '\n'
        ofh.write(outline)

    fh.close()
    ofh.close()

    return 0

def main() -> int:
    """entry point to this little subprogram"""

    # get arguments
    get_arguments()

    # get the files
    gpMap   = get_maps()
    sppList = get_spp()

    reformat_finder(gpMap, sppList)

    return 0

if __name__ == "__main__":
    main()
