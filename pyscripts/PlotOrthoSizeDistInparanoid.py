#!/bin/env python3

import argparse
import os
import glob
import gzip
import matplotlib.pyplot as plt
from  collections import defaultdict
from  math import log

class Stats:
    def __init__(self, spp1: str, spp2: str) -> None:
        self.spp1  = spp1
        self.spp2  = spp2
        self.sizes = list()
        self.freqs = list()
        self.cnter = defaultdict(int)
        self.mx    = 0

    def add_size(self, size: int) -> None:
        self.cnter[size] += 1

    def make_arrays(self, ulg: bool) -> None:
        for s, f in self.cnter.items():
            self.sizes.append(s)
            if (ulg):
                f = log(f)
            self.freqs.append(f)
            self.mx = max(self.mx, f)

        # add spacing when plotting
        self.mx = (self.mx + 100)

    def clear_arrays(self) -> None:
        self.sizes.clear()
        self.freqs.clear()
        self.cnter.clear()

def get_arguments() -> tuple:
    """get the arguments"""

    d = "script to draw histogram of ortholog group size counts from the table* generated from inparanoid"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-d", "--dir",  help="directory to the table* files from inparanoid", required=True)
    parser.add_argument("-l", "--log", help="Plot log(Frequency) instead of abs(Frequency)", action="store_true")

    args = parser.parse_args()
    idir = args.dir
    ulg  = args.log # use log

    assert os.path.isdir(idir), f"Could not locate {idir}"

    return (idir, ulg)

def get_itables(idir: str) -> list:
    """return a list of each table* generated from inparanoid"""

    itables = list()

    if (idir[-1] == '/'):
        itables = glob.glob(f"{idir}table*-*")
    else:
        itables = glob.glob(f"{idir}/table*-*")

    assert len(itables) > 0, f"Could not locate any table* files at {idir}"

    return itables

def get_spp_name(itable: str) -> tuple:
    """get the species names from the file name"""

    bname = os.path.basename(itable)
    sects = bname.split('.')
    spp1  = sects[1]
    spp2  = sects[2].split('-')[1]

    return (spp1, spp2)

def count_genes(column: str) -> int:
    """count the number genes in an individual column"""

    ngenes = 0
    column = column.split() # remove whitespaces
    for i in range(0, len(column), 2):
        ngenes += 1

    return ngenes

def parse_itable(itable: str) -> Stats:
    """parse a table* file from inparanoid to get the size distribution of orthogroups"""

    spp1, spp2 = get_spp_name(itable)

    sobj = Stats(spp1, spp2)
    fh   = gzip.open(itable, "rt") if itable.endswith(".gz") else open(itable, 'r')

    for line in fh:
        if (line[0] == '#'): continue
        fields  = line.split('\t')
        ngenes  = 0
        column1 = fields[2].strip()
        column2 = fields[3].strip()
        ngenes += count_genes(column1)
        ngenes += count_genes(column2)
        sobj.add_size(ngenes)

    fh.close()

    return sobj

def parse_methods(itables: list[str]) -> list[Stats]:
    """return a list of Stats obj - one per comparison"""

    statObjs = []

    for itable in itables:
        statObjs.append(parse_itable(itable))

    return statObjs

def plot_methods(statsObjs: list[Stats], ulg: bool) -> None:
    """create a 2 x 5 figure of bar plots of orthogroup size frequencies"""
    
    if (ulg):
        ylim   = log(17_500)
        ytitle = "log(Frequency)"
    else:
        ylim   = 16_000
        ytitle = "Frequency"

    fig, axes = plt.subplots(nrows=2, ncols=5, figsize=(20, 12))
    cols = [[0, 0], [0, 1], [0, 2], [0, 3], [0,4],
            [1, 0], [1, 1], [1, 2], [1, 3], [1,4]]

    for i, sobj in enumerate(statsObjs):
        title = f"{sobj.spp1} vs {sobj.spp2}"
        r, c  = cols[i][0], cols[i][1]
        sobj.make_arrays(ulg)
        axes[r, c].set_title(title, fontsize = 16)
        axes[r, c].set_xlabel("OrthoGroup Size", fontsize = 14)
        axes[r, c].set_ylabel(ytitle, fontsize = 14)
        axes[r, c].bar(sobj.sizes, sobj.freqs, width = 1.0, color = "#008B8B", edgecolor = "black")
        axes[r, c].set_ylim(0, ylim) # prev was sobj.mx
        sobj.clear_arrays()

    plt.tight_layout()

    if (ulg):
        plt.savefig("InparanoidOrthoGroupSizesLogged.png")
    else:
        plt.savefig("InparanoidOrthoGroupSizes.png")

def main() -> int:
    """entry point to this little subprogram"""

    # get arguments
    idir, ulg = get_arguments()

    # get the input files
    itables = get_itables(idir)

    # get the counts
    statObjs = parse_methods(itables)

    # plot the counts
    plot_methods(statObjs, ulg)

    return 0

if __name__ == "__main__":
    main()
