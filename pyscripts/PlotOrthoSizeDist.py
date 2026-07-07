#!/bin/env python3

import argparse
import os
import gzip
import matplotlib.pyplot as plt
from  collections import defaultdict
from  math import log

class Stats:
    def __init__(self) -> None:
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

    d = "script to draw histogram of ortholog group size counts from synolog, orthomcl, & orthofinder"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-s", "--synolog", help="orthologs.tsv file generated from synolog", required=True)
    parser.add_argument("-m", "--orthomcl", help="Reformated orthomcl final output", required=True)
    parser.add_argument("-f", "--orthofnd", help="Reformated orthofinder final output", required=True)
    parser.add_argument("-l", "--log", help="Plot log(Frequency) instead of abs(Frequency)", action="store_true")
    args = parser.parse_args()
    syn  = args.synolog
    mcl  = args.orthomcl
    ofn  = args.orthomcl
    ulg  = args.log # use log

    assert os.path.isfile(syn), f"Could not locate {syn}"
    assert os.path.isfile(mcl), f"Could not locate {mcl}"
    assert os.path.isfile(ofn), f"Could not locate {ofn}"

    return (syn, mcl, ofn, ulg)

def parse_reformated(rfile: str) -> Stats:
    """parse the reformated file to get the size distribution of orthogroups"""

    sobj = Stats()
    fh   = gzip.open(rfile, "rt") if rfile.endswith(".gz") else open(rfile, 'r')

    for line in fh:
        if (line[0] == '#'): continue
        fields = line.strip().split('\t')
        ngenes = 0
        for i in range(1, len(fields)):
            if (fields[i] == ''): continue
            genes   = fields[i].split(", ")
            ngenes += len(genes)
        sobj.add_size(ngenes)

    fh.close()

    return sobj

def parse_synolog(syn: str) -> Stats:
    """parse the orthologs.tsv file to get the size distribution of orthogroups"""

    sobj   = Stats()
    fh     = gzip.open(syn, "rt") if syn.endswith(".gz") else open(syn, 'r')
    cur    = None
    ngenes = 0

    for line in fh:
        if (line[0] == '#'): continue
        grp = line.split('\t')[0]
        if (cur == None):
            cur = grp
        if (cur == grp):
            ngenes += 1
        else:
            sobj.add_size(ngenes)
            ngenes = 1 # counting current gene
            cur    = grp

    # add last bit
    sobj.add_size(ngenes)

    fh.close()

    return sobj

def parse_methods(syn: str, mcl: str, ofn: str) -> list:
    """return a list of Stats obj - one per method"""

    # order Synolog, OrthoMCL, OrthoFinder
    statObjs = [None, None, None]

    statObjs[0] = parse_synolog(syn)
    statObjs[1] = parse_reformated(mcl)
    statObjs[2] = parse_reformated(ofn)

    return statObjs

def plot_methods(statsObjs: list[Stats], ulg: bool) -> None:
    """create a 1 x 3 figure of bar plots of orthogroup size frequencies"""

    titles = ["Synolog", "OrthoMCL", "OrthoFinder"]
    
    if (ulg):
        ylim   = log(8500)
        ytitle = "log(Frequency)"
    else:
        ylim   = 8000
        ytitle = "Frequency"

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(20, 6))

    for i, sobj in enumerate(statsObjs):
        sobj.make_arrays(ulg)
        axes[i].set_title(titles[i], fontsize = 16)
        axes[i].set_xlabel("OrthoGroup Size", fontsize = 14)
        axes[i].set_ylabel(ytitle, fontsize = 14)
        axes[i].bar(sobj.sizes, sobj.freqs, width = 1.0, color = "#008B8B", edgecolor = "black")
        axes[i].set_ylim(0, ylim) # prev was sobj.mx
        sobj.clear_arrays()

    plt.tick_params()
    if (ulg):
        plt.savefig("OrthoGroupSizesLogged.png")
    else:
        plt.savefig("OrthoGroupSizes.png")

def main() -> int:
    """entry point to this little subprogram"""

    # get arguments
    syn, mcl, ofn, ulg = get_arguments()

    # parse the input files
    statObjs = parse_methods(syn, mcl, ofn)

    # plot the counts
    plot_methods(statObjs, ulg)

    return 0

if __name__ == "__main__":
    main()
