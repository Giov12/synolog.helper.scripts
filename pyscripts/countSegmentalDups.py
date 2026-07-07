#!/bin/env python3

import argparse
import os
import glob
import gzip
import matplotlib.pyplot as plt
from  collections import defaultdict

def get_arguments() -> tuple:
    """get the arguments"""

    d = "script to draw histogram of segmental duplication counts from synolog"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-d", "--dirc", help="directory containing *-SegmentalDuplications.tsv files from `synolog`", required=True)
    parser.add_argument("-t", "--target", help="target species", required=True)
    args = parser.parse_args()
    dirc = args.dirc
    spp  = args.target

    assert os.path.isdir(dirc), f"Could not locate directory {dirc}"

    return (dirc, spp)

def get_tsvs(dirc: str, spp: str) -> list:
    """get the segmental dups files for the target spp"""

    segFiles = list()

    if (dirc[-1] == '/'):
        segFiles = glob.glob(f"{dirc}*{spp}*-SegmentalDuplcations.tsv")
    else:
        segFiles = glob.glob(f"{dirc}/*{spp}*-SegmentalDuplcations.tsv")

    assert len(segFiles) > 0, f"Could not locate any *-SegmentalDuplcations.tsv files for {spp}"

    return segFiles

def make_array(freqDict: dict, vals: list) -> None:
    """place the frequencies & sizes into a 2D array"""

    for size, freq in freqDict.items():
        vals[0].append(size)
        vals[1].append(freq)

def other_spp(segFile: str, spp: str) -> str:
    """return the name of the other spp being analyzed"""

    bname  = os.path.basename(segFile)
    fields = bname.split('-')

    return fields[0] if fields[0] != spp else fields[1]


def read_segFile(segFile: str, spp: str, freqDict: dict) -> None:
    """fill the dict with the size of each segmental duplication"""

    fh   = gzip.open(read_segFile, "rt") if segFile.startswith(".gz") else open(segFile, 'r')
    prev = '-1'

    for line in fh:
        if (line[0] == '#'): continue
        fields = line.split('\t')
        if (fields[0] == prev) or (fields[2] != spp): continue
        freqDict[(int(fields[1]))] += 1
        prev = fields[0]

    fh.close()


def get_counts(dirc: str, spp: str) -> None:
    """loop through and plot the figures"""

    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 10))
    rowCols   = [(0, 0), (0, 1), (1, 0), (1, 1)]
    segFiles  = get_tsvs(dirc, spp)
    i         = 0

    for segFile in segFiles:
        freqDict = defaultdict(int)
        read_segFile(segFile, spp, freqDict)
        sppB = other_spp(segFile, spp)
        # create a 2D array with the values
        vals = [[], []] # x- & y-vals
        # now to plot
        make_array(freqDict, vals)
        n = sum(vals[1])
        r, c = rowCols[i][0], rowCols[i][1]
        axes[r, c].set_title(sppB, fontsize=14)
        axes[r, c].set_xlabel("Size (#anchoring genes)", fontsize=14)
        axes[r, c].set_ylabel("Frequency", fontsize=14)
        axes[r, c].bar(vals[0], vals[1], width=1.0, color="#008B8B", edgecolor="black")
        axes[r, c].set_ylim(0, 150)
        axes[r, c].text(0.95, 0.95, f"n = {n}", ha="right", va="top", transform=axes[r,c].transAxes, fontsize = 14)
        i += 1

    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    plt.savefig(f"{spp}-segDups.png")


def main() -> int:
    """entry point to this little subprogram"""

    # get arguments
    dirc, spp = get_arguments()

    # let get_counts() handle everything else
    get_counts(dirc, spp)

    return 0

if __name__ == "__main__":
    main()
