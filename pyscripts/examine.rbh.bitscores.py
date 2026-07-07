#!/bin/env python3

import argparse
import os
import gzip
import math
from   collections import defaultdict
import matplotlib.pyplot as plt
from   scipy import stats

def get_arguments() -> tuple[str, str]:
    """get the arguments"""

    d = "Plot a linear regression of the bitscore and length products b/t two blast hits"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-b", "--btable", help="product.rbhs.bitscore.tsv file", required=True)
    parser.add_argument("-d", "--dups", help="orthologs.tsv file with duplicates", required=False, default='')
    args  = parser.parse_args()
    bfile = args.btable
    dfile = args.dups

    assert os.path.isfile(bfile), f"Could not locate file {bfile}"
    if (dfile != ''): assert os.path.isfile(dfile), f"Could not locate file {dfile}"

    return (bfile, dfile)

def get_duplicates(dfile: str) -> dict[str, set[str]]:
    """return a dict storing the duplicated genes for each spp"""

    if (dfile == ''):
        return {}
    
    dups = defaultdict(set)
    fh   = gzip.open(dfile, "rt") if dfile.endswith(".gz") else open(dfile, 'r')

    for line in fh:
        if (line[0] == '#'): continue
        fields = line.split('\t')
        if (fields[8] != 'D'): continue
        spp = fields[4]
        ge  = fields[5]
        dups[spp].add(ge)

    fh.close()

    return dups

def parse_bfile(bfile: str, dups: dict[str, set[str]]) -> dict[tuple[str, str], list[list[int], list[int], list[str]]]:
    """return a dictionary containing the bitscores & product length"""
    bpMap = {}
    fh    = gzip.open(bfile, "rt") if bfile.endswith(".gz") else open(bfile, 'r')
    teal  = "#008080"
    red   = "#FF2800"

    for line in fh:
        if (line[0] == '#'):
            continue
        fields = line.strip('\n').split('\t')
        spp1   = fields[0]
        spp2   = fields[3]
        g1     = fields[1]
        g2     = fields[4]
        l1     = int(fields[2]) # length of gene 1
        l2     = int(fields[5]) # length of gene 2
        pl     = math.log10(l1 * l2) # product length
        # colors will be either teal or red
        c = teal
        if (spp1 in dups and g1 in dups[spp1]) and (spp2 in dups and g2 in dups[spp2]):
            c = red
        # bitscr = float(fields[7])
        bitscr = math.log10(float(fields[6]))
        t      = (spp1, spp2)
        if (t not in bpMap):
            bpMap[t] = [[], [], []] # list of lists
        bpMap[t][0].append(bitscr)
        bpMap[t][1].append(pl)
        bpMap[t][2].append(c)

    fh.close()

    return bpMap

def make_linear_regression(bpMap: dict[tuple[str, str], list[list[int], list[int], list[str]]]) -> None:
    """calculate the linear regression and plot it"""
    
    fig, axes = plt.subplots(ncols=5, nrows=4, figsize=(24, 16))
    i         = 0
    ylab      = "Bitscore"
    xlab      = "Sequence Product Lengths"

    # positions in fig
    cols = [[0, 0], [0, 1], [0, 2], [0, 3], [0, 4], [1, 0], [1, 1], [1, 2], [1, 3], [1, 4],
            [2, 0], [2, 1], [2, 2], [2, 3], [2, 4], [3, 0], [3, 1], [3, 2], [3, 3], [3, 4]]

    for sppT, vals in bpMap.items():
        r, c       = cols[i][0], cols[i][1]
        spp1, spp2 = sppT[0], sppT[1]
        bits       = vals[0] # bitscores
        pls        = vals[1] # product lengths
        colors     = vals[2]   
        slope, intercept, rval, p, std_err = stats.linregress(pls, bits)
        preds      = [slope * x + intercept for x in pls]
        rval       = rval ** 2
        slope      = round(slope, 4)
        intercept  = round(intercept, 2)
        rval       = round(rval, 2)
        t          = f"Slope: {slope}\nIntercept: {intercept}\nr^2: {rval}"
        axes[r, c].set_title(f"{spp1} against {spp2}", fontsize = 16)
        axes[r, c].set_ylabel(ylab, fontsize = 14)
        axes[r, c].set_xlabel(xlab, fontsize = 14)
        axes[r, c].scatter(pls, bits, color = colors, s=10)
        axes[r, c].plot(pls, preds, "--k")
        axes[r, c].text(.01, 0.99, t, ha="left", va="top", transform=axes[r, c].transAxes, fontsize=14)
        i += 1

    plt.tight_layout()
    plt.savefig("UnCorrected.Bitscore.Product.png")

def main() -> int:
    """entry point to this little subprogram"""

    # get the product rbh file & dups if provided
    bfile, dfile = get_arguments()

    # get the duplicates if present
    dups = get_duplicates(dfile)

    # parse out the gene product lengths & corresponding bitscores
    bpMap = parse_bfile(bfile, dups)

    # now plot
    make_linear_regression(bpMap)

    return 0

if __name__ == "__main__":
    main()
