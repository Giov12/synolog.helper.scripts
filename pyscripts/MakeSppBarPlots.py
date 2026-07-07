#!/bin/env python3

import matplotlib.pyplot as plt
from  math import log


def plot_methods() -> None:
    """hard coded function to plot precomputed ortholog counts per spp"""

    cesoPer = round((19489 / 20636) * 100, 2)
    cgunPer = round((20220 / 21659) * 100, 2)
    emacPer = round((19495 / 20633) * 100, 2)
    gacuPer = round((18499 / 20787) * 100, 2)
    locuPer = round((15697 / 18341) * 100, 2)

    full = False

    if (full):
        cesoOF = 95.90
        cgunOF = 95.90
        emacOF = 96.40
        gacuOF = 94.40
        locuOF = 96.90
        cesoOM = 92.57
        cgunOM = 92.41
        emacOM = 94.00
        gacuOM = 92.18
        locuOM = 88.73
    else:
        cesoOF = 95.01
        cgunOF = 94.92
        emacOF = 95.57
        gacuOF = 92.71
        locuOF = 91.64
        cesoOM = 91.96
        cgunOM = 91.27
        emacOM = 93.27
        gacuOM = 89.95
        locuOM = 85.75


    sppMap = {"ceso" : [["Total", "Synolog", "OrthoFinder", "OrthoMCL"], [100.0, cesoPer, cesoOF, cesoOM]],
              "cgun" : [["Total", "Synolog", "OrthoFinder", "OrthoMCL"], [100.0, cgunPer, cgunOF, cgunOM]],
              "emac" : [["Total", "Synolog", "OrthoFinder", "OrthoMCL"], [100.0, emacPer, emacOF, emacOM]],
              "gacu" : [["Total", "Synolog", "OrthoFinder", "OrthoMCL"], [100.0, gacuPer, gacuOF, gacuOM]],
              "locu" : [["Total", "Synolog", "OrthoFinder", "OrthoMCL"], [100.0, locuPer, locuOF, locuOM]]}

    colors = ["darkcyan", "steelblue", "cadetblue", "lightblue"]
    
    i = 0
    fig, axes = plt.subplots(nrows=1, ncols=5, figsize=(20, 6))

    for spp, orthoinfo in sppMap.items():
        titles = orthoinfo[0]
        percs  = orthoinfo[1]
        axes[i].set_title(spp, fontsize = 16)
        axes[i].set_xticklabels(titles, fontsize = 9)
        if (i == 0):
            axes[i].set_ylabel("Percent of genes (%)", fontsize = 14)
        axes[i].bar(titles, percs, width = 0.95, color = colors, edgecolor = "black")
        axes[i].set_ylim(0, 105)
        for j in range(len(percs)):
            axes[i].text(j, percs[j] + 0.5 , f"{percs[j]}%", ha = "center", fontsize = 12) 
        i += 1


    plt.tight_layout()
    plt.savefig("OrthoSppBarPlots.png")

def main() -> int:
    """entry point to this little subprogram"""

    # plot the counts
    plot_methods()

    return 0

if __name__ == "__main__":
    main()
