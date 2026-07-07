#!/bin/env python3

import argparse
import os
import gzip
import glob

"""
    compare the orthologs.tsv file with the table files 
    from inparanoid
"""

class OrthoRes:
    def __init__(self) -> None:
        self.geneMap   = dict()
        self.orthogrps = list()
    
    def get_geneMap(self) -> dict:
        return self.geneMap
    
    def get_orthogrps(self) -> list:
        return self.orthogrps
    
    def add_geneMap(self, geneMap: dict) -> None:
        self.geneMap = geneMap

    def add_orthogrps(self, ogrps: list) -> None:
        self.orthogrps = ogrps
    
class FileManager:
    def __init__(self, spp1: str, spp2: str) -> None:
        self.spp1 = spp1
        self.spp2 = spp2
    
    def create_file_handles(self) -> None:
        self.fh1 = open(f"{self.spp1}-{self.spp2}.NotInMethod.txt", 'w')
        self.fh2 = open(f"{self.spp1}-{self.spp2}.Subset.txt", 'w')
        self.fh3 = open(f"{self.spp1}-{self.spp2}.Superset.txt", 'w')
        self.fh4 = open(f"{self.spp1}-{self.spp2}.Equal.txt", 'w')
        self.fh5 = open(f"{self.spp1}-{self.spp2}.Diff.txt", 'w')
        self.fh6= open(f"{self.spp1}-{self.spp2}.Split.txt", 'w')
        self.fh7 = open(f"{self.spp1}-{self.spp2}.NotInSynolog.txt", 'w')

    def close_file_handles(self) -> None:
        for f in [self.fh1, self.fh2, self.fh3, self.fh4, self.fh5, self.fh6, self.fh7]:
            f.close()
        
def get_arguments() -> tuple:
    """get the arguments"""

    d = "Some python code to compare the tables generated from inparanoid with synolog"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-d", "--dir",  help="directory to the table* files from inparanoid", required=True)
    parser.add_argument("-g", "--gmaps",  help="directory to the *gpMap.tsv files", required=True)
    parser.add_argument("-s", "--synolog", help="orthologs.tsv file from Synolog",  required=True)

    args = parser.parse_args()
    idir = args.dir
    gdir = args.gmaps
    osyn = args.synolog

    assert os.path.isdir(idir), f"Could not locate directory: {idir}"
    assert os.path.isdir(gdir), f"Could not locate directory: {gdir}"
    assert os.path.isfile(osyn), f"Could not locate file: {osyn}"
    
    return (idir, gdir, osyn)

def load_geneMap(gmap: str) -> dict:
    """load in the transcript to gene id map"""

    fh = gzip.open(gmap, "rt") if gmap.endswith(".gz") else open(gmap, 'r')

    tr2gene = {} # transcript to gene

    for line in fh:
        fields = line.strip().split('\t')
        gene   = fields[0]
        trncrp = fields[1]
        tr2gene[trncrp] = gene

    fh.close()

    return tr2gene

def get_itables(idir: str) -> list:
    """return a list of each table* generated from inparanoid"""

    itables = list()

    if (idir[-1] == '/'):
        itables = glob.glob(f"{idir}table*-*")
    else:
        itables = glob.glob(f"{idir}/table*-*")

    assert len(itables) > 0, f"Could not locate any table* files at {idir}"

    return itables

def get_tmaps(gdir: str) -> dict:
    """return a dict of dict for each file to get the transcript to gene conversion"""

    gmaps = list()

    if (gdir[-1] == '/'):
        gmaps = glob.glob(f"{gdir}*gpMap.tsv*")
    else:
        gmaps = glob.glob(f"{gdir}/*gpMap.tsv*")

    assert len(gmaps) > 0, f"Could not locate any table* files at {gdir}"

    tmap = dict()

    for gmap in gmaps:
        spp = os.path.basename(gmap).split('.')[0]
        tmap[spp] = load_geneMap(gmap)

    return tmap

def load_synolog(osyn: str) -> OrthoRes:
    """load the orthologs.tsv file into a OrthoRes structure"""

    # open the file (gzip compression supported)
    fh = gzip.open(osyn, "rt") if osyn.endswith(".gz") else open(osyn, 'r')

    prev = None
    og   = set()
    idx  = 0 # ortholog group id for class

    # structs to give OrthoRes
    orthogrps = list()
    geneMap   = dict()

    for line in fh:
        if (line[0] == '#'):
            continue
        fields = line.split('\t')
        gi     = fields[0] # group id
        spp    = fields[4].split('.')[0]
        ge     = fields[5]
        gene   = f"{spp}:{ge}"
        if (prev == None):
            prev = gi
        if (gi != prev):
            orthogrps.append(og)
            og   = set()
            idx += 1
            prev = gi
        og.add(gene)
        geneMap[gene] = idx
        
    # add last group
    orthogrps.append(og)
    fh.close()

    # now to initialize the obj and load the data
    orObj = OrthoRes()
    orObj.add_geneMap(geneMap)
    orObj.add_orthogrps(orthogrps)

    return orObj

def get_spp_name(itable: str) -> tuple:
    """get the species names from the file name"""

    bname = os.path.basename(itable)
    sects = bname.split('.')
    spp1  = sects[1]
    spp2  = sects[2].split('-')[1]

    return (spp1, spp2)

def get_genes_from_column(column: str, spp: str, tmaps: dict) -> list:
    """get the genes & attach spp name from an individual column"""

    genes  = []
    column = column.split() # remove whitespaces
    for i in range(0, len(column), 2):
        gene = column[i].split('.')[0]
        if (spp in tmaps):
            gene = tmaps[spp][gene]
        genes.append(f"{spp}:{gene}")

    return genes

def compare_synolog_to_inparanoid(orObj: OrthoRes, itables: list, tmaps: dict) -> None:
    """compare 1 by 1 the results from inparanoid to synolog"""

    # get the synolog orthologs
    ogrps = orObj.get_orthogrps()
    gmap  = orObj.get_geneMap()

    # now to loop through each file & compare row by row
    for itable in itables:
        
        # set up the metrics
        splitcnt = 0
        notinmtd = 0 # not in method
        eqlCnter = 0
        subcnter = 0
        supcnter = 0
        missCnt  = 0
        revisit  = 0
        ogrpinp  = 0
        sgrpcnt  = set()
        tot      = 0
        
        spp1, spp2 = get_spp_name(itable)
        
        # count the number of orthologs for the two spp
        stot = 0
        for gene in gmap:
            if ((gene.startswith(spp1)) or (gene.startswith(spp2))):
                stot += 1

        fh = gzip.open(itable, "rt") if itable.endswith(".gz") else open(itable, 'r')
        fm = FileManager(spp1, spp2)
        fm.create_file_handles() # create the file handles to write to

        for line in fh:
            if (line[0] == 'O'): continue
            ogrpinp += 1
            fields   = line.split('\t')
            genes1   = fields[2].strip()
            genes1   = get_genes_from_column(genes1, spp1, tmaps)
            genes2   = fields[3].strip()
            genes2   = get_genes_from_column(genes2, spp2, tmaps)
            genes1.extend(genes2)
            # collect the group nums for synolog ortho groups
            syngrp = set()
            mcnt   = 0 # match count
            tot   += len(genes1)
            for gene in genes1:
                if (gene not in gmap):
                    fm.fh7.write(gene + '\n')
                    missCnt += 1
                    continue
                idx = gmap[gene]
                syngrp.add(idx)
                sgrpcnt.add(idx)
                mcnt += 1
            if (mcnt == len(genes1)):
                # all genes belong to a single group
                if (len(syngrp) == 1):
                    eqlCnter += 1
                # all genes found but split across groups
                else:
                    splitcnt += 1
            # not all genes found in synolog
            else:
                # what genes were found at least belong to a single group
                if (len(syngrp) == 1):
                    subcnter += 1
                # not all genes are found and in multiple synolog ortho groups
                elif (mcnt > 0):
                    revisit += 1
        fh.close()
        fm.close_file_handles()

        eper = round(((eqlCnter / ogrpinp) * 100), 2)
        sper = round(((splitcnt / ogrpinp) * 100), 2)
        mper = round(((subcnter / ogrpinp) * 100), 2)
        Mper = round(((revisit / ogrpinp) * 100), 2)
        aper = round(((missCnt / tot) * 100), 2)


        print("Comparison between", spp1, "and", spp2)
        print("Number of orthogroups in synolog:", len(sgrpcnt))
        print("Number of orthologs in synolog:", stot)
        print("Number of orthogroups in inparanoid:", ogrpinp)
        print("Number of orthologs in inparanoid:", tot)
        print("Equal Count:", eqlCnter, f"({eper}%)")
        print("All found but across 2+ groups:", splitcnt, f"({sper}%)")
        print("Total number of genes in inparanoid but absent in synolog:", missCnt, f"({aper}%)")
        print("Missing some but remaining in a single group:", subcnter, f"({mper}%)")
        print("Missing some but remaining in 2+ groups:", revisit, f"({Mper}%)\n")

def main() -> int:
    """Entry point to compare orthologs.tsv to all table* from the inparanoid software"""

    # get arguments
    idir, gdir, osyn = get_arguments()

    # load synolog
    orObj = load_synolog(osyn)

    # get the inparanoid tables
    itables = get_itables(idir)

    # get the transcript to gene map
    tmaps = get_tmaps(gdir)

    # run the comparisons
    compare_synolog_to_inparanoid(orObj, itables, tmaps)

    return 0

if __name__ == "__main__":
    main()
