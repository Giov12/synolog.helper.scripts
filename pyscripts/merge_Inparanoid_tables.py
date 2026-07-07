#!/bin/env python3

import argparse
import os
import gzip
import glob
from collections import defaultdict

def get_arguments() -> tuple:
    """get the arguments"""

    d = "merge all the table*.faa files from inparanoid into a tsv where all the genes agree with one another"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-d", "--dir", help="directory containing all the table*.faa files", required=True)
    parser.add_argument("-m", "--maps", help="directory containing any map gene->protein id mappping", default='')
    args = parser.parse_args()
    idir = args.dir
    mdir = args.maps

    assert os.path.isdir(idir), f"Could not locate directory {idir}"

    if (mdir != ''):
        assert os.path.isdir(mdir), f"Could not locate directory {mdir}"

    return (idir, mdir)

def get_tables(idir: str) -> list:
    """return the list of input files"""

    ifiles = list()

    if (idir.endswith('/')):
        ifiles = glob.glob(f"{idir}table.*")
    else:
        ifiles = glob.glob(f"{idir}/table.*")

    assert len(ifiles) > 0, f"Could not locate any table.* files at {idir}"

    return ifiles

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

def get_spp_names(fname: str) -> tuple:
    """get the spp names from the file name"""

    bname = os.path.basename(fname)
    bname = bname.replace("table.", '')
    bname = bname.replace(".faa", '')
    spp   = bname.split('-')

    return (spp[0], spp[1])


def create_imap(ifiles: list, gpMap: dict) -> tuple:
    """go file by file, building upon the ortholog group"""

    # sort so we encounter the first spp >1 time
    ifiles.sort()

    imap = {} # iparanoid mapping
    taxa = []

    # next, read in all the data before we go in & merge
    for f in ifiles:
        spp  = get_spp_names(f)
        
        if (spp[0] not in imap):
            imap[spp[0]] = {}
            taxa.append(spp[0])
        if (spp[1] not in imap):
            imap[spp[1]] = {}
            taxa.append(spp[1])

        fh = open(f, 'r')
        for line in fh:
            if (line[0] == 'O'):
                continue
            fields = line.split('\t')
            genes1 = fields[2].strip()
            genes2 = fields[3].strip()
            fields = [genes1, genes2]
            add    = [[], []]

            for i, genes in enumerate(fields):
                genes = genes.split(' ')
                org  = spp[i]
                for g in genes:
                    if (g[0].isnumeric()): continue
                    g = g.split('.')[0]
                    if (org in gpMap):
                        g = gpMap[org][g]
                    if (g not in imap[org]):
                        imap[org][g] = defaultdict(set)
                    add[i].append(g)

            for i, genes in enumerate(add):
                j = 1 if (i == 0) else 0
                a = spp[i]
                b = spp[j]
                for g1 in genes:
                    for g2 in genes:
                        if (g1 != g2):
                            imap[a][g1][a].add(g2)
                    for g2 in add[j]:
                        imap[a][g1][b].add(g2)

        fh.close()

    return (imap, taxa)

def make_parent_map(imap: dict) -> dict:
    """create a mapping of a gene to its assigned parent"""

    parent = {}

    for spp, gene_map in imap:
        parent[spp] = {}
        for g in gene_map.keys():
            parent[spp][g] = g

    return parent

def examine_imap(imap: dict, taxa: list) -> None:
    """check the level of agreement per gene"""

    examined = {spp : set() for spp in taxa}
    parent   = make_parent_map(imap)

    NewImap  = {}

    for spp in taxa:
        genes_dict = imap[spp]
        for gene, gene_map in genes_dict.items():
            if (gene in examined[spp]):
                continue
            # get parent
            gene   = parent[spp][gene]
            newmap = {}
            for spp2 in taxa:
                if (spp2 not in gene_map):
                    continue
                spp2Genes = gene_map[spp2]
                for gene2 in spp2Genes:
                    if (gene2 in imap[spp2]):
                        gene2_map = imap[spp2][gene2]
                    else:
                        gene2_map = {}
                    


            examined[spp].add(gene)




def main() -> int:
    """get the gtf file and create a tsv file to serve as a map"""

    # get arguments
    idir, mdir = get_arguments()

    # get the files
    ifiles = get_tables(idir)
    gpMap  = get_maps(mdir)

    # start constructing
    imap, spp = create_imap(ifiles, gpMap)
    examine_imap(imap, spp)


    return 0

if __name__ == "__main__":
    main()
