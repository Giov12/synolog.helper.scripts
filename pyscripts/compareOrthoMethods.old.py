#!/bin/env python3

import argparse
import os
import glob
import sys
import gzip
from collections import defaultdict

"""
go through every gene and if shared is empty, skip
"""

#
# a global list that will store the indices per spp
#
sppKey = list()

class Gene:
    def __init__(self) -> None:
        """
        struct:
                {[set1, set2, etc..], [set1, set2, etc..], [set1, set2, etc..]}
        """
        self.orthoKey = {i : defaultdict(set) for i in range(len(sppKey))}
        self.group    = {i : set() for i in range(len(sppKey))}
        self.shared   = {s : set() for s in sppKey}

    def add_entry(self, software: int, spp: str, entries: set) -> None:

        if (len(entries) == 0):
            return

        for e in entries:
            self.orthoKey[software][spp].add(e)
            mem = spp + '|' + e
            self.group[software].add(mem)

    def in_a_group(self, genes : set) -> bool:

        for group in self.group.values():
            remaining = genes.difference(group)
            if (len(remaining) == 0):
                return True
        return False

    def compareSoftware(self, visited: dict) -> None:

        OFgroup = self.orthoKey[0]
        IPgroup = self.orthoKey[1]
        OMgroup = self.orthoKey[2]

        for spp in sppKey:
            OFgenes = OFgroup[spp]
            IPgenes = IPgroup[spp]
            OMgenes = OMgroup[spp]
            SHgenes = set() # shared genes

            # for genes in [OFgenes, IPgenes, OMgenes]:
            #     if (len(SHgenes) == 0) and (len(genes) > 0):
            #         SHgenes = genes
            #     elif (len(genes) > 0):
            #         SHgenes = SHgenes.intersection(genes)

            # get sizes for readability
            OFcnt = len(OFgenes)
            IPcnt = len(IPgenes)
            OMcnt = len(OMgenes)

            if (OFcnt > 0) and (IPcnt > 0) and (OMcnt > 0):
                SHgenes = OFgenes.intersection(IPgenes)
                SHgenes = SHgenes.intersection(OMgenes)
            elif (OFcnt > 0) and (IPcnt > 0):
                SHgenes = OFgenes.intersection(IPgenes)
            elif (OFcnt > 0) and (OMcnt > 0):
                SHgenes = OFgenes.intersection(OMgenes)
            elif (IPcnt > 0) and (OMcnt > 0):
                SHgenes = IPgenes.intersection(OMgenes)
            elif (OFcnt > 0):
                SHgenes = OFgenes
            elif (IPcnt > 0):
                SHgenes = IPgenes
            elif (OMcnt > 0):
                SHgenes = OMgenes
            
            # mark as visited
            for genes in [OFgenes, IPgenes, OMgenes]:
                for gene in genes:
                    visited[spp].add(gene)
            # for gene in SHgenes:
            #     visited[spp].add(gene)

            self.shared[spp] = SHgenes

    def write_group(self, grpID: int, gmap: dict, ogrp: set, outFH) -> None:
            for spp, genes in self.shared.items():
                if (len(genes) == 0): continue
                cnt   = len(genes)
                genes = list(genes)
                geout = ','.join(genes) # genes out
                line  = f"{grpID}\t{spp}\t{geout}\t{cnt}\n"
                outFH.write(line)
                for gene in genes:
                    mem = spp + '|' + gene
                    ogrp.add(mem)
                    gmap[mem] = grpID


    def __str__(self) -> str:
        s = f"""
             OrthoFinder:
                {sppKey[0]}: {self.orthoKey[0][sppKey[0]]} ({len(self.orthoKey[0][sppKey[0]])} genes)
                {sppKey[1]}: {self.orthoKey[0][sppKey[1]]} ({len(self.orthoKey[0][sppKey[1]])} genes)
                {sppKey[2]}: {self.orthoKey[0][sppKey[2]]} ({len(self.orthoKey[0][sppKey[2]])} genes)

            InParanoid:
                {sppKey[0]}: {self.orthoKey[1][sppKey[0]]} ({len(self.orthoKey[1][sppKey[0]])} genes)
                {sppKey[1]}: {self.orthoKey[1][sppKey[1]]} ({len(self.orthoKey[1][sppKey[1]])} genes)
                {sppKey[2]}: {self.orthoKey[1][sppKey[2]]} ({len(self.orthoKey[1][sppKey[2]])} genes)

            OrthoMCL:
                {sppKey[0]}: {self.orthoKey[2][sppKey[0]]} ({len(self.orthoKey[2][sppKey[0]])} genes)
                {sppKey[1]}: {self.orthoKey[2][sppKey[1]]} ({len(self.orthoKey[2][sppKey[1]])} genes)
                {sppKey[2]}: {self.orthoKey[2][sppKey[2]]} ({len(self.orthoKey[2][sppKey[2]])} genes)
            """
        return s
        
def get_arguments() -> tuple:
    """get the arguments"""

    """
        software key:
            OrthoFinder: 0
            InParanoid:  1
            OrthoMCL:    2
    """

    d = "Some python code to unify the results of OrthoFinder, InParanoid, & OrthoMCL"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-F", "--orthofinder", help="output N0.tsv from orthofinder", required=True)
    parser.add_argument("-I", "--inparanoid",  help="directory containing all pairwise output tables from inparnoid", required=True)
    parser.add_argument("-M", "--orthoMCL",    help="output file of orthologs from the OrthoMCL pipeline", required=True)
    parser.add_argument("-O", "--orthologs",    help="orthologs.tsv file from Synolog", required=True)

    args       = parser.parse_args()
    ofinderOut = args.orthofinder
    inParDir   = args.inparanoid
    oMClOut    = args.orthoMCL
    orthotsv   = args.orthologs

    assert os.path.isfile(ofinderOut), f"Could not locate file: {ofinderOut}"
    assert os.path.isfile(oMClOut),    f"Could not locate file: {oMClOut}"
    assert os.path.isfile(orthotsv),   f"Could not locate file: {orthotsv}"
    assert os.path.isdir(inParDir),    f"Could not locate directory: {inParDir}"
    
    return (ofinderOut, inParDir, oMClOut, orthotsv)

def fill_key(ofinderOut: str) -> None:
    """file sppKey with the species names"""

    fh   = gzip(ofinderOut, "rt") if (ofinderOut.endswith(".gz")) else open(ofinderOut, 'r')
    cols = fh.readline().strip().split('\t')[3:]
    cols.sort()

    for c in cols:
        sppKey.append(c)

    fh.close()

def strip_transcript_tag(geneID: str) -> str:
    """remove transcript ids from gene ids"""

    i = len(geneID) - 1

    while (i > 0):
        while (geneID[i].isdigit()):
            i -= 1
        if (geneID[i].lower() != 't'): return geneID
        if (i - 1 > 0 and (geneID[i - 1] == '.')):
            return geneID[:i - 1]
        return geneID
    return geneID

def parse_inparanoid_files(inParDir: str, Genomes: list) -> None:
    """scan the contents of the inparanoid output directory for the pairwise outfiles"""
    
    if (inParDir[-1] == '/'):
        inParDir = inParDir[:-1]

    tablesList = glob.glob(f"{inParDir}/table.*")

    if (len(tablesList) == 0):
        m = f"Could not locate any table.* files in {inParDir}"
        sys.exit(m)
    elif (len(tablesList) != len(sppKey)):
        m = f"Number of tables.* files does not match the number of expected ({len(sppKey)}) species. Check {inParDir}"
        sys.exit(m)

    for tableFile in tablesList:
        parse_inparanoid_table(tableFile, Genomes)

def parse_inparanoid_table(tableFile: str, Genomes: dict) -> None:
    """parse the table-formatted result of inparanoid"""

    bname = os.path.basename(tableFile)
    comps = bname.replace("table.", '').split('-')
    sidx  = 1
    orgs  = []

    for org in comps:
        if (org.endswith(".aa")):
            org = org.replace(".aa", '')
        orgs.append(org)

    # now to actually go through the file
    fh = gzip.open(tableFile, "rt") if (tableFile.endswith(".gz")) else open(tableFile, 'r')

    for line in fh:
        if (line[0].isdigit() == False):
            continue
        fields = line.strip().split('\t')[2:]
        entries = [set()] * 2
        for i in range(2):
            subfields = fields[i].strip().split(' ')
            for entry in subfields:
                if (entry == '') or (entry[0].isdigit()): continue
                gene = strip_transcript_tag(entry)
                entries[i].add(gene)
        # add entries to self & other
        for i in range(2):
            j      = 1 if (i == 0) else 0
            orgA   = orgs[i]
            orgB   = orgs[j]
            genes  = entries[i]
            ogenes = entries[j]
            for gene in genes:
                Genomes[orgA][gene].add_entry(sidx, orgA, genes)
                Genomes[orgA][gene].add_entry(sidx, orgB, ogenes)

        
    fh.close()

def parse_OrthoFinder(ofinderOut: str, Genomes: dict) -> None:
    """add the ortholog groups from OrthoFinder"""

    sidx = 0 # software index
    fh   = gzip.open(ofinderOut, "rt") if (ofinderOut.endswith(".gz")) else open(ofinderOut, 'r')

    for line in fh:
        if (line[0] != 'N'):
            header = line.strip().split('\t')
            orgs   = header[3:]
            continue

        cols    = line.strip().split('\t')[3:]
        entries = defaultdict(set)

        for i in range(len(cols)):
            col = cols[i]
            if (col == ''): continue
            org   = orgs[i]
            genes = col.split(',')
            for g in genes:
                g = g.strip(' ')
                g = strip_transcript_tag(g)
                entries[org].add(g)

        for i in range(3):
            if (i == 0):
                j, k = 1, 2
            elif (i == 1):
                j, k = 0, 2
            else:
                j, k = 0, 1
            sppA     = orgs[i]
            sppB     = orgs[j]
            sppC     = orgs[k]
            for g in entries[sppA]:
                Genomes[sppA][g].add_entry(sidx, sppA, entries[sppA])
                Genomes[sppA][g].add_entry(sidx, sppB, entries[sppB])
                Genomes[sppA][g].add_entry(sidx, sppC, entries[sppC])

    fh.close()

def parse_orthoMCL(oMCLOut: str, Genomes: dict) -> None:
    """parse the space delimited OrthoMCL output file"""

    sidx = 2
    fh   = gzip.open(oMCLOut, "rt") if (oMCLOut.endswith(".gz")) else open(oMCLOut, 'r')

    for line in fh:
        fields  = line.strip().split(' ')
        fields  = fields[1:] # skip ortho group id1
        entries = defaultdict(set)
        for field in fields:
            subfields = field.split('|')
            org       = subfields[0]
            gene      = strip_transcript_tag(subfields[1])
            entries[org].add(gene)

        for org, genes in entries.items():
            for gene in genes:
                Genomes[org][gene].add_entry(sidx, org, genes)
                for orgB, genesB in entries.items():
                    if (org == orgB): continue
                    Genomes[org][gene].add_entry(sidx, orgB, genesB)
                
    fh.close()

def compare_softwares(genomes: dict, gmap: dict, ogrps: list) -> None:
    """find only the shared genes in ortholog groups"""

    visited = {spp: set() for spp in sppKey}

    for spp, genesDict in genomes.items():
        for gname, gene in genesDict.items():
            if (gname in visited[spp]): continue
            gene.compareSoftware(visited)

    outFH  = open("CommonOrthoGroups.tsv", 'w')
    grpcnt = 0

    for spp, genesDict in genomes.items():
        for gname, gene in genesDict.items():
            if (len(gene.shared) > 0):
                ogrp = set()
                gene.write_group(grpcnt, gmap, ogrp, outFH)
                ogrps.append(ogrp)
                grpcnt += 1         

    outFH.close()

def load_orthologsTSV(orthotsv: str, OrthoGroups: list) -> None:
    """read in the orthologs.tsv with duplicates into Ortholog groups"""

    fh = gzip.open(orthotsv, "rt") if (orthotsv.endswith(".gz")) else open(orthotsv, 'r')

    curId  = 0
    curgrp = set()

    for line in fh:
        if (line[0] == '#'): continue
        fields = line.strip('\n').split('\t')
        spp    = fields[4].split('.')[0]
        gene   = fields[5]
        mem    = spp + '|' + gene
        oID    = int(fields[0])
        if (oID == curId):
            curgrp.add(mem)
        else:
            OrthoGroups.append(curgrp)
            curgrp = set()
            curId  = oID
            curgrp.add(mem)

    if (len(curgrp) > 0): OrthoGroups.append(curgrp)

    fh.close()

def compareOrthoGroups(gmap: dict, genomes: dict, ogrps: list, orthoGroups: list) -> None:
    """compare synolog orthogroups to other software"""

    fh = open("Differences.txt", 'w')
    gd = 0

    """
    Need to find where synolog is adding genes not present in shared by other software

    from the one's that are not shared, does any other software have them?
    """

    for i, og in enumerate(orthoGroups):
        idx = None 
        cnt = len(og)
        dis = list() # disagrees
        for mem in og:
            if (mem in gmap):
                if (idx == None):
                    idx  = gmap[mem]
                    cnt -= 1
                elif (gmap[mem] == idx):
                    cnt -= 1
                else:
                    dis.append(mem)
            else:
                dis.append(mem)
        if (cnt != 0):
            diff = ','.join(dis)
            line = str(i) + '\t' + diff + '\n'
            fh.write(line)
        else:
            gd += 1

    th = open("Testing.txt", 'w')
    for i, og in enumerate(orthoGroups):
        for mem in og:
            if (mem in gmap):
                omems = ogrps[gmap[mem]]
                smems = og.difference(omems)
                if (len(smems) > 0):
                    spp, gen = mem.split('|')
                    if (genomes[spp][gen].in_a_group(og)):
                        continue
                    smems = ','.join([m for m in smems])
                    line  = str(i) + '\t' + smems + '\n'
                    th.write(line)
                break

    fh.close()
    th.close()

    print("Number of good groups:", gd)

def main() -> int:
    """Parse and merge results across ortholog-calling software"""

    # get arguments
    ofinderOut, inParDir, oMClOut, orthotsv = get_arguments()

    # get the focal species
    fill_key(ofinderOut)

    # genomes will be represented as dicts of genes
    genomes = {s : defaultdict(Gene) for s in sppKey}

    # now parse the orthofinder file
    parse_OrthoFinder(ofinderOut, genomes)

    # parse the pairwise inparanoid tables
    parse_inparanoid_files(inParDir, genomes)

    # lastly, parse the orthomcl output
    parse_orthoMCL(oMClOut, genomes)

    # compare the softwares & map out the genes to their groups
    gmap  = dict()
    ogrps = list()
    compare_softwares(genomes, gmap, ogrps)
    print(len(ogrps))

    # now set collect the ortholog groups from synolog
    orthoGroups = list()
    load_orthologsTSV(orthotsv, orthoGroups)
    print(len(orthoGroups))

    # now compare
    compareOrthoGroups(gmap, genomes, ogrps, orthoGroups)

    # print(genomes["ceso"]["g17983"])
    # print('=' * 50)
    # print(genomes["ceso"]["g17988"])

    return 0

if __name__ == "__main__":
    main()
