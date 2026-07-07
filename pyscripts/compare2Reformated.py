#!/bin/env python3

import argparse
import os
import gzip
import glob
from   collections import defaultdict

"""
    compare the orthologs.tsv file with the reformated file of
    orthofinder or orthomcl
"""

class OrthoRes:
    def __init__(self) -> None:
        self.geneMap   = dict()
        self.orthogrps = list()
        self.syogMap   = dict()
        self.paralogs  = set()

    def get_name(self) -> str:
        return self.name
    
    def get_geneMap(self) -> dict[str, int]:
        return self.geneMap
    
    def get_orthogrps(self) -> list[set[str]]:
        return self.orthogrps
    
    def get_synoMap(self) -> dict[str, str]:
        return self.syogMap
    
    def get_paralogs(self) -> set[str]:
        return self.paralogs
    
    def parse_reformat_tsv(self, tsv: str, fltr: bool) -> None:
        """fill the geneMap dict and create a list of the ortholog groups"""

        sppList = list()
        fh      = gzip.open(tsv, "rt") if tsv.endswith(".gz") else open(tsv, 'r')
        pgrps   = 0 # groups of single species paralogs
        tot     = 0 # total orthogroups

        for line in fh:
            fields = line.strip().split('\t')
            if (line[0] == '#'):
                sppList = fields[1:]
                continue
            orthogroup = set()
            idx        = len(self.orthogrps)
            genes      = list()
            sppcnt     = 0
            tot       += 1
            for i in range(1, len(fields)):
                if (fields[i] == ''): continue
                sppcnt  += 1
                spp      = sppList[i - 1]
                subfield = fields[i].split(", ")
                for gene in subfield:
                    gene = f"{spp}:{gene}"
                    genes.append(gene)
                    # orthogroup.add(gene)
                    # self.geneMap[gene] = idx
            # if not filtering or at least 2 spp. present
            if ((fltr == False) or (sppcnt > 1)):
                for gene in genes:
                    orthogroup.add(gene)
                    self.geneMap[gene] = idx
            else:
                pgrps += 1
                for gene in genes:
                    self.paralogs.add(gene)
            
            self.orthogrps.append(orthogroup)
        fh.close()

        p = round(((pgrps / tot) * 100), 2)
        print(f"Number of single-species orthogroups filtered: {pgrps} ({p}%)\n")

    def load_synolog(self, osyn: str) -> None:
        """"this function does the actual reading and loading"""

        fh = gzip.open(osyn, "rt") if osyn.endswith(".gz") else open(osyn, 'r')

        prev = None
        og   = set()
        idx  = 0 # ortholog group id for class
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
                self.orthogrps.append(og)
                og   = set()
                idx += 1
                prev = gi
            og.add(gene)
            self.geneMap[gene] = idx
            self.syogMap[gene] = gi # gi's can skip nums
        
        # add last group
        self.orthogrps.append(og)

        fh.close()
        
def get_arguments() -> tuple[bool, str, str, str]:
    """get the arguments"""

    d = "Some python code to compare either orthomcl or orthofinder to synolog"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-f", "--filter",  help="filter for orthogroups with 2 or more species", action="store_true")
    parser.add_argument("-m", "--method",  help="reformated tsv of either orthofinder or orthomcl", required=True)
    parser.add_argument("-g", "--gtfs",    help="directory containing gtf files for focal species", required=True)
    parser.add_argument("-s", "--synolog", help="orthologs.tsv file from Synolog",                  required=True)

    args = parser.parse_args()
    fltr = args.filter
    omtd = args.method
    gdir = args.gtfs
    osyn = args.synolog

    assert os.path.isfile(omtd), f"Could not locate file: {omtd}"
    assert os.path.isfile(osyn), f"Could not locate file: {osyn}"
    assert os.path.isdir(gdir),  f"Could not locate directory: {gdir}"
    
    return (fltr, omtd, gdir, osyn)

def get_gtfs(gdir: str) -> list[str]:
    """helper function to get the gtf files in the provided directory"""

    gtfs = list()

    if (gdir[-1] == '/'):
        gtfs = glob.glob(f"{gdir}*.gtf.gz")
        if (len(gtfs) == 0):
            gtfs = glob.glob(f"{gdir}*.gtf")
    else:
        gtfs = glob.glob(f"{gdir}/*.gtf.gz")
        if (len(gtfs) == 0):
            gtfs = glob.glob(f"{gdir}/*.gtf")

    assert len(gtfs) > 0, f"Could not locate any gtf files in {gdir}"

    return gtfs

def get_gene_id(column: str) -> str:
    """parse the attributes column of a gtf to get the gene_id or transcript id"""

    fields = column.split(';')
    geid   = None

    if (len(fields) > 1):
        for field in fields:
            field = field.strip()
            field = field.split(' ')
            if (field[0] != "gene_id"): continue
            geid  = field[1].replace('"', '')
            break
    else:
        geid = fields[0].strip('\n')

    assert geid != None, f"Error processing attributes: {column}"

    return geid

def make_gene_map(gdir: str) -> dict[str, tuple[str, int]]:
    """create a mapping of each gene to its chr & index on that chromosome"""

    # get gtf files
    gtfs = get_gtfs(gdir)

    # NOTE: assumes gtf is properly sorted
    genemap = dict()

    for gtf in gtfs:
        spp = os.path.basename(gtf).split('.')[0]
        fh  = gzip.open(gtf, "rt") if gtf.endswith(".gz") else open(gtf, 'r')
        idx = 0 # index
        cur = None
        tot = 0

        for line in fh:
            if (line[0] == '#'): continue
            fields = line.split('\t')
            if (fields[2] != "gene"): continue
            Chr = fields[0]
            atr = fields[8] # attributes
            gID = get_gene_id(atr)
            gID = f"{spp}:{gID}"
            if (cur == None):
                cur = Chr
            if (Chr == cur):
                idx += 1
            else:
                idx = 1 # restart to new gene
                cur = Chr
            genemap[gID] = (Chr, idx)
            tot += 1

        fh.close()
        print(f"Loaded {tot} genes from {gtf}")

    return genemap


def parse_methods(omthd: str, osyn: str, fltr: bool) -> list[OrthoRes]:
    """return a list of each methods results in the OrthoRes data structure"""

    res = list()

    # first instance will be for the other method
    ORstruct = OrthoRes()
    ORstruct.parse_reformat_tsv(omthd, fltr)
    res.append(ORstruct)

    # second instance will be for synolog
    ORstruct = OrthoRes()
    ORstruct.load_synolog(osyn)
    res.append(ORstruct)
 
    return res

def make_gChr_map(ogrp: set[str], gene_map: dict[str, tuple[str, str]]) -> dict[str, set[str]]:
    """helper function to create a spp -> chrs map"""

    gChr_map = defaultdict(set)

    for gene in ogrp:
        spp = gene.split(':')[0]
        Chr = gene_map[gene][0]
        gChr_map[spp].add(Chr)

    return gChr_map

def decrement_count_spp(spp: str, remain: set[str]) -> int:
    """count the number of genes belonging to this spp"""

    # avoid repeated code below
    cnt = 0

    for ge in remain:
        if (ge.startswith(spp)):
            cnt += 1

    return cnt

def decrement_count_chr(spp: str, chrSet: set[str], remain: set[str], gene_map: dict[str, tuple[str, str]]) -> int:
    """count the number of genes belonging to this spp"""

    # get the chr from the set
    Chr = ''
    for c in chrSet:
        Chr = c
        break

    # avoid repeated code below
    cnt = 0

    for ge in remain:
        if (ge.startswith(spp)):
            oChr = gene_map[ge][0]
            if oChr == Chr:
                reason = "Same chrom"
            else:
                cnt += 1
            
    return cnt

def compare_subsets(sortho: set[str], remain: set[str], gene_map: dict[str, tuple[str, str]]) -> int:
    """check the number of chromosomes and dist between the two sets"""

    # let's create a mapping of spp -> chr
    sgChr_map = make_gChr_map(sortho, gene_map)
    rgChr_map = make_gChr_map(remain, gene_map)

    # is there a better way of keeping count instead of brute force?
    rcnt = len(remain)

    # check if we just don't have this spp in our orthogroup
    for spp in rgChr_map:
        if (spp not in sgChr_map):
            reason = "Missing species" # TODO need to figure this out
            rcnt  -= decrement_count_spp(spp, remain)

    if (rcnt == 0):
        return rcnt

    # for the spp that we have, check the chr
    for spp, chrSet in sgChr_map.items():
        if (spp not in rgChr_map): continue
        oChrSet = rgChr_map[spp]
        if (len(oChrSet) == 1):
            if (chrSet != oChrSet):
                reason = "Different chromosome"
                rcnt  -= decrement_count_spp(spp, remain)
        else:
            # let's check if our chr is found in their set of chrs
            if (chrSet.issubset(oChrSet)):
                reason = "Our chr is in here, need dist"
                rcnt  -= decrement_count_chr(spp, chrSet, remain, gene_map)
            else:
                # they have genes from multiple chrs that are not in our group
                reason = "Our chr is not in here, dist is inf"
                rcnt  -= decrement_count_spp(spp, remain)

    return rcnt

def compare_Subsets(sorthos: set[int], remain: set[str], gene_map: dict[str, tuple[str, str]], orthogroups: list[set[str]], synoMap: dict[str, str], fh) -> int:
    """compare the ge to the all the groups; TODO - this is n * (m*g) algorithm, could be improved"""

    for rem in remain:
        if (rem not in gene_map):
            fh.write(f"{rem} - Not in annotations\n")
            continue
        rChr, rIdx = gene_map[rem] # chr & index of the remaining ge
        spp        = rem.split(':')[0]
        sppFound   = False # assume spp not in any group
        sppDifChr  = True # assume different chr
        rdist      = float("inf") # will minimize dist
        reason     = ''
        ogids      = set()
        dist       = 15 # furthest a dup can be located
        for idx in sorthos:
            ogrp = orthogroups[idx]
            # now, loop through each of our groups and break early if need be
            for ge in ogrp:
                if (ge.startswith(spp)):
                    sppFound   = True
                    sChr, sIDX = gene_map[ge]
                    if (sChr == rChr):
                        sppDifChr = False
                        rdist = min(rdist, abs(sIDX - rIdx))
            ogids.add(synoMap[ge])
        if (sppFound == False):
            reason = "Species not found"
        elif (sppDifChr == True):
            reason = "On a different Chr"
        elif (rdist > dist):
            reason = f"Distance exceeds {dist} genes ({rChr} : {rdist} genes away)"
        else:
            reason = "Sequence similarity"
        # write the possible reason
        ogids   = list(ogids)
        outline = ','.join(ogids) + f" {rem} - {reason}\n"
        fh.write(outline)


def process_subsets(subsets: set[int], res: list[OrthoRes], gene_map: dict[str, tuple[str, str]]) -> None:
    """explain why the genes not in synolog orthogroup are in the other group"""

    # get the synolog orthogroups
    Sor   = res[1]
    ogrps = Sor.get_orthogrps()
    gmap  = Sor.get_geneMap()
    giMap = Sor.get_synoMap()

    # get the other method orthogroups
    Oor   = res[0]
    mgrps = Oor.get_orthogrps()
    mmap  = Oor.get_geneMap()

    # reasons for descrepancies
    micnt = 0 # missing & in a single group
    Smcnt = 0 # missing & in separate groups
    spcnt = 0 # all present in separate groups
    dfcnt = 0 # we have genes grouped where they have them in different groups
    otcnt = 0 # subsets that has a mixture of missing & splits in our and their software
    rcnt  = 0 # count of remaining genes per orthogroup of other method
    n     = len(ogrps)
    m     = len(subsets)
    nper  = round(((m / n) * 100), 2)
    print(f"Number of subsets evaluating: {m} ({nper}%)")
    fh    = open("TestingSubsets.txt", 'w')

    # determine the reason for these based on location
    dlist = list()

    for idx in subsets:
        ogrp = ogrps[idx]
        # get the other group
        for g in ogrp:
            if (g in mmap):
                jdx = mmap[g]
                break
        mgrp = mgrps[jdx]
        # get the members not in our group
        mems = mgrp.difference(ogrp)
        remn = set() # remaining genes
        oids = set([idx]) # other indexes
        miss = 0
        splt = 0

        # now bin the genes
        while (len(mems) > 0):

            tmp = list() # remove these genes
            
            for m in mems:
                if (m not in gmap):
                    miss += 1
                    tmp.append(m)
                    remn.add(m)
                    continue
                kdx  = gmap[m]
                kgrp = ogrps[kdx]
                for g2 in kgrp:
                    if (g2 in mems):
                        tmp.append(g2)
                    # we have a member in a different orthogroup
                    elif (g2 in mmap):
                        splt += 1
                        tmp.append(g2)
                        remn.add(g2)
                oids.add(kdx)
            
            for m in tmp:
                mems.discard(m)

        ##### end of while loop

        # flag to see if need to dig deeper
        flag = False
        # case 1, all genes are in separate groups
        if ((len(oids) > 1) and (miss == 0) and (splt == 0)):
            spcnt += 1
            dlist.append(oids)
        # case 2, we are missing some and our genes are in separate orthogroups
        elif ((len(oids) > 1) and (miss > 0) and (splt == 0)):
            flag   = True
            Smcnt += 1
            dlist.append(oids)
        # case 3, all genes are present in separate groups but the other method has them grouped in other orthogroups
        elif ((len(oids) > 1) and (miss == 0) and (splt > 0)):
            dfcnt += 1
        # case 4, all genes are present in a single group and we are just missing the rest
        elif ((len(oids) == 1) and (miss > 0) and (splt == 0)):
            flag   = True
            micnt += 1
            dlist.append(oids)
        # case 5, there is a mixture of splits in both softwares & we are missing some
        else:
            otcnt += 1
        
        # do we dig deeper?
        if (flag):
            rcnt += len(remn)
            compare_Subsets(oids, remn, gene_map, ogrps, giMap, fh)

    # calc percentages
    miper = round(((micnt / n) * 100), 2)
    Smper = round(((Smcnt / n) * 100), 2)
    spper = round(((spcnt / n) * 100), 2)
    dfper = round(((dfcnt / n) * 100), 2)
    otper = round(((otcnt / n) * 100), 2)

    print(f"Number of subsets where we are missing the remaining members: {micnt} ({miper}%)")
    print(f"Number of subsets where we are split and missing the rest: {Smcnt} ({Smper}%)")
    print(f"Number of subsets where we split the members: {spcnt} ({spper}%)")
    print(f"Number of subsets where we are split and they have them in different groups: {dfcnt} ({dfper}%)")
    print(f"Number of subsets where there is discrepancies in assignments and we are missing members: {otcnt} ({otper}%)")
    print(f"Number of remaining genes investigate: {rcnt}\n")

    # process the subsets that do not have a contradiction within groups
    process_dlist(dlist, res, gene_map)

    fh.close()


def process_superset(supersets: set[int], res: list[OrthoRes], gene_map: dict[str, tuple[str, str]]) -> None:
    """explain why the genes in the synolog orthogroup are not in the other group"""

    # get the synolog orthogroups
    Sor   = res[1]
    ogrps = Sor.get_orthogrps()
    gmap  = Sor.get_geneMap()
    imap  = Sor.get_synoMap()

    # get the other method orthogroups
    Oor   = res[0]
    mgrps = Oor.get_orthogrps()
    mmap  = Oor.get_geneMap()
    para  = Oor.get_paralogs()

    # reasons for descrepancies
    pcase = 0 # perfect case when they were just split
    Pcase = 0 # case where they have paralogs that we grouped and they kept it separate
    scase = 0 # case where synolog has extra
    mcase = 0 # case where other method has extra
    lcase = 0 # case where we have less but different genes remaining
    gcase = 0 # case where we have greater but different genes remaining
    ecase = 0 # case where we have an equal count of differing genes
    n     = len(mgrps)
    m     = len(supersets)
    nper  = round(((m / n) * 100), 2)
    print(f"Number of supersets from the other method : {m} ({nper}%)")
    fh    = open("Superset.txt", 'w')

    for idx in supersets:
        mgrp = mgrps[idx]
        # get our genes & their indices
        for g in mgrp:
            if (g in gmap):
                jdx = gmap[g]
                break
        sgrp = ogrps[jdx]
        # create a merged set of the other methods genes to see if they are split
        mdxs = set()
        for sg in sgrp:
            if (sg in mmap):
                mdxs.add(mmap[sg])
        mgrp = set()
        for mdx in mdxs:
            mgrp.update(mgrps[mdx])

        # get the members not in each other's group
        Mems = mgrp.difference(sgrp)
        Sems = sgrp.difference(mgrp)
        tmp  = set()
        c    = 0
    

        if ((len(Mems) == 0) and (len(Sems) == 0)):
            pcase += 1
        elif (len(Mems) == 0):
            scase += 1
            for g in Sems:
                is_para = (g in para)
                sidx    = imap[g]
                if (is_para):
                    c += 1
                    fh.write(f"{sidx} {g} - Paralog\n")
                elif (g not in mmap):
                    fh.write(f"{sidx} {g} - Extra (Missing)\n")
                else:     
                    fh.write(f"{sidx} {g} - Extra\n")
        elif (len(Sems) == 0):
            mcase += 1
        elif (len(Mems) > len(Sems)):
            lcase += 1
            for g in Sems:
                is_para = (g in para)
                sidx    = imap[g]
                if (is_para):
                    c += 1
                    fh.write(f"{sidx} {g} - Paralog\n")
                elif (g not in mmap):
                    fh.write(f"{sidx} {g} - Less (Missing)\n")
                else:    
                    fh.write(f"{sidx} {g} - Less\n")
        elif (len(Mems) < len(Sems)):
            gcase += 1
            for g in Sems:
                is_para = (g in para)
                sidx    = imap[g]
                if (is_para):
                    c += 1
                    fh.write(f"{sidx} {g} - Paralog\n")
                elif (g not in mmap):
                    fh.write(f"{sidx} {g} - Greater (Missing)\n")
                else:    
                    fh.write(f"{sidx} {g} - Greater\n")
        else:
            ecase += 1
            for g in Sems:
                is_para = (g in para)
                sidx    = imap[g]
                if (is_para):
                    c += 1
                    fh.write(f"{sidx} {g} - Paralog\n")
                elif (g not in mmap):
                    fh.write(f"{sidx} {g} - Equal (Missing)\n")
                else:    
                    fh.write(f"{sidx} {g} - Equal\n")
        if (c > 0):
            Pcase += 1

    # calc percentages
    pper = round(((pcase / m) * 100), 2)
    Pper = round(((Pcase / m) * 100), 2)
    sper = round(((scase / m) * 100), 2)
    mper = round(((mcase / m) * 100), 2)
    lper = round(((lcase / m) * 100), 2)
    gper = round(((gcase / m) * 100), 2)
    eper = round(((ecase / m) * 100), 2)

    print(f"Number of cases where the synolog orthogroup is split: {pcase} ({pper}%)")
    print(f"Number of cases where the synolog orthogroup grouped a gene they have in a paralog group: {Pcase} ({Pper}%)")
    print(f"Number of cases where the synolog orthogroup has extra genes and the other method is split: {scase} ({sper}%)")
    print(f"Number of cases where the synolog orthogroup has less genes than the combined orthogroups of the other method: {mcase} ({mper}%)")
    print(f"Number of cases where the synolog orthogroup has less (but >0) remaining genes than the combined orthogroups of the other method: {lcase} ({lper}%)")
    print(f"Number of cases where the synolog orthogroup has greater remaining genes than the combined orthogroups of the other method (>0 remaining): {gcase} ({gper}%)")
    print(f"Number of cases where the synolog orthogroup has an equal but differing set of remaining genes to the combined orthogroups of the other method: {ecase} ({eper}%)")

    fh.close()

    return None

def check_missing(missing: set[str], res: list[OrthoRes], genemap: dict[str, tuple[str, int]]) -> None:
    """try to come up with an explanation on why we are missing these genes"""

    fh    = open("MissingSynologGenes.txt", 'w')
    fh2   = open("missingSynologGenes_idxs.txt", 'w')
    # get the synolog orthogroups
    Sor   = res[1]
    ogrps = Sor.get_orthogrps()
    gmap  = Sor.get_geneMap()
    imap  = Sor.get_synoMap()

    # get the other method orthogroups
    Oor   = res[0]
    mgrps = Oor.get_orthogrps()
    mmap  = Oor.get_geneMap()

    for missG in missing:
        if (missG not in genemap):
            jidx = mmap[missG]
            fh.write(f"{missG} (Method group: {jidx}) - Not in annotations\n")
            continue
        mChr, mIdx = genemap[missG] # chr & index of the missing ge
        idxs       = set() # hold indices for our orthogroups
        jidx       = mmap[missG]
        mgrp       = mgrps[jidx]
        spp        = missG.split(':')[0]
        sppFound   = False # assume spp not in any group
        sppDifChr  = True # assume different chr
        rdist      = float("inf") # will minimize dist
        reason     = ''
        ogids      = set()
        dist       = 15 # furthest a dup can be located
        for mG in mgrp:
            if (mG in gmap):
                idx = gmap[mG]
                idxs.add(idx)
        # now loop through the indices we collected
        for idx in idxs:
            sgrp = ogrps[idx]
            for ge in sgrp:
                if (ge.startswith(spp) == False): continue
                sppFound   = True
                sChr, sIDX = genemap[ge]
                if (sChr == mChr):
                    sppDifChr = False
                    rdist = min(rdist, abs(sIDX - mIdx))
            ogids.add(imap[ge])
        if (sppFound == False):
            reason = "Species not found"
        elif (sppDifChr == True):
            reason = "On a different Chr"
        elif (rdist > dist):
            reason = f"Distance exceeds {dist} genes ({mChr} : {rdist} genes away)"
        else:
            reason = "Sequence similarity"
        # write the possible reason
        ogids   = list(ogids)
        outline = ','.join(ogids) + f" {missG} - {reason}\n"
        fh.write(outline)
        fh2.write(f"{jidx} {missG} {mChr} {mIdx}\n")
        

    fh.close()
    fh2.close()

def process_dlist(dlist: list[set[int]], res: list[OrthoRes], genemap: dict[str, tuple[str, int]]) -> None:
    """determine if it is synteny or sequence similarity that is causing these subsets"""

    # get the synolog orthogroups
    Sor   = res[1]
    ogrps = Sor.get_orthogrps()
    gmap  = Sor.get_geneMap()

    # get the other method orthogroups
    Oor   = res[0]
    mgrps = Oor.get_orthogrps()
    mmap  = Oor.get_geneMap()

    synr  = 0 # synteny is the reason
    seqr  = 0 # sequence similarity is the reason
    mbsn  = 0 # missing b/c synteny
    mbss  = 0 # missing b/c seq similarity

    for oids in dlist:
        gMap = defaultdict(list) # org -> [g1, g2, etc..]
        for oid in oids:
            ogrp = ogrps[oid]
            # only adding the ones present (missing will be dealt with later)
            for g in ogrp:
                spp = g.split(':')[0]
                pos = genemap[g]
                gMap[spp].append(pos)
        # sort by chr & then index
        for l in gMap.values():
            l.sort(key = lambda x: (x[0], x[1]))
        
        diffChrs = 0
        dist     = 0
        # now to check the localities of all the genes in oids
        for spp, l in gMap.items():
            if (len(l) > 1):  # single gene, no neighbors to compare to
                chrs = 1
                mdis = 0 # max distance
                prev = l[0][1]
                pchr = l[0][0]
                for i in range(1, len(l)):
                    if (pchr == l[i][0]):
                        mdis = max(mdis, l[i][1] - prev)
                    else:
                        chrs += 1
                        mdis  = 0 
                        pchr  = l[i][0]
                    prev = l[i][1]
                if (chrs > 1):
                    diffChrs += 1
                if (mdis > 150): # this is just a test
                    dist += 1
        if ((diffChrs > 1) or (dist > 1)):
            synr += 1
        else:
            seqr += 1 # Not sure how well founded this logic is
        
        # now to handle the missing genes
        mgenes = list()
        for oid in oids:
            ogrp = ogrps[oid]
            for g in ogrp:
                # get the other orthogroup
                if (g in mmap):
                    jdx  = mmap[g]
                    mgrp = mgrps[jdx]
                    break
            for g in mgrp:
                # not in synolog
                if (g not in gmap):
                    mgenes.append(g)
        
        # now to see where this gene falls within our list of genes
        for g in mgenes:
            if (g not in genemap):
                continue
            spp = g.split(':')[0]
            l   = gMap[spp]
            pos = genemap[g] # (chr, index)
            # see if we match any of the chrs
            if (len(l) > 1):
                i = 0
                f = False
                while (i < len(l)):
                    if (l[i][0] != pos[0]):
                        i += 1
                    # i is the first place we saw this chr, j will be the last
                    f = True
                    j = i
                    while ((j + 1 < len(l)) and (l[j + 1][0] == pos[0])):
                        j += 1
                    # get the dist from both ends
                    ldis = abs(l[i][1] - pos[1])
                    rdis = abs(l[j][1] - pos[1])
                    # just a test for now
                    if (min(ldis, rdis) > 150):
                        mbsn += 1 
                    else: # if it's close enough
                        mbss += 1
                    break
                if (f == False): # not found, so not in synteny block
                    mbsn += 1
            elif (len(l) == 1):
                if (l[0][0] != pos[0]):
                    dis = float("inf")
                else:
                    dis = abs(l[0][1] - pos[1])
                if (dis > 150):
                    mbsn += 1
                else:
                    mbss += 1

    print(f"Number of orthogroups split because of synteny: {synr}")
    print(f"Number of orthogroups split because of sequence similarity: {seqr}")
    print(f"Number of genes missing because of synteny: {mbsn}")
    print(f"Number of genes missing because of sequence similarity: {mbss}\n")


def examine_discrepancies(res: list[OrthoRes], subsyno: set[int], supsyno: set[int], genemap: dict[str, tuple[str, int]]) -> None:
    """examine why we have discrepancies with out sub & supersets"""

    # first, the subsets
    process_subsets(subsyno, res, genemap)

    # next, the supersets
    process_superset(supsyno, res, genemap)


def compare_methods(res: list[OrthoRes], genemap: dict[str, tuple[str, str]]) -> tuple[set[int], set[int]]:
    """"compare the synolog results to the other method & return the super- & sub-sets & missing genes found"""

    Oor = res[0] # other ortho method
    Sor = res[1] # synolog ortho method

    # other res
    mgrps = Oor.get_orthogrps()
    mmap  = Oor.get_geneMap()

    # synolog res
    ogrps = Sor.get_orthogrps()
    gmap  = Sor.get_geneMap()
    imap  = Sor.get_synoMap()
    
    # some synolog stats
    print("\nNumber of synolog orthologs", len(gmap), "across", len(ogrps), "orthogroups")
    fh  = open("NotInMethod.txt", 'w')
    cnt = 0

    for gene in gmap:
        if (gene not in mmap):
            fh.write(gene + '\n')
            cnt += 1
    print("Number of synolog specific orthologs:", cnt)
    fh.close()

    # other method stats
    print("Number of other orthologs", len(mmap), "across", len(mgrps), "orthogroups")
    
    fh  = open("NotInSynolog.txt", 'w')
    cnt = 0
    mgs = set() # missing genes
    for gene in mmap:
        if (gene not in gmap):
            fh.write(gene + '\n')
            mgs.add(gene)
            cnt += 1
    print("Number of other method specific orthologs:", cnt)
    fh.close()

    check_missing(mgs, res, genemap)

    #
    # j = indexes for other method ; i = indexes for synolog
    # j indexes will be used to print stats
    # i indexes will be used to examine & try to categorize the reasoning
    #       behind the differing orthogroups
    #
    splitcnt = 0
    eqlCnter = set() # stores i
    subsyno  = set() # stores i
    supcnter = set() # stores i
    supsyno  = set() # stores j
    visited  = [-1] * len(ogrps)

    # collection of file handlers
    
    fh2 = open("Subset.txt", 'w') # genes in subset category
    fh4 = open("Equal.txt", 'w') # genes in equal orthogroup category
    fh5 = open("Diff.txt", 'w') # genes in orthogroups of same size but different members
    fh6 = open("Split.txt", 'w') # genes that were binned in separate orthogroups

    for gene, i in gmap.items():
        if (visited[i] != -1): continue
        visited[i] = 1
        ogrp = ogrps[i]
        for g in ogrp:
            if (g not in mmap):
                continue
            ### continued logic
            j     = mmap[g]
            ogrp2 = mgrps[j]
            # check subset approach
            if (len(ogrp2) > len(ogrp)):
                if (ogrp.issubset(ogrp2)):
                    subsyno.add(i)
                    for sg in ogrp:
                        fh2.write(sg + '\n')
                    break
             # check superset
            elif (len(ogrp2) < len(ogrp)):
                if (ogrp2.issubset(ogrp)):
                    supcnter.add(i)
                    supsyno.add(j)
                    break                                               
            elif (ogrp2 == ogrp):
                eqlCnter.add(i)
                for sg in ogrp:
                    fh4.write(sg + '\n')
                break
            else:
                fh5.write(g + '\n')
                # check for left overs / remainders
            # should be call cases

        oids = set() # other ids
        for g in ogrp:
            if (g in mmap):
                oids.add(mmap[g])

        if (len(oids) > 1):
            splitcnt += 1
            for g in ogrp:
                fh6.write(g + '\n')

    eper = round(((len(eqlCnter) / len(ogrps)) * 100), 2)
    sper = round(((len(subsyno) / len(ogrps)) * 100), 2)
    Sper = round(((len(supcnter) / len(ogrps)) * 100), 2)
    tper = round(((splitcnt / len(ogrps)) * 100), 2)
    
    print("\nCounts and percentages of synolog orthogroups")
    print("Equal count:", len(eqlCnter), f"({eper}%)")
    print("Subset count:", len(subsyno), f"({sper}%)")
    print("Superset count:", len(supcnter), f"({Sper}%)")
    print("Number of orthogroups that are split:", splitcnt, f"({tper}%)\n")

    for f in [fh2, fh4, fh5, fh6, fh6]:
        f.close()

    return (subsyno, supsyno)

def main() -> int:
    """Compare the orthologs.tsv file to either the orthofinder or orthomcl software"""

    # get arguments
    fltr, omtd, gdir, osyn = get_arguments()

    # load the results
    res = parse_methods(omtd, osyn, fltr)

    # load the genemap
    genemap = make_gene_map(gdir)

    # compare this method to synolog
    subsyno, supsyno = compare_methods(res, genemap)

    # tease out the discrepancies
    examine_discrepancies(res, subsyno, supsyno, genemap)

    return 0

if __name__ == "__main__":
    main()
