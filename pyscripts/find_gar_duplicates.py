#!/bin/env python3

import argparse
import os
import gzip
from collections import defaultdict, Counter
from typing import TextIO

# this code is custom made for identifying orthogroups
# where the group contains genes from spotted gar, zebrafish
# & medaka, and locates distant paralogs (presumed by OrthoDB v12)
# which may serve as distant ohnologs

class INFO:
    #
    # a simple container to hold each genes information
    # such as synolog copy state, chromosome, etc..
    #
    __slots__ = ("copy_state", "chrom", "orthodb_grp")
    def __init__(self, copy_state: str, chrom: str) -> None:
        self.copy_state  = copy_state
        self.chrom       = chrom
        self.orthodb_grp = '' # determine later

class GarGroup:
    #
    # a synolog orthogroup with exactly one gar gene and at least
    # one zebrafish and one medaka gene
    #
    __slots__ = ("id", "gar_gene", "zf_genes", "md_genes", "cands")
    
    def __init__(self, id_: str, gar_gene: str, zf_genes: list[str], md_genes: list[str]) -> None:
        self.id       = id_
        self.gar_gene = gar_gene
        self.zf_genes = zf_genes
        self.md_genes = md_genes
        self.cands    = list()

class OrthoDB:
    __slots__ = ("gar_counts", "synolog_grps")

    def __init__(self):
        self.gar_counts   = 0
        self.synolog_grps = defaultdict(list) # synolog group -> members

class Candidate:
    #
    # a synolog orthogroup whose zebrafish and medaka genes are
    # paralogs (same orthodb orthogroup) of a gar orthogroup's fish genes
    #
    __slots__ = ("id", "category", "one_gar", "zf_genes", "md_genes")

    def __init__(self, id_: str, category: str, one_gar: bool) -> None:
        self.id       = id_
        self.category = category # teleost-only, other species, no gar, or also has gar
        self.one_gar  = one_gar  # the shared orthodb orthogroup has exactly 1 gar gene
        self.zf_genes = list()   # zebrafish paralogs in this orthogroup
        self.md_genes = list()   # medaka paralogs in this orthogroup


synolog_file = ''
orthodb_file = ''
genes_info   = dict() # gene -> INFO
orthodb_ogs  = dict() # gene ->
synolog_grps = defaultdict(lambda: defaultdict(list)) # grp -> spp -> [genes]
odb_mapping  = defaultdict(OrthoDB)
gar_groups   = list()    # list of GarGroup objects
counts       = Counter() # for some stats
categ_counts = Counter() # category counts
categories   = ("teleost-only", "other species, no gar", "also has gar")

# constant globals
GAR                  = "locu"
MEDAKA               = "olat"
ZEBRAFISH            = "drer"
ONE_GAR_TWO_FISH     = "1 gar + >=1 zebrafish + >=1 medaka"
EXCLUDED_COPIES      = "excluded: extra fish copies"
GAR_ANCHORED         = "gar-anchored groups"
OFF_CHROMS           = "gar-anchored, but a fish gene off assembled chromosomes"
WITH_PARALOG_GROUP   = "...with a paralog group"
PARALOG_NO_GAR       = "...paralog group has no gar gene"
PARALOG_NO_GAR_ODB   = "...and OrthoDB says gar has 1 gene"
PARALOG_TELEOST_ONLY = "...and paralog group is zebrafish+medaka only"

def get_arguments() -> int:
    """get the arguments"""

    global orthodb_file, synolog_file

    desc  = "Using the renamed genes from OrthoDB v12, locate orthogroups that may have fish ohnologs based on spotted gar"
    shelp = "orthogroups.tsv file generated from synolog"
    ohelp = "A 4 column file generated with rename_odb_genes for a orthodb dataset"

    parser = argparse.ArgumentParser(description = desc)
    parser.add_argument("-s", "--synolog", help=shelp, required=True, type=str)
    parser.add_argument("-o", "--orthodb", help=ohelp, required=True, type=str)

    args          = parser.parse_args()
    synolog_file  = args.synolog
    orthodb_file  = args.orthodb

    assert os.path.isfile(synolog_file), f"Could not locate file: {synolog_file}"
    assert os.path.isfile(orthodb_file), f"Could not locate file: {orthodb_file}"

    return 0

def open_file(infile: str) -> TextIO:
    """a simple helper function to handle gzipped or flat files"""
    return gzip.open(infile, "rt") if infile.endswith(".gz") else open(infile, 'r')


def load_synolog() -> int:
    """create load the synolog orthogroups"""

    global synolog_file, genes_info, synolog_grps

    # synolog orthologs.tsv has chrom position
    
    fh  = open_file(synolog_file)
    tot = 0

    # header
    # [GroupID][Type] [ElementCount] [OrgCount]
    # [Org] [ElementId] [Name] [CopyNumState] [ParalogCount] 
    # [Chr] [StartBp] [EndBp]
    #

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.split('\t')
        assert len(fields) == 12, f"Malformed line encountered in {synolog_file}:\n{line}"
        tot  += 1
        grpID = fields[0]
        spp   = fields[4].lower()
        idx   = spp.find('.')
        if (idx != -1): # only keep org ID
            spp = spp[:idx]
        gene  = fields[5].lower()
        mem   = f"{spp}_{gene}"
        chrom = fields[9]
        cpyst = fields[7]
        genes_info[mem] = INFO(cpyst, chrom)
        synolog_grps[grpID][spp].append(gene)
    
    fh.close()

    print(f"Loaded {tot} genes from synolog")

    return 0

def load_orthodb() -> int:
    """create the reference orthodb orthogroups"""

    global odb_mapping, genes_info, GAR, MEDAKA, ZEBRAFISH

    fh   = open_file(orthodb_file)
    gcnt = 0 # number of gar genes loaded
    zcnt = 0 # number of zebrafish genes loaded
    mcnt = 0 # number of medaka genes loaded

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line[:-1].split('\t')
        assert len(fields) == 4, f"Malformed line encountered in {orthodb_file}:\n{line}" 
        if (fields[3] != "renamed"): # this gene was not in any of the annotations
            continue
        grpID = fields[0]
        spp   = fields[2].lower()
        if (spp != GAR and spp != MEDAKA and spp != ZEBRAFISH):
            continue # we only care about these species
        gene  = fields[1].lower()
        mem   = f"{spp}_{gene}"
        if (spp == GAR):
            odb_mapping[grpID].gar_counts += 1 # be conservative and count regardless
        if (mem not in genes_info):
            continue # we only care about genes that are found in Synolog
        genes_info[mem].orthodb_grp = grpID
        if (spp == GAR):
            gcnt += 1
        elif (spp == ZEBRAFISH):
            zcnt += 1
        else:
            mcnt += 1

    fh.close()

    print(f"Loaded {gcnt} spotted gar genes, {zcnt} zebrafish genes, and {mcnt} medaka genes in {os.path.basename(orthodb_file)}")

    return 0

def map_syn_to_odb() -> int:
    """add synolog genes to their orthodb counterparts"""

    global synolog_grps, genes_info, odb_mapping, ZEBRAFISH, MEDAKA

    fishes = [ZEBRAFISH, MEDAKA]

    for grpID, grp in synolog_grps.items():
        for spp in fishes:
            if (spp not in grp):
                continue
            genes = grp[spp]
            for gene in genes:
                mem = f"{spp}_{gene}"
                odb = genes_info[mem].orthodb_grp
                if (odb == ''):
                    continue # gene not in orthodb, so not useful for this application
                mems = odb_mapping[odb].synolog_grps[grpID]
                if (mem not in mems):
                    mems.append(mem)

def acceptable(spp: str, genes: list[str]) -> bool:
    """helper function to determine if the gene(s) are single-copy or tandem paralogs"""

    global genes_info

    non_tandem = 0
    for gene in genes:
        mem         = f"{spp}_{gene}"
        non_tandem += 1 if (genes_info[mem].copy_state != 'T') else 0
        if (non_tandem > 1):
            return False
    return True # we never got past 1

def on_diff_chrom(mem: str, used: set[str]) -> bool:
    """helper function if this member is on a different chrom"""

    global genes_info

    try:
        chrom = genes_info[mem].chrom
        # only chromosomal level sequences
        return (chrom.startswith("NC_") and chrom not in used)
    except:
        return False # not found in Synolog, so not useful

def all_on_chroms(spp: str, genes: list[str]) -> bool:
    """helper function to see if all genes are on all chroms and not unplaced scaffolds"""

    global genes_info

    for gene in genes:
        mem = f"{spp}_{gene}"
        if (genes_info[mem].chrom.startswith("NC_") == False):
            return False

    return True

def all_on_diff_chroms(mems: list[str], used: dict[str, set[str]]) -> bool:
    """helper function to see if all paralogs are on chroms the gar group does not use"""

    for mem in mems:
        idx = mem.find('_') # find first _ char
        if (on_diff_chrom(mem, used[mem[:idx]]) == False):
            return False
    return True

def assign_candidate_category(grp: dict[str, list[str]]) -> str:
    """assign a candidate a category"""

    global GAR, MEDAKA, ZEBRAFISH, categories

    if (GAR in grp):
        return categories[2]
    spp_mems = set(grp.keys())
    if (spp_mems <= {ZEBRAFISH, MEDAKA}):
        return categories[0]
    return categories[1]

def add_candidates(gar_group: GarGroup) -> int:
    """give this gar group a list of candidate orthodb orthogroups"""

    global synolog_grps, genes_info, odb_mapping, ZEBRAFISH, MEDAKA

    # ensure paralogs are on separate chroms
    used = {ZEBRAFISH : set(), MEDAKA : set()}
    odbs = set()

    # collect the chroms and orthodb ids
    for spp, genes in ((ZEBRAFISH, gar_group.zf_genes), (MEDAKA, gar_group.md_genes)):
        for gene in genes:
            info = genes_info[f"{spp}_{gene}"]
            used[spp].add(info.chrom)

            if (info.orthodb_grp != ''):
                odbs.add(info.orthodb_grp)

    # enforce consistency
    for odb in sorted(odbs):
        syn_grps = odb_mapping[odb].synolog_grps
        for other, mems in syn_grps.items():
            if (other == gar_group.id):
                continue # this is the same synolog orthogroup
            grp = synolog_grps[other]
            zf_ = grp.get(ZEBRAFISH, [])
            md_ = grp.get(MEDAKA, [])
            if (len(zf_) == 0 or len(md_) == 0):
                continue # this other orthogroup is missing genes
            tag1 = f"{ZEBRAFISH}_"
            tag2 = f"{MEDAKA}_"

            # collect the members from orthodbs perspective
            zf = list()
            md = list()
            for mem in mems:
                if (mem.startswith(tag1)):
                    zf.append(mem)
                elif (mem.startswith(tag2)):
                    md.append(mem)

            if (len(zf) == 0 or len(md) == 0):
                continue # we are missing genes

            all_mems = zf + md
            if (all_on_diff_chroms(all_mems, used) == False):
                continue # paralogs must be not be linked to gar groups' syntenic region

            cand = Candidate(other, assign_candidate_category(grp), odb_mapping[odb].gar_counts == 1)
            for mem in zf:
                idx = mem.find('_')
                cand.zf_genes.append(mem[idx + 1:])
            for mem in md:
                idx = mem.find('_')
                cand.md_genes.append(mem[idx + 1:])
            gar_group.cands.append(cand)

    return 0

def select_best(cands: list[Candidate]) -> Candidate:
    """return the best candidate based off our ranking scheme"""

    global categories

    best = None
    for cand in cands:
        if (best == None):
            best = cand
            continue
        rank1 = categories.index(cand.category)
        rank2 = categories.index(best.category)
        if (rank1 < rank2 or (rank1 == rank2 and int(cand.id) < int(best.id))):
            best = cand # take lower rank or first in orthologs.tsv

    return best

def find_gar_groups() -> int:
    """find orthogroups that have 1 gar and >= 1 genes for other fish"""

    global synolog_grps, genes_info, counts, gar_groups, GAR, ZEBRAFISH, MEDAKA, \
           ONE_GAR_TWO_FISH, EXCLUDED_COPIES, GAR_ANCHORED, OFF_CHROMS

    for grpID, grp in synolog_grps.items():

        # use get to avoid inserts with default dict
        gar = grp.get(GAR, [])
        zf  = grp.get(ZEBRAFISH, [])
        md  = grp.get(MEDAKA, [])

        if (len(gar) != 1 or len(zf) == 0 or len(md) == 0):
            continue # not useful for finding shared ohnologs

        # count this category
        counts[ONE_GAR_TWO_FISH] += 1

        if (acceptable(ZEBRAFISH, zf) == False or acceptable(MEDAKA, md) == False):
            counts[EXCLUDED_COPIES] += 1
            continue

        counts[GAR_ANCHORED] += 1 # this is acceptable

        for spp, genes in ((ZEBRAFISH, zf), (MEDAKA, md)):
            if (all_on_chroms(spp, genes) == False):
                counts[OFF_CHROMS] += 1
                break # only count once

        gar_groups.append(GarGroup(grpID, gar[0], zf, md))

    return 0

def count_groups() -> int:
    """count the gar-anchored groups that have a paralog group, & write pairs"""

    global gar_groups, counts, categories, WITH_PARALOG_GROUP, PARALOG_NO_GAR, \
           PARALOG_NO_GAR_ODB, PARALOG_TELEOST_ONLY

    hdr = "GarGroup\tGarGene\tCorrespondingGroup\tCorrespondingType\tGarGroup_drer\t" + \
          "GarGroup_olat\tCorr_drer\tCorr_olat\n"
    ofh = open("candidate_ohnolog_groups.tsv", 'w')
    ofh.write(hdr)


    for gar_group in gar_groups:
        if (len(gar_group.cands) == 0):
            continue # no paralog group

        counts[WITH_PARALOG_GROUP] += 1

        no_gar = list()
        for cand in gar_group.cands:
            if (cand.category != categories[2]):
                no_gar.append(cand) # candidate doesn't have gar mem

        if (len(no_gar) > 0):
            counts[PARALOG_NO_GAR] += 1
            best = select_best(no_gar)
            row  = [gar_group.id, gar_group.gar_gene, best.id, best.type,
                   ','.join(gar_group.zf_genes), ','.join(gar_group.md_genes),
                   ','.join(best.zf_genes), ','.join(best.md_genes)]
            ofh.write('\t'.join(row) + '\n')

        # do some finer scale counting

    ofh.close()

def main() -> int:
    """entry point to benchmarking"""

    # get arguments
    get_arguments()

    # load datasets
    load_synolog()
    load_orthodb()
    

    return 0

if __name__ == "__main__":
    main()
