#!/bin/env python3

import argparse
import os
import gzip
from collections import defaultdict
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
    __slots__ = ("copy_state", "chrom")
    def __init__(self, copy_state: str, chrom: str) -> None:
        self.copy_state = copy_state
        self.chrom      = chrom

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
    __slots__ = ("id", "type", "one_gar", "zf_genes", "md_genes")

    def __init__(self, id_: str, type_: str, one_gar: bool) -> None:
        self.id       = id_
        self.type     = type_   # teleost-only, other species, no gar, or also has gar
        self.one_gar  = one_gar # the shared orthodb orthogroup has exactly 1 gar gene
        self.zf_genes = list()  # zebrafish paralogs in this orthogroup
        self.md_genes = list()  # medaka paralogs in this orthogroup


synolog_file = ''
orthodb_file = ''
genes_info   = dict() # gene -> INFO
orthodb_ogs  = dict() # gene ->
synolog_grps = defaultdict(lambda: defaultdict(list)) # grp -> spp -> [genes]
odb_mapping  = defaultdict(OrthoDB)

# constant globals
GAR = "locu", MEDAKA = "olat", ZEBRAFISH = "drer"

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

    global orthodb_file, odb_mapping

    fh  = open_file(orthodb_file)
    tot = 0
    cnt = 0

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line[:-1].split('\t')
        assert len(fields) == 4, f"Malformed line encountered in {orthodb_file}:\n{line}" 
        if (fields[3] != "renamed"): # this gene was not in any of the annotations
            cnt += 1
            continue
        tot  += 1
        grpID = fields[0]
        spp   = fields[2].lower()
        gene  = fields[1].lower()
        mem   = f"{spp}_{gene}"
        

    fh.close()


    print(msg)

    return 0

def summarize() -> int:
    """print the results of the comparisons"""

    global 

    # calculate precision as TP / (TP + FP)
    pSum = (tot_tp + tot_fp)
    if (pSum == 0):
        precision = 0
    else:
        precision = tot_tp / pSum

    # calculate recall as TP / (TP + FN)
    rSum = (tot_tp + tot_fn)
    if (rSum == 0):
        recall = 0
    else:
        recall = tot_tp / rSum

    # calculate F1 as (2 * Precision * Recall) / (Precision + Recall)
    tot = precision + recall
    if (tot == 0):
        f1 = 0
    else:
        f1 = (2 * precision * recall) / tot

    # print a high level view of the results
    print(f"\nGlobal: TP={tot_tp} FP={tot_fp} FN={tot_fn}")
    print(f"Global: precision={precision:.4f} recall={recall:.4f} f1={f1:.4f}")
    print(f"Skipped (gene not in orthodb reference): {pairsSkipODB}")
    print(f"Skipped (gene not in synolog output):    {pairsSkipSyn}")

    if (checked_cnt == 0):
        gene_precision = 0.0
    else:
        gene_precision = 1 - (misplaced_cnt / checked_cnt)

    print(f"\nGene-level placement: {checked_cnt} genes checked,", 
          f"{misplaced_cnt} misplaced ({gene_precision:.4f} agreement rate)\n")

    # print species view of things
    print(f"{'Species A':<10}{'Species B':<10}{'TP':>8}{'FP':>8}{'FN':>8}{'Precision':>12}{'Recall':>10}")
    keys = sorted(pair_stats.keys())

    for key in keys:
        sppA  = key[0]
        sppB  = key[1]
        stats = pair_stats[key]
        tp    = stats["tp"]
        fp    = stats["fp"]
        fn    = stats["fn"]
        p     = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        r     = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        print(f"{sppA:<10}{sppB:<10}{tp:>8}{fp:>8}{fn:>8}{p:>12.4f}{r:>10.4f}")
    
    return 0


def main() -> int:
    """entry point to benchmarking"""

    # get arguments
    get_arguments()

    # load datasets
    load_orthodb()
    load_synolog()

    # precompute any discrepancies due to synteny
    global discrepancies, discreps_orig
    compute_discrepancies(discrepancies, use_orig=False)
    compute_discrepancies(discreps_orig, use_orig=True)

    # calculate precision & recall
    calc_prec()
    calc_recall()

    # ensure we count penalities only once
    calc_misplaced()

    # summarize
    summarize()
    summarize_para()

    # write some genes for manual inspection
    write_para_sample()

    return 0

if __name__ == "__main__":
    main()
