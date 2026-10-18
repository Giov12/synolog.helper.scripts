#!/bin/env python3

import argparse
import os
import gzip
from collections import defaultdict
from itertools import combinations

class Index:
    __slots__ = ("synolog_idx", "orthodb_idx", "synolog_grp", "orig_idx")

    def __init__(self) -> None:
        self.synolog_idx = -1 # default value
        self.orthodb_idx = -1
        self.synolog_grp = -1 # group ID
        self.orig_idx    = -1 # if this is ever changed, this mem belonged to a single-species group

class Discrepancy:
    #
    # class that carries the summary of when
    # two genes in synolog are sepearate
    # orthogroups in synolog, but in
    # a single one in orthodb
    #
    __slots__ = ("main_idx", "purity", "is_subset",
                 "n_members", "n_mapped")
    def __init__(self, n: int) -> None:
        self.main_idx    = -1
        self.purity      = 0.0
        self.is_subset   = False
        self.n_members   = n  # synolog member count
        self.n_mapped    = 0  # orthoDB member count

class ParaGene:
    #
    # evidence collected for a gene that was excluded because
    # can be found in a single-species OrthoDB group
    #
    __slots__ = ("pairs", "exact", "is_loc", "spp_in_main")

    def __init__(self, is_loc: bool, spp_in_main: bool) -> None:
        self.pairs       = 0
        self.exact       = False       # same gene symbol as co-member
        self.is_loc      = is_loc
        self.spp_in_main = spp_in_main # main orthogroup

# singleton objects for this analysis
orthodb_grps   = list() # list of sets of orthogroups
synolog_grps   = list()
orthodb_file   = ''
synolog_file   = ''
tot_tp         = 0
tot_fp         = 0
tot_fn         = 0
misplaced_cnt  = 0 # misplaced synolog genes
checked_cnt    = 0 # how many misplaced genes we checked
pairsSkipODB   = 0 # pairs skipped if not in orthodb
pairsSkipSyn   = 0 # pairs skipped if not in synolog
gene_indices   = defaultdict(Index)
pair_stats     = defaultdict(lambda: {"tp": 0, "fp": 0, "fn": 0})
para_stats     = {"would_tp": 0, "would_fp": 0, "no_ref": 0} # would be if gene wasn't paralogous
para_fp_by_grp = defaultdict(int) # count the number of synolog groups with orthodb paralogs
discrepancies  = list()
discreps_orig  = list()           # using original idxs (when dealing with paralogs)
misplaced      = list()
odb_spp_cache  = defaultdict(set) # orthodb idx -> set of species in that orthogroup
para_genes     = dict()           # member -> Paragene

def get_arguments() -> int:
    """get the arguments"""

    global orthodb_file, synolog_file

    desc  = "Benchmark synolog's orthogroups.tsv file to a renamed orthodb dataset"
    shelp = "orthogroups.tsv file generated from synolog"
    ohelp = "A 4 column file generated with find_ensembl_ids for a orthodb dataset"

    parser = argparse.ArgumentParser(description = desc)
    parser.add_argument("-s", "--synolog", help=shelp, required=True, type=str)
    parser.add_argument("-o", "--orthodb", help=ohelp, required=True, type=str)

    args           = parser.parse_args()
    synolog_file = args.synolog
    orthodb_file = args.orthodb

    assert os.path.isfile(synolog_file), f"Could not locate file: {synolog_file}"
    assert os.path.isfile(orthodb_file), f"Could not locate file: {orthodb_file}"

    return 0

def load_orthodb() -> int:
    """create the reference orthodb orthogroups"""

    global orthodb_grps, orthodb_file, gene_indices, odb_spp_cache

    cur = "-1"
    grp = set()
    fh  = gzip.open(orthodb_file, "rt") if orthodb_file.endswith(".gz") else open(orthodb_file, 'r')
    tot = 0
    cnt = 0
    tmp = dict()

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
        if (cur == "-1"):
            # first group
            cur = grpID

        if (cur != grpID):
            orthodb_grps.append(grp)
            cur = grpID
            grp = set()
        idx = len(orthodb_grps)
        gene_indices[mem].orthodb_idx = idx
        grp.add(mem)
        odb_spp_cache[idx].add(spp)
        tmp[mem] = grpID

    fh.close()

    # add the last group
    orthodb_grps.append(grp)

    para = set()
    for grpID, spp_set in tmp.items():
        if (len(spp_set) == 1):
            para.add(grpID)

    if (len(para) > 0):
        excluded = 0

        print(f"\nExcluding {len(para)} Orthodb orthogroups containing only members from a single species")
        for mem_id, mem in gene_indices.items():
            if (tmp[mem_id] in para):
                mem.orig_idx    = mem.orthodb_idx
                mem.orthodb_idx = -1
                excluded        += 1

        print(f"Exclude {excluded} genes found within single-species OrthoDB groups\n")

    msg = f"Orthodb: Skipped {cnt} orthologs not found in the annotations\n" + \
          f"Orthodb: {tot} orthologs were loaded into {len(orthodb_grps)} orthogroups"

    print(msg)

    return 0

def load_synolog() -> int:
    """create load the synolog orthogroups"""

    global synolog_grps, synolog_file

    # similarly to loading orthodb members, load synolog members
    
    cur = "-1"
    grp = set()
    fh  = gzip.open(synolog_file, "rt") if synolog_file.endswith(".gz") else open(synolog_file, 'r')
    tot = 0
    
    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.split('\t')
        assert len(fields) in [12, 10], f"Malformed line encountered in {synolog_file}:\n{line}" 
        tot  += 1
        grpID = fields[0]
        spp   = fields[4].lower()
        idx   = spp.find('.')
        if (idx != -1): # only keep org ID
            spp = spp[:idx]
        gene  = fields[5].lower()
        mem   = f"{spp}_{gene}"
        if (cur == "-1"):
            # first group
            cur = grpID
        if (cur != grpID):
            synolog_grps.append(grp)
            cur = grpID
            grp = set()
        gene_indices[mem].synolog_idx = len(synolog_grps)
        gene_indices[mem].synolog_grp = int(grpID)
        grp.add(mem)
    
    fh.close()
    
    # add the last group
    synolog_grps.append(grp)
    msg = f"Synolog: {tot} orthologs were loaded into {len(synolog_grps)} orthogroups"
    
    print(msg)

    return 0

def get_spp_id(mem: str) -> str:
    """get the species ID from a member ID"""
    idx = mem.find('_')
    return mem[:idx]

def is_para(mem: str) -> bool:
    """return whether or not this member was in a single-species OrthoDB orthogroup"""
    global gene_indices
    return gene_indices[mem].orig_idx != -1

def ref_idx(mem: str) -> int:
    """return the original idx that this member belonged too"""
    global gene_indices

    idx = gene_indices[mem]
    return idx.orig_idx if idx.orig_idx != -1 else idx.orthodb_idx


def create_spp_pair_key(sppA: str, sppB: str) -> tuple[str, str]:
    """return a consistent tuple to use for hashing later"""
    return tuple(sorted((sppA, sppB))) # inner tuple is sorted, then converted to tuple

def compute_discrepancies(targets: list[Discrepancy], use_orig: bool = False) -> int:
    """precompute all the possible discrepancies since our FN rate
       will be higher due to subsetting via synteny"""

    global synolog_grps, gene_indices

    # prefill

    for grp in synolog_grps:
        discrp       = Discrepancy(len(grp))
        n_odb_mems   = 0
        odb_idx_cnts = defaultdict(int)

        # count the frequency of each orthogroup from orthoDB
        for mem in sorted(grp):
            odb_idx = ref_idx(mem) if use_orig else gene_indices[mem].orthodb_idx
            if (odb_idx == -1):
                continue
            n_odb_mems += 1
            odb_idx_cnts[odb_idx] += 1

        # edge case
        # synolog orthogroup not at all found in
        # odb
        
        if (n_odb_mems == 0):
            targets.append(discrp)
            continue

        # get orthoDB idx that shows up
        # most frequently
        cur  = -1
        best = 0
        for idx, cnt in odb_idx_cnts.items():
            if (cnt > best):
                best = cnt
                cur  = idx

        discrp.main_idx  = cur
        discrp.n_mapped  = n_odb_mems
        discrp.purity    = best / n_odb_mems
        discrp.is_subset = discrp.purity > 0.9 # we can change this later
        targets.append(discrp)

    return 0

def check_para_pair(syn_idx: int, memA: str, memB: str, ofh) -> int:
    """check a synolog paired gene doublet where one of the members is in a single-species orthodb group"""

    global para_stats, discreps_orig, para_fp_by_grp, gene_indices, para_genes

    #
    # grab original orthodb indices
    #
    idxA = ref_idx(memA)
    idxB = ref_idx(memB)

    if (idxA == -1 or -idxB == -1):
        para_stats["no_ref"] += 1 # was not within orthodb
        return 0

    discrp = discreps_orig[syn_idx]
    synGrp = gene_indices[memA].synolog_grp

    if ((discrp.is_subset == True and discrp.main_idx in (idxA, idxB)) == False):
        para_stats["would_fp"] += 1
        para_fp_by_grp[synGrp] += 1
        return 0

    #
    # this would have been a true positive, so let's record
    # for later justification
    #
    para_stats["would_tp"] += 1
    ofh.write(f"{synGrp}\t{memA}\t{memB}\n")

    # get gene symbols
    symA    = memA.split('_')[1]
    symB    = memB.split('_')[1]
    entries = ((memA, symA, symB), (memA, symA, symB))

    for entry in entries:
        mem       = entry[0]
        sym       = entry[1]
        other_sym = entry[2]

        if (is_para(mem) == False):
            continue

        if (mem not in para_genes):
            exct    = sym.startswith("loc")
            spp     = get_spp_id(mem)
            in_main = spp in odb_spp_cache(discrp.main_idx)
            para_genes[mem] = ParaGene(exct, in_main)

        rec        = para_genes[mem]
        rec.pairs += 1
        rec.exact  = rec.exact or (sym == other_sym)

    return 0


def calc_recall() -> int:
    """Caclulate Synolog's Recall"""

    # 
    # recall: amongst all the gene pairs within OrthoDB
    # what fraction does synolog also recover?
    # pair -> within the same orthogroup
    #

    global tot_fn, pairsSkipSyn, gene_indices, \
           orthodb_grps, pair_stats

    for grp in orthodb_grps:
            members = sorted(grp) # ensure deterministism
            for memA, memB in combinations(members, 2):
                sppA = get_spp_id(memA)
                sppB = get_spp_id(memB)
    
                if (sppA == sppB):
                    # we can skip paralogs in this case
                    continue

                idxA = gene_indices[memA].synolog_idx
                idxB = gene_indices[memB].synolog_idx
                key  = create_spp_pair_key(sppA, sppB)

                if (idxA == -1 or idxB == -1):
                    # one of them was not in synolog
                    pairsSkipSyn += 1
                    tot_fn       += 1
                    pair_stats[key]["fn"] += 1
                    continue    
                elif (idxA != idxB):
                    tot_fn += 1
                    pair_stats[key]["fn"] += 1 # false negatives, synolog placed them separately
    
    return 0

def calc_prec() -> int:
    """Calculate Synolog's Precision"""

    # 
    # precision: amongst all the gene pairs within synolog
    # what fraction does orthodb have "established"
    # pair -> within the same orthogroup
    #

    global tot_tp, tot_fp, pairsSkipODB, gene_indices, \
           synolog_grps, pair_stats, discrepancies

    ofh = open("false_positives.tsv", 'w')
    pfh = open("para_pairings.tsv", 'w')
    pfh.write("#SynologGroup\tMemA\tMemB\n")

    for syn_idx, grp in enumerate(synolog_grps):
        members = sorted(grp) # ensure deterministism
        for memA, memB in combinations(members, 2):
            sppA = get_spp_id(memA)
            sppB = get_spp_id(memB)

            #
            # we do not need to check tandem duplicates directly
            # we can do see if the paralogs are grouped
            # with the other orthologs (i.e., other species' genes)
            # in the orthogroup
            #
            if (sppA == sppB):   
                continue

            idxA = gene_indices[memA].orthodb_idx
            idxB = gene_indices[memB].orthodb_idx

            if (idxA == -1 or idxB == -1):
                # one of them is not in orthodb or in a paralogous group
                pairsSkipODB += 1

                if (is_para(memA) or is_para(memB)):
                    check_para_pair(syn_idx, memA, memB, pfh)
                continue

            key = create_spp_pair_key(sppA, sppB)

            if (idxA == idxB):
                tot_tp += 1
                pair_stats[key]["tp"] += 1 # true positive
            else:
                # true false positive or are they subset groups?
                discp = discrepancies[syn_idx]
                if (discp.is_subset and discp.main_idx in [idxA, idxB]):
                    tot_tp += 1
                    pair_stats[key]["tp"] += 1 # true positive
                else:
                    tot_fp += 1
                    pair_stats[key]["fp"] += 1 # false positive
                    grp = gene_indices[memA].synolog_grp
                    ofh.write(f"{grp}\t{memA}\t{memB}\n")

    ofh.close()
    pfh.close()

    return 0

def calc_misplaced() -> int:
    """count every synolog gene that orthodb does not believe fits its orthogroup"""

    #
    # the premise behind this function is that if a single
    # synolog gene is within an orthoDB orthogroup, the
    # pairwise comparisons will flag it as N-1 FP hits
    # which skews the statistics
    #

    global misplaced_cnt, checked_cnt, misplaced, \
           gene_indices, discrepancies, synolog_grps

    for syn_idx, grp in enumerate(synolog_grps):
        discp = discrepancies[syn_idx]

        if (discp.main_idx == -1):
            continue # no orthodb group was linked

        for mem in sorted(grp):
            odb_idx = gene_indices[mem].orthodb_idx
            if (odb_idx == -1):
                continue # not in orthodb
            checked_cnt += 1
            if (odb_idx != discp.main_idx):
                misplaced_cnt += 1
                misplaced.append((mem, syn_idx, odb_idx, discp.main_idx))
    return 0

def summarize() -> int:
    """print the results of the comparisons"""

    global tot_tp, tot_fp, tot_fn, pair_stats, \
           pairsSkipODB, pairsSkipSyn, misplaced_cnt, \
           checked_cnt

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

def summarize_para() -> int:
    """helper function to write out the overall paralogous stats"""

    global para_genes, para_stats

    if (len(para_genes) == 0):
        return 1 # nothing to do hear

    s = para_stats
    n = len(para_genes)

    print("\nPairs skipped because of the single-species exclusion (other gene in reference):")
    print(f"  would have been rescued TP: {s['would_tp']}")
    print(f"  would have been FP:         {s['would_fp']}")
    print(f"  sum:                        {s['would_tp'] + s['would_fp']}")
    print(f"  other gene not in reference (already skipped): {s['no_ref']}")

    if (n == 0):
        return 0
    print(f"\nExcluded genes in rescued pairs: {n}")
    print(f"  same symbol as a partner:           {sum(r['exact'] for r in para_genes.values()) / n:.3f}")
    print(f"  its species already in the main OG: {sum(r['spp_in_main'] for r in para_genes.values()) / n:.3f}")
    print(f"  LOC placeholder symbol:             {sum(r['loc'] for r in para_genes.values()) / n:.3f}")

    by_spp = defaultdict(int)
    for mem in para_genes:
        by_spp[get_spp_id(mem)] += 1
    print("  by species:", dict(sorted(by_spp.items())))

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

    return 0

if __name__ == "__main__":
    main()
