#!/bin/env python3

import argparse
import os
import gzip
import glob
from   collections import defaultdict

"""
    compare the orthologs.tsv file with the reformated file of
    orthofinder or orthomcl while filtering out multi-transcript IDs
"""

synologGenes  = defaultdict(list) # gene -> orthogroups
otherGenes    = defaultdict(list)
synologGroups = list()
otherGroups   = list()
paralogs      = set()
multiGrped    = set()
missingSynGes = set()
splitSynGes   = set()
equalSynGes   = set()
geneLocations = dict()
keptCount     = 0
nSpp          = 0

class OrthoGroup():
    def __init__(self, id_: str) -> None:
        self.id        = id_
        self.spp       = set()
        self.members   = set()
        self.singleSpp = False
        self.hasUniq   = False # the synolog orthogroups
        self.isDup     = False
        self.otherCnt  = 0     # to get some summary stats on
        
class Summary:
    def __init__(self) -> None:
        self.orthoID   = ''     # other group id
        self.status    = ''     # either Missing/Equal/Superset/Subset/MultiGroups/Conflict
        self.othCount  = 0      # num of other members
        self.synCount  = 0      # num of synolog members
        self.groupCnt  = 0      # num of synolog orthogroups
        self.groupIDs  = ''     # ids of synolog groups
        self.singleCpy = False  # full single copy for all orgs
        self.isDup     = False  # if the OrthoGroup is a duplicate
        self.resolved  = False  # true only if no extras/differences
        self.isParalog = False  # other method is single-spp group
        self.MultiGrp  = 0      # num of this orthogroup mems that are multi-grouped
        self.missSyn   = list() # IDs synolog are missing
        self.missReas  = list() # predicted reasons why a gene is missing in Synolog
        self.missOther = list() # IDs that this method didn't find
        self.synExtras = list() # IDs not in this specific orthogroup
        self.othExtras = list() # IDs not in Synolog's orthogroups
        self.orthogrp  = None   # reference to the OrthoGroup() object

    def add_missSyn(self, gene: str, reason: str) -> None:
        # record the genes synolog didn't group but are in this orthogroup
        self.missSyn.append(gene)
        self.missReas.append(reason)

    def add_synExtra(self, entry: str) -> None:
        # notes why this gene is not in the other method's group
        self.synExtras.append(entry)
    
    def add_missOther(self, gene) -> None:
        # genes not grouped by other method
        self.missOther.append(gene)

    def add_othExtra(self, entry) -> None:
        # notes on why this gene is in a different orthogroup
        self.othExtras.append(entry)

    def has_missing(self) -> bool:
        return len(self.missSyn) > 0
    
    def get_missing(self) -> str:
        # construst the lines for missing genes & their reasons
        outlines = [''] * len(self.missSyn)
        para     = "Paralogous - " if (self.isParalog) else ''

        for i in range(len(outlines)):
            outlines[i] = f"{para}{self.orthoID}: {self.missSyn[i]} Missing - {self.missReas[i]}"

        return '\n'.join(outlines) + '\n'

    def outline(self) -> str:
        # construct the entry to write to the output

        if (self.resolved):
            resolved = 'T'
            if (self.groupCnt > 1):
                self.status = "Split"
        else:
            resolved = 'F'

        if (self.isParalog):
            paralog = "Paralogous"
        else:
            paralog = "Orthologous"

        if (self.isDup):
            paralog = paralog + "-Duplicate"

        synMissCnt = len(self.missSyn)
        othMissCnt = len(self.missOther)
        othMissGes = ','.join(self.missOther)
        synMissGes = ','.join(self.missSyn)
        multigrp   = f"{self.MultiGrp} Multi-Grouped"
        synNotes   = ','.join(self.synExtras)
        othNotes   = ','.join(self.othExtras)

        outline = f"{self.orthoID}\t{self.status}\t{paralog}\t" + \
                  f"{multigrp}\t{resolved}\t{self.groupCnt}\t" + \
                  f"{self.groupIDs}\t{self.othCount}\t{self.synCount}\t" + \
                  f"{othMissCnt}\t{othMissGes}\t{synMissCnt}\t" + \
                  f"{synMissGes}\t{othNotes}\t{synNotes}\n"

        return outline

def get_arguments() -> tuple[bool, str, str, str]:
    """get the arguments"""

    d = "Some python code to compare either orthomcl or orthofinder to synolog"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-f", "--filter",  help="filter for orthogroups with 2 or more species",    action="store_true")
    parser.add_argument("-m", "--method",  help="reformated tsv of either orthofinder or orthomcl", required=True)
    parser.add_argument("-g", "--gtfs",    help="directory containing gtf files for focal species", required=True)
    parser.add_argument("-s", "--synolog", help="orthologs.tsv file from Synolog",                  required=True)
    parser.add_argument("-d", "--dist",    help="distance used for sliding window",                 default=25)

    args = parser.parse_args()
    fltr = args.filter
    omtd = args.method
    gdir = args.gtfs
    osyn = args.synolog
    dist = args.dist

    assert os.path.isfile(omtd), f"Could not locate file: {omtd}"
    assert os.path.isfile(osyn), f"Could not locate file: {osyn}"
    assert os.path.isdir(gdir),  f"Could not locate directory: {gdir}"
    assert dist > 0, "--dist must be >=1"
    
    return (fltr, omtd, gdir, osyn, dist)

def get_gtfs(gdir: str) -> list[str]:
    """helper function to get the gtf files in the provided directory"""

    gtfs = list()
    gffs = list()

    if (gdir[-1] == '/'):
        gtfs = glob.glob(f"{gdir}*.gtf.gz")
        if (len(gtfs) == 0):
            gtfs = glob.glob(f"{gdir}*.gtf")
    else:
        gtfs = glob.glob(f"{gdir}/*.gtf.gz")
        if (len(gtfs) == 0):
            gtfs = glob.glob(f"{gdir}/*.gtf")

    if (gdir[-1] == '/'):
        gffs = glob.glob(f"{gdir}*.gff.gz")
        if (len(gtfs) == 0):
            gffs = glob.glob(f"{gdir}*.gff")
    else:
        gffs = glob.glob(f"{gdir}/*.gff.gz")
        if (len(gtfs) == 0):
            gffs = glob.glob(f"{gdir}/*.gtf")
    
    gtfs.extend(gffs)

    assert len(gtfs) > 0, f"Could not locate any gtf files in {gdir}"

    return gtfs

def get_gene_id(column: str) -> str:
    """parse the attributes column of a gtf to get the gene_id or transcript id"""

    fields  = column.split(';')
    gene_id = ''
    id_     = ''

    if (len(fields) == 1):
        if ("gene_id" in fields[0]):
            gene_id = fields[0].replace("gene_id", '')
            gene_id = gene_id.strip(' "\n')
        elif ("ID" in fields[0]):
            gene_id = fields[0].replace("ID", '')
            gene_id = gene_id.strip(' "\n=')
        else:
            gene_id = fields[0].strip(' "\n')
        return gene_id
    
    for field in fields:
        field = field.strip(' "')
        if ((field.startswith("gene_id") == False) and (field.startswith("ID") == False)):
            continue
        subfields = field.strip(' "\n,=').split(' ')
        if (len(subfields) == 1 and '=' in subfields[0]):
            subfields = field.strip(' "\n,=').split('=')
        record_id = subfields[-1]
        record_id = record_id.strip(' "\n')

        # hold onto this id if no gene_id found
        if (field[0] == 'I'):
            id_     = record_id
        else:
            gene_id = record_id
            break
            
    if (gene_id == ''):
        gene_id = id_ # assume an ID= was found

    return gene_id


def make_gene_map(gdir: str) -> int:
    """create a mapping of each gene to its chr & index on that chromosome"""

    global geneLocations, synologGenes, otherGenes, multiGrped, paralogs

    #
    # hold all the genes so we can quantify how many are
    # in the annotations
    #
    annGenes = set()
    sgenes   = set(synologGenes.keys())
    ogenes   = set(otherGenes.keys())

    # get gtf files
    gtfs = get_gtfs(gdir)

    for gtf in gtfs:
        spp = os.path.basename(gtf).split('.')[0]
        fh  = gzip.open(gtf, "rt") if gtf.endswith(".gz") else open(gtf, 'r')
        idx = 0 # index
        cur = None
        tot = 0
        tmp = set()

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
            geneLocations[gID] = (Chr, idx)
            tot += 1
            tmp.add(gID)
            annGenes.add(gID)

        fh.close()

        # calculate 
        synGenes = sgenes.intersection(tmp)
        othGenes = ogenes.intersection(tmp)
        synCount = len(synGenes)
        othCount = len(othGenes)
        synPercn = round(((synCount / tot) * 100), 2)
        othPercn = round(((othCount / tot) * 100), 2)
        paraCnt  = len(paralogs.intersection(tmp))
        mgrpCnt  = len(multiGrped.intersection(tmp))
        parpercn = round(((paraCnt / tot) * 100), 2)
        mgppercn = round(((mgrpCnt / tot) * 100), 2)
        remcnt   = len(othGenes.difference(paralogs))
        rempercn = round(((remcnt / tot) * 100), 2)


        msg = f"\nLoaded {tot} genes from {gtf}\n" + \
              f"{synCount} ({synPercn}%) found in Synolog\n" + \
              f"{othCount} ({othPercn}%) found in Other Method where:\n" + \
              f"\t{paraCnt} ({parpercn}%) are in paralogous groups\n" + \
              f"\t{remcnt} ({rempercn}%) not in paralogous group\n" + \
              f"\t{mgrpCnt} ({mgppercn}%) are multi-grouped"

        print(msg)

    sgenes = set(synologGenes.keys())
    ogenes = set(otherGenes.keys())

    # synolog print
    n = len(sgenes.difference(annGenes))
    p = round(((n / len(synologGenes)) * 100), 2)
    print("\nTotal Number of genes in Synolog not found in annotations:", n, f"({p}%)")

    n = len(annGenes.intersection(sgenes))
    p = round(((n / len(annGenes)) * 100), 2)
    print("Total Number of annotation genes grouped by Synolog:", n, f"({p}%)")

    # synolog does not report paralgous groups & is isoform aware
    print("Total Number of annotation non-multigrouped and non-paralgous genes grouped by Synolog:", n, f"({p}%)")

    # other method print
    n = len(ogenes.difference(annGenes))
    p = round(((n / len(otherGenes)) * 100), 2)
    print("Total Number of genes in Other Method not found in annotations:", n, f"({p}%)")

    n = len(annGenes.intersection(ogenes))
    p = round(((n / len(annGenes)) * 100), 2)
    print("Total Number of annotation genes grouped by Other Method:", n, f"({p}%)")

    n = len(annGenes.intersection(otherGenes).difference(paralogs).difference(multiGrped))
    p = round(((n / len(annGenes)) * 100), 2)
    print("Total Number of annotation non-multigrouped and non-paralgous genes grouped by Other Method:", n, f"({p}%)")

    return 0

def predict_reason_for_missing(missG: str, otherMems: set[str], dist: int) -> str:
    """try to come up with an explanation on why we are missing this gene for this orthogroup"""

    global geneLocations, synologGenes, otherGenes, multiGrped

    reason = ''

    if (missG not in geneLocations):
        reason = "Not in annotations"
    elif (missG in multiGrped):
        reason = "Multi-Grouped"
    else:
        mChr, mIdx = geneLocations[missG] # chr & index of the missing ge
        spp        = missG.split(':')[0]
        sppFound   = False # assume spp not in any group
        sppDifChr  = True # assume different chr
        rdist      = float("inf") # will minimize dist

        # loop through the members for this other group
        for otherMem in otherMems:
            if (otherMem.startswith(spp) == False): 
                continue
            sppFound   = True
            if (otherMem not in geneLocations):
                continue
            sChr, sIDX = geneLocations[otherMem]

            if (sChr == mChr):
                sppDifChr = False
                rdist = min(rdist, abs(sIDX - mIdx))

        if (sppFound == False):
            reason = "Species not found"
        elif (sppDifChr == True):
            reason = "On a different Chr"
        elif (rdist > dist):
            reason = f"Distance exceeds {dist} genes ({mChr} : {rdist} genes away)"
        else:
            reason = "Sequence similarity"

    return reason

def load_synolog(osyn: str) -> int:
    """load the synolog gene into memory"""

    global synologGenes, synologGroups

    fh   = gzip.open(osyn, "rt") if osyn.endswith(".gz") else open(osyn, 'r')
    OG   = None

    for line in fh:
        if (line[0] == '#'):
            continue
        fields = line.split('\t')
        gi     = fields[0] # group id
        spp    = fields[4].split('.')[0]
        ge     = fields[5]
        gene   = f"{spp}:{ge}"
        if (OG == None):
            OG  = OrthoGroup(gi)
        if (gi != OG.id):
            synologGroups.append(OG)
            OG = OrthoGroup(gi)
        OG.members.add(gene)
        OG.spp.add(spp)
        synologGenes[gene].append(OG)
        
    # add last group
    synologGroups.append(OG)

    fh.close()

    msg = f"{len(synologGenes)} genes across {len(synologGroups)} orthogroups loaded " + \
           "from Synolog results\n"
    print(msg)

    return 0

def load_other_method(omtd: str, fltr: bool) -> int:
    """read in the genes/groups from the other method"""

    global otherGenes, otherGroups, paralogs, multiGrped, \
           keptCount, synologGenes, nSpp

    sppList  = list()
    totges   = set()
    duplics  = defaultdict(list)
    dups     = set() 
    fh       = gzip.open(omtd, "rt") if omtd.endswith(".gz") else open(omtd, 'r')
    ofh      = open("paralogs.genes.txt", 'w')
    ofh2     = open("multigrouped.genes.txt", 'w')
    pgrps    = 0 # groups of single species paralogs
    tot      = 0 # total orthogroups
    totrep   = 0 # total reported "orthologs"
    redSing  = 0 # total number of reduced count to single copy
    redMult  = 0 # same as above but excluding all spp
        
    for line in fh:
        fields = line.strip().split('\t')
        if (line[0] == '#'):
            sppList = fields[1:]
            nSpp    = len(sppList)
            continue
        grpId      = fields[0]
        orthogroup = OrthoGroup(grpId)
        seen       = set()
        gecnter    = dict() # org -> ge -> counts
        sppcnt     = 0
        tot       += 1
        for i in range(1, len(fields)):
            if (fields[i] == ''):
                continue

            sppcnt      += 1
            spp          = sppList[i - 1]
            subfield     = fields[i].split(", ")
            gecnter[spp] = defaultdict(int)

            for gene in subfield:
                totrep             += 1
                gecnter[spp][gene] += 1
                if ((gene, i) in seen):
                    continue # multi-transcript
                seen.add((gene, i))
                gene = f"{spp}:{gene}"
                if (gene in otherGenes):
                    multiGrped.add(gene)
                    ofh2.write(gene + '\n')
                orthogroup.members.add(gene)
                orthogroup.spp.add(spp)
                otherGenes[gene].append(orthogroup)
                totges.add(gene)

        # if not filtering or at least 2 spp. present
        kept = False
        if ((fltr == False) or (sppcnt > 1)):
            kept       = True
            keptCount += 1
        else:
            pgrps += 1
            orthogroup.singleSpp = True
            for gene in orthogroup.members:
                paralogs.add(gene)
                ofh.write(gene)
                ofh.write('\n')
        if (kept):
            reduced = False
            for cnter in gecnter.values():
                ngenes = len(cnter)
                if (ngenes == 1):
                    # only one gene found
                    for cnt in cnter.values():
                        reduced = cnt > 1
                else:
                    for cnt in cnter.values():
                        if (cnt > 1):
                            reduced = True
            if (reduced):
                if (len(gecnter) == len(sppList)):
                    redSing  += 1
                else:
                    redMult += 1
        okey = ','.join(sorted(orthogroup.members))
        if (okey in duplics):
            dups.add(okey) # store index
        duplics[okey].append(len(otherGroups))
        otherGroups.append(orthogroup)

    fh.close()
    ofh.close()
    ofh2.close()

    sgenes = set(synologGenes.keys())
    para   = sgenes.intersection(paralogs)
    mgrp   = sgenes.intersection(multiGrped)

    # mark the duplicates
    dupCnt = 0
    for dup in dups:
        idxs = duplics[dup]
        dupCnt += len(idxs)
        for i in idxs:
            otherGroups[i].isDup = True

    print(f"Number of other method \"genes\": {totrep}")
    print(f"Number of other method identified genes: {len(totges)}")
    print(f"Number of other method orthogroups: {tot}")
    print(f"Number of other method duplicated orthogroups: {dupCnt}")

    p = round(((pgrps / tot) * 100), 2)
    print(f"Number of single-species orthogroups filtered: {pgrps} ({p}%)")

    p = round(((len(paralogs) / len(totges)) * 100), 2)
    print(f"Number of single-species orthologs filtered: {len(paralogs)} ({p}%)")

    p = round(((len(para) / len(paralogs)) * 100), 2)
    print(f"Number of single-species orthologs grouped by Synolog: {len(para)} ({p}%)")

    r = round(((redSing / tot) * 100), 2)
    print(f"Number of orthogroups reduced to single-copy full: {redSing} ({r}%)")

    r = round(((redMult / tot) * 100), 2)
    print(f"Number of orthogroups reduced to single-copy partial: {redMult} ({r}%)")

    r = round(((len(multiGrped) / len(totges))) * 100, 2)
    print(f"Number of multigrouped genes: {len(multiGrped)} ({r}%)")

    if (len(multiGrped) > 0):
        r = round(((len(mgrp) / len(multiGrped))) * 100, 2)
    else:
        r = 0.0
    print(f"Number of multigrouped genes grouped by Synolog: {len(mgrp)} ({r}%)\n")

    ogenes = set(otherGenes.keys())
    shared = sgenes.intersection(ogenes)

    r1 = round(((len(shared) / len(synologGenes))) * 100, 2)
    r2 = round(((len(shared) / len(otherGenes))) * 100, 2)
    print(f"Number of overlapping genes between methods: {len(shared)} ({r1}% in Synolog; {r2}% in Other Method)")

    shared = shared.difference(paralogs)
    nonpar = ogenes.difference(paralogs)
    r1     = round(((len(shared) / len(synologGenes))) * 100, 2)
    r2     = round(((len(shared) / len(nonpar))) * 100, 2)
    print(f"Number of overlapping genes between methods excluding paralogous groups: {len(shared)} ({r1}% in Synolog; {r2}% in Other Method)")

    suniq = sgenes.difference(totges)
    sunip = sgenes.difference(totges.difference(paralogs))
    ouniq = totges.difference(sgenes)
    ounip = totges.difference(paralogs).difference(sgenes)
    print(f"Number of Unique Genes in Synolog not found in Other Method: {len(suniq)} ({len(sunip)} excluding paralogs)")
    print(f"Number of Unique Genes in Other Method not found in Synolog: {len(ouniq)} ({len(ounip)} excluding paralogs)")

    return 0

def summarize_other_orthogroup(orthogroup: OrthoGroup, orthoSum: Summary, dist: int) -> int:
    """helper function to summarize the given orthogroup"""

    global synologGenes, synologGroups, otherGenes, missingSynGes, \
           splitSynGes, equalSynGes, multiGrped, nSpp

    orthoSum.orthoID   = orthogroup.id
    otherMems          = orthogroup.members
    orthoSum.othCount  = len(otherMems)
    orthoSum.isParalog = orthogroup.singleSpp
    orthoSum.isDup     = orthogroup.isDup
    orthoSum.singleCpy = (len(otherMems) == nSpp and len(orthogroup.spp) == nSpp)
    missing            = set() # genes that Synolog did not group
    synGrpsIds         = set()
    synGrps            = list()

    for otherMem in otherMems:
        if (otherMem not in synologGenes):
            missing.add(otherMem)
            continue
        if (otherMem in multiGrped):
            orthoSum.MultiGrp += 1
        synGrp = synologGenes[otherMem][0]
        if (synGrp.id not in synGrpsIds):
            synGrpsIds.add(synGrp.id)
            synGrps.append(synGrp)
    
    # continually grab the genes Synolog did not get
    missingSynGes.update(missing)

    #
    # we have collected all the synolog groups
    # now to check the possible cases by
    # exhaustively checking the possible cases
    #
    orthoSum.groupCnt = len(synGrpsIds)
    if (orthoSum.groupCnt == 0):
        orthoSum.status = "Missing"
        return 0

    #
    # will will grab all the members irrespective
    # of group count
    #
    synMems            = set()
    orthoSum.groupIDs  = ",".join(list(synGrpsIds))
    for syngrp in synGrps:
        synMems.update(syngrp.members)
    
    orthoSum.synCount = len(synMems)

    # set the differences
    synDiff  = set()
    othDiff  = set()

    if (len(synMems) == len(otherMems)):
        if (synMems == otherMems):
            orthoSum.resolved = True
            orthoSum.status   = "Equal"
            if (orthoSum.groupCnt == 1):
                equalSynGes.update(synMems)
            else:
                splitSynGes.update(synMems)
        else:
            synDiff  = synMems.difference(otherMems)
            othDiff  = otherMems.difference(synMems)
    elif (len(otherMems) > len(synMems)):
        othDiff = otherMems.difference(synMems)
        if (synMems.issubset(otherMems)):
            orthoSum.status   = "Superset"
        else:
            synDiff = synMems.difference(otherMems)
    else:
        synDiff = synMems.difference(otherMems)
        if (otherMems.issubset(synMems)):
            orthoSum.status   = "Subset"
        else:
            othDiff = otherMems.difference(synMems)

    if (orthoSum.status != "Equal" and len(othDiff) > 0 and len(synDiff) == 0):
        if (othDiff.issubset(multiGrped)):
            orthoSum.status = "MultiGroups"
    
    if (orthoSum.status == ''):
        orthoSum.status = "Conflict"

    #
    # make a best guess on why there are differences
    # in gene assignments
    #

    if (len(synDiff) == 0 and len(othDiff) == 0):
        # no differences are left to investigate
        return 0
    if (len(synDiff) > 0):
        #
        # let's figure out why synolog has members outside of this method's
        # orthogroup
        #
        for gene in synDiff:
            if (gene not in otherGenes):
                orthoSum.add_missOther(gene)
                continue
            # cases
            # this gene is a multi-mapper in different groups
            # this gene is in a paralog group
            # this gene is in a single different group
            #
            geneOtherGrps = otherGenes[gene]
            geneGrpCount  = len(geneOtherGrps)
            geneGrpIds    = list()
            inParalogous  = False

            for geneGrp in geneOtherGrps:
                if (geneGrp.id != orthogroup.id):
                    geneGrpIds.append(geneGrp.id)
                if (geneGrp.singleSpp):
                    inParalogous = True    
            geneGrpIds   = ','.join(geneGrpIds)
            inParalogous = "In Paralog Group" if inParalogous else ''
            multiGrp     = "MultiGroup" if (geneGrpCount > 1) else ''
            out          = f"({gene} {geneGrpCount} Groups:{multiGrp}:{geneGrpIds}:{inParalogous})"
            orthoSum.add_synExtra(out)
    if (len(othDiff) > 0):
        #
        # do the opposite here and figure out why this method
        # has genes that Synolog does not have
        #
        for gene in othDiff:
            if (gene in missing):
                reason = predict_reason_for_missing(gene, otherMems, dist)
                orthoSum.add_missSyn(gene, reason)
                out    = f"({gene} Missing - {reason})"
                orthoSum.add_othExtra(out)
                continue 
                
            geneSynGrp = synologGenes[gene][0]
            grpID      = geneSynGrp.id

            # no shared members in this group
            otherIds  = set([orthogroup.id])
            otherGrps = otherGenes[gene]
            for otherGrp in otherGrps:
                otherIds.add(otherGrp.id)
            otherIds = ','.join(list(otherIds))
            out      = f"({gene} SynologGroup {grpID}: Found in {otherIds})"
            orthoSum.add_othExtra(out)

    orthoSum.orthogrp = orthogroup

    return 0

def print_counts(equalCnt: int, splitCnt: int, splitSynCnt: int, subsetCnt: int, 
                 supersetCnt: int, conflictCnt: int, confParaCnt: int,
                 missingCnt: int, missParaCnt: int, multiGrpRes: int,
                 singleCnts: int, synSingle: int, paraCnt: int) -> int:
    """print some percentages to the console before the program terminates"""

    global otherGroups, synologGroups, missingSynGes, paralogs, equalSynGes, \
           splitSynGes

    numSyn   = len(synologGroups)
    numOther = len(otherGroups)

    # equal counts
    pSyn = round(((equalCnt / numSyn) * 100), 2)
    pOth = round(((equalCnt / numOther) * 100), 2)
    print(f"\nNumber of equal orthogroups: {equalCnt} ({pSyn}% in Synolog; {pOth}% in Other Method)")

    # split counts
    pSyn = round(((splitSynCnt / numSyn) * 100), 2)
    pOth = round(((splitCnt / numOther) * 100), 2)
    print(f"Number of split orthogroups: {splitSynCnt} into {splitCnt} ({pSyn}% in Synolog; {pOth}% in Other Method)")

    # subset counts
    pSyn = round(((subsetCnt / numSyn) * 100), 2)
    pOth = round(((subsetCnt / numOther) * 100), 2)
    print(f"Number of subset orthogroups: {subsetCnt} ({pSyn}% in Synolog; {pOth}% in Other Method)")

    # superset counts
    pSyn = round(((supersetCnt / numSyn) * 100), 2)
    pOth = round(((supersetCnt / numOther) * 100), 2)
    print(f"Number of superset orthogroups: {supersetCnt} ({pSyn}% in Synolog; {pOth}% in Other Method)")

    # number of groups that would be resolved if not due to multi-grouped genes
    pOth = round(((multiGrpRes / numOther) * 100), 2)
    print(f"Number of conflicted orthogroups due to multigrouped genes: {multiGrpRes} ({pOth}% in Other Method)")

    # conflict counts
    pSyn = round(((conflictCnt / numSyn) * 100), 2)
    pOth = round(((conflictCnt / numOther) * 100), 2)
    poth = round(((confParaCnt / conflictCnt) * 100), 2)
    print(f"Number of conflicted orthogroups: {conflictCnt} ({pSyn}% in Synolog; {pOth}% in Other Method; {confParaCnt} ({poth}%) are paralogous)")

    ccnt = conflictCnt - confParaCnt
    poth = round(((ccnt / numOther) * 100), 2)
    print(f"Number of conflicted orthogroups excluding paralogous groups: {ccnt} ({pOth}% in Other Method)")

    # missing counts
    pOth = round(((missingCnt / numOther) * 100), 2)
    poth = round(((missParaCnt / numOther) * 100), 2)
    print(f"Number of missing orthogroups: {missingCnt} ({pOth}% in Other Method); {missParaCnt} are paralogous ({poth}%)")

    mcnt = missingCnt - missParaCnt
    poth = round(((mcnt / numOther) * 100), 2)
    print(f"Number of missing orthogroups excluding paralogous groups: {mcnt} ({poth}% in Other Method)")
    
    mPar = missingSynGes.intersection(paralogs)
    n    = len(missingSynGes)
    poth = round(((n / len(otherGenes)) * 100), 2)
    pOth = round(((len(mPar) / len(otherGenes)) * 100), 2)
    peql = round(((len(equalSynGes) / len(synologGenes)) * 100), 2)
    pSpl = round(((len(splitSynGes) / len(synologGenes)) * 100), 2)
    print(f"Number of Other Method Specific Genes: {n} ({poth}%); {len(mPar)} ({pOth}%) found in paralogous groups")
    print(f"Number of genes in equal orthogroups: {len(equalSynGes)} ({peql}%)")
    print(f"Number of genes in split orthogroups: {len(splitSynGes)} ({pSpl}%)")

    # now for the single copy orthogroups
    pSyn = round(((synSingle / numSyn) * 100), 2)
    pOth = round(((singleCnts / numOther) * 100), 2)
    poth = round(((singleCnts / (numOther - paraCnt)) * 100), 2)
    print(f"Number of full single-copy orthogroups in Synolog: {synSingle} ({pSyn}%)")
    print(f"Number of full single-copy orthogroups in Other Method: {singleCnts} ({pOth}%; {poth}% excluding paralogous groups)")

    return 0

def process_superset(summary: Summary, fh) -> int:
    """func() to reason the assignments outside the synolog subset"""

    global synologGenes, missingSynGes, geneLocations

    orthogroup = summary.orthogrp
    nsynGrps   = summary.groupCnt # number of synolog groups
    syngenes   = defaultdict(set)
    others     = defaultdict(set)
    othersCnt  = 0
    missingCnt = 0

    # isolate the members outside our group
    for mem in orthogroup.members:
        spp = mem.split(':')[0]
        if (mem in synologGenes):
            syngenes[spp].add(mem)
        else:
            others[spp].add(mem)
            othersCnt += 1
        if (mem in missingSynGes):
            missingCnt += 1

    if (othersCnt == missingCnt):
        miss = "All Missing"
    elif (missingCnt > 0):
        miss = "Partially Missing"
    else:
        miss = "None Missing"

    # create the reason keys
    rkeys = ["Not in Annotation", "Different Chrom", "Sequence Similarity",
             "Distance", "Spp Missing"]
    reasons = {key : list() for key in rkeys}
    # indices for reasons
    NIA = rkeys[0]
    DC  = rkeys[1]
    SS  = rkeys[2]
    D   = rkeys[3]
    SM  = rkeys[4]
        
    # we will do this species by species
    for spp, othGenes in others.items():
        if (spp not in syngenes):
            reasons[SM].extend(list(othGenes))
            continue
        synMems = syngenes[spp]

        # we will bin them by chrom locations
        synLocs = defaultdict(list) 
        othLocs = defaultdict(list) # [(idx, gene)]

        for ogene in othGenes:
            if (ogene not in geneLocations):
                reasons[NIA].append(ogene)
                continue
            chrom, idx = geneLocations[ogene]
            othLocs[chrom].append((idx, ogene))
        for sgene in synMems:
            chrom, idx = geneLocations[sgene]
            synLocs[chrom].append((idx, sgene))

        # now to iterate chrom by chrom
        for chrom, ogenes in othLocs.items():
            if (chrom not in synLocs):
                reasons[DC].extend(othGenes)
                continue
            sgenes = synLocs[chrom]

            #
            # sort them to ensure our two pointer
            # approach will work
            #
            ogenes.sort()
            sgenes.sort()            

            # cases
            synwindow = 25 # default for synolog
            i         = 0 # index for synolog

            for ogene in ogenes:
                idx = ogene[0]
                mem = ogene[1]

                # see if we can find a synolog member close to this mem
                while (i < len(sgenes) and sgenes[i][0] < idx - synwindow):
                    i += 1
                if (i < len(sgenes) and abs(idx - sgenes[i][0]) <= synwindow):                
                    # if close, sequence similarity
                    reasons[SS].append(mem)
                else:
                    # if far, it is a distance
                    reasons[D].append(mem)
    # now to write out the reasons
    outline = list()
    for reason, rlist in reasons.items():
        if (len(rlist) == 0):
            continue
        mems = ','.join(rlist)
        out  = reason + f" {len(rlist)}" + ": " + mems
        outline.append(out)
    
    tmpline = "; ".join(outline) if (len(outline) > 0) else ''
    count   = len(orthogroup.members)
    ortho   = "Paralogous" if summary.isParalog else "Orthologous"
    ortho   = ortho + "-Duplicate" if summary.isDup else ortho
    outline = f"{summary.orthoID} ({count}) {ortho}: "
    outline = outline + f"{summary.groupIDs} ({nsynGrps}); "
    outline = outline + miss + "; "
    outline = outline + "Total: " + f"{othersCnt}; "
    outline = outline + tmpline
    outline = outline + '\n'

    fh.write(outline)
    
    return 0 

def process_subset(summary: Summary, fh) -> int:
    """func() to reason the assignments outside the synolog superset"""

    global synologGenes, otherGenes, missingSynGes, geneLocations, paralogs, \
           multiGrped

    orthogroup = summary.orthogrp
    orthoMems  = orthogroup.members
    synGroups  = list()
    seen       = set()

    # isolate the members outside our group
    for mem in orthogroup.members:
        if (mem in synologGenes):
            synGrp = synologGenes[mem][0]
            if (synGrp.id not in seen):
                seen.add(synGrp.id)
                synGroups.append(synGrp)

    # now collect all the members
    synMembers = set()

    for synGrp in synGroups:
        synMembers.update(synGrp.members)

    # now to pull out the ones not in this orthogroup
    synSpecific = synMembers.difference(orthoMems)

    # the possible cases
    rdict  = dict() # reason dict 

    # these are paralogs
    if (synSpecific.issubset(paralogs)):
        count      = len(synSpecific)
        paras      = ','.join(synSpecific)
        key        = "All paralogs"
        value      = (count, paras)
        rdict[key] = value
    elif (synSpecific.issubset(multiGrped)):
        count      = len(synSpecific)
        multi      = ','.join(synSpecific)
        key        = "All multigrouped"
        value      = (count, multi)
        rdict[key] = value
    # they are missing (not grouped)
    else:
        missing = list()
        for gene in synSpecific:
            if (gene not in otherGenes):
                missing.append(gene)
        if (len(missing) == len(synSpecific)):
            count      = len(synSpecific)
            miss       = ','.join(synSpecific)
            key        = "All ungrouped"
            value      = (count, miss)
            rdict[key] = value

    # case when it's a mixture of the 3
    if (len(rdict) == 0):
        binned    = set()
        # count the number of paralogs
        synPara   = paralogs.intersection(synSpecific)
        count     = len(synPara)
        paraMulti = set()
        if (count > 0):
            paras      = ','.join(synPara)
            key        = "In paraglous groups"
            value      = (count, paras)
            rdict[key] = value
            # check how many of these paralogs are multigrouped
            paraMulti  = synPara.intersection(multiGrped)
            binned.update(synPara)
        #
        # let's separate out any genes that are also
        # in paralogous groups
        #
        if (len(paraMulti) > 0):
            count      = len(paraMulti)
            paraM      = ','.join(paraMulti)
            key        = "Both paralgous and multigrouped"
            value      = (count, paraM)
            rdict[key] = value
            binned.update(paraMulti)
        # count the number of multigrouped
        synMulti = multiGrped.intersection(synSpecific).difference(paraMulti)
        count    = len(synMulti)
        if (count > 0):
            multi      = ','.join(synMulti)
            key        = "Is multigrouped"
            value      = (count, multi)
            rdict[key] = value
            binned.update(synMulti)
        # count the number of missing
        if (len(missing) > 0):
            count      = len(missing)
            miss       = ','.join(missing)
            key        = "Not grouped"
            value      = (count, miss)
            rdict[key] = value
            binned.update(missing)
        # the rest are in separate orthogroups
        if (len(binned) < len(synSpecific)):
            otherGrp   = synSpecific.difference(binned)
            count      = len(otherGrp)
            othMems    = ','.join(otherGrp)
            key        = "Assigned another group"
            value      = (count, othMems)
            rdict[key] = value

    # now to put everything into a single line
    outline = list()

    for reason, values in rdict.items():
        count  = values[0]
        genes  = values[1]
        reason = reason + f" ({count}): "
        reason = reason + genes
        outline.append(reason)
    
    synIDs  = ','.join(seen)
    outline = "; ".join(outline)
    count   = len(seen)
    total   = len(synMembers)
    line    = synIDs + f" (Groups {count} ; Members {total}) " + outline + '\n'
    memCnt  = len(orthogroup.members)
    ortho   = "Paralogous" if summary.isParalog else "Orthologous"
    ortho   = ortho + "-Duplicate" if summary.isDup else ortho
    line    = summary.orthoID + f" ({memCnt}) {ortho}: " + line
    fh.write(line)

    return 0

def process_conflict(summary: Summary, fh) -> int:
    """func() to reason the assignments differentiating the two groups"""

    #
    # since this orthogroup conflicts with 1 or more of our groups
    # this will be a combination of process_subset() and process_superset()
    #

    global synologGenes, otherGenes, missingSynGes, geneLocations, paralogs, \
           multiGrped
    
    orthogroup  = summary.orthogrp
    orthoMems   = orthogroup.members
    othSppdict  = defaultdict(set)
    othSpecific = set()
    synGroups   = list()
    seen        = set()

    reasons = ["All synolog additions are paralogs",
               "All synolog additions are multigrouped",
               "All synolog additions are not grouped",
               "Synolog mix of paralogs and not grouped",
               "Synolog mix of multigroup and not grouped",
               "Synolog mix of paralogs and multigroups",
               "Synolog mix of paralogs, multigroup, and not grouped",
               "All other additions are multigroup",
               "All other additions are missed",
               "Other mix of multigroup and missed",
               ]

    # isolate the members outside our group
    for mem in orthoMems:
        spp = mem.split(':')[0]
        othSppdict[spp].add(mem)
        if (mem in synologGenes):
            synGrp = synologGenes[mem][0]
            if (synGrp.id not in seen):
                seen.add(synGrp.id)
                synGroups.append(synGrp)
        else:
            othSpecific.add(mem)

    # now collect all the members
    synSppdict  = defaultdict(set)
    synMembers  = set()
    synSpecific = set()
    synOnly     = set() # only found in synolog

    for synGrp in synGroups:
        synMembers.update(synGrp.members)
        for mem in synGrp.members:
            spp = mem.split(':')[0]
            synSppdict[spp].add(mem)
            if (mem not in orthoMems):
                synSpecific.add(mem)
                if (mem not in otherGenes):
                    synOnly.add(mem)

    # first, let's handle the genes synolog failed to group
    # create a bucket to store genes
    missReasons = {
        "Not in Annotation":   list(),
        "Different Chrom":     list(),
        "Sequence Similarity": list(),
        "Distance":            list(),
        "Spp Missing":         list()
    }

    missCount = 0

    for mem in othSpecific:
        spp = mem.split(':')[0]
        if (spp not in synSppdict):
            missReasons["Spp Missing"].append(mem)
            missCount += 1
            continue
        synLocs = defaultdict(list)
        othLocs = defaultdict(list)

        for gene in synSppdict[spp]:
            chrom, idx = geneLocations[gene]
            synLocs[chrom].append((idx, gene))
        
        for gene in othSppdict[spp]:
            if (gene not in geneLocations):
                missCount += 1
                missReasons["Not in Annotation"].append(gene)
                continue
            chrom, idx = geneLocations[gene]
            othLocs[chrom].append((idx, gene))

        # now to figure out why we are missing these
        for chrom, ogenes in othLocs.items():
            if (chrom not in synLocs):
                genes      = [ogene[1] for ogene in ogenes]
                missCount += 1
                missReasons["Different Chrom"].extend(genes)
                continue
            sgenes = synLocs[chrom]
            ogenes.sort()
            sgenes.sort()

            synwindow = 25
            i         = 0

            for ogene in ogenes:
                idx = ogene[0]
                mem = ogene[1]

                # see if we can find a synolog member close to this mem
                while (i < len(sgenes) and sgenes[i][0] < idx - synwindow):
                    i += 1
                if (i < len(sgenes) and abs(idx - sgenes[i][0]) <= synwindow):                
                    # if close, sequence similarity
                    missReasons["Sequence Similarity"].append(mem)
                else:
                    # if far, it is a distance
                    missReasons["Distance"].append(mem)
                missCount += 1


    # synolog specific classifications
    synPara   = synSpecific.intersection(paralogs)   # these are in paralog groups
    synMulti  = synSpecific.intersection(multiGrped) # these are multi grouped
    synbinned = synPara.union(synMulti).union(synOnly)
    synUnkwn  = synSpecific.difference(synbinned)    # in theory, this should always be empty 

    # if (len(synUnkwn) > 0 and summary.isParalog == False):
    #     print("I am missing some case with ", orthogroup.id, ','.join(synUnkwn))

    # other specific cases
    othPara   = othSpecific.intersection(paralogs)
    othMulti  = othSpecific.intersection(multiGrped)

    # determine the different cases
    reason  = ''
    parts   = list()
    if (len(synSpecific) == len(synPara)):
        reason = reasons[0]
    elif (len(synSpecific) == len(synMulti)):
        reason = reasons[1]
    elif (len(synSpecific) == len(synOnly)):
        reason = reasons[2]
    if (reason != ''):
        count = len(synSpecific)
        part  = f"{reason} ({count}): "
        part += ','.join(synSpecific)
        parts.append(part)  
    else:
        hasPara  = (len(synPara)  > 0)
        hasMulti = (len(synMulti) > 0)
        hasUniq  = (len(synOnly)  > 0)

        if (hasPara == True and hasMulti == False and hasUniq == True):
            reason = reasons[3]
        elif (hasPara == False and hasMulti == True and hasUniq == True):
            reason = reasons[4]
        elif (hasPara == True and hasMulti == True and hasUniq == False):
            reason = reasons[5]
        elif (hasPara == True and hasMulti == True and hasUniq == True):
            reason = reasons[6]

        if (hasPara and hasMulti):
            ParaMulti = synPara.intersection(synMulti)
            count     = len(ParaMulti)
            if (count > 0):
                part  = f"Marked as Paralog + Multigrouped ({count}): "
                part += ','.join(ParaMulti)
                parts.append(part)             
        else:
            ParaMulti = set()
        
        if (hasPara):
            tmp  = synPara.difference(ParaMulti)
            count = len(tmp)
            if (len(tmp) > 0):
                part  = f"Marked as Paralogs ({count}): "
                part += ','.join(tmp)
                parts.append(part)
        if (hasMulti):
            tmp  = synMulti.difference(ParaMulti)
            count = len(tmp)
            if (len(tmp) > 0):
                part  = f"Marked as Multigrouped ({count}): "
                part += ','.join(tmp)
                parts.append(part)            
        if (hasUniq):
            count = len(synOnly)
            part  = f"Synolog Only ({count}): "
            part += ','.join(synOnly)
            parts.append(part)

    # now process the other orthogroup
    othPara   = othSpecific.intersection(paralogs)
    othMulti  = othSpecific.intersection(multiGrped)
    ParaMulti = othPara.intersection(othMulti)

    # filter if needed
    if (len(ParaMulti) > 0):
        othPara  = othPara.difference(ParaMulti)
        othMulti = othMulti.difference(ParaMulti)

    hasPara  = len(othPara) > 0
    hasMulti = len(othMulti) > 0
    hasBoth  = len(ParaMulti) > 0

    if (hasBoth):
        count = len(ParaMulti)
        part  = f"Other has multigrouped paralogs ({count}): "
        part += ','.join(ParaMulti)
        parts.append(part)
    if (hasPara):
        count = len(othPara)
        part  = f"Other has paralogs ({count}): "
        part += ','.join(othPara)
        parts.append(part)
    if (hasMulti):
        count = len(othMulti)
        part  = f"Other has multigrouped ({count}): "
        part += ','.join(othMulti)
        parts.append(part)

    # now to put everything together into a single line
    count = len(orthoMems)
    ortho = "Paralogous" if (summary.isParalog) else "Orthologous" 
    ortho = ortho + "-Duplicate" if (summary.isDup) else ortho
    line  = f"{summary.orthoID} {ortho} ({count}): "
    line += ','.join(seen) + f" ({len(synMembers)}); "
    line += "; ".join(parts)

    if (missCount > 1):
        parts = list()
        tag   = "Missing reasons:" 
        for reason, genes in missReasons.items():
            if (len(genes) == 0):
                continue
            count = len(genes)
            if (count == missCount):
                tag = "All missing because of"
            part  = f"{reason} ({count}): "
            part += ','.join(genes)
            parts.append(part)
        part  = "; ".join(parts)
        line += f"; {tag} " + part

    line += '\n'
    fh.write(line)

    return 0

def find_duplicated_groups() -> int:
    """write out a file with the duplicated groups"""

    global otherGroups

    # first, bucket the groups by mem count
    buckets = defaultdict(list)

    for orthogroup in otherGroups:
        mems   = orthogroup.members
        sppSet = set()
        for mem in mems:
            spp = mem.split(':')[0]
            sppSet.add(spp)
        sppKey = '#'.join(sorted(sppSet)) + '#' + str(len(mems))
        buckets[sppKey].append(orthogroup)

    # now to construct the file handler
    fh    = open("Duplicated.OrthoGroups.txt", 'w')
    odups = set()
    ogens = set()
    pdups = set()
    pgens = set()

    for orthogroups in buckets.values():
        for i in range(len(orthogroups)):
            ogrp1 = orthogroups[i]
            mems1 = ogrp1.members
            id_1  = ogrp1.id
            cnt   = len(mems1)
            out   = ','.join(mems1)
            ortho = "Paralogous" if ogrp1.singleSpp else "Ortholgous"
            for j in range(i + 1, len(orthogroups)):
                ogrp2 = orthogroups[j]
                mems2 = ogrp2.members
                if (mems1 == mems2):
                    id_2 = ogrp2.id
                    line = f"{id_1} and {id_2}: {ortho} {cnt} members - {out}\n"
                    fh.write(line)
                    if (ortho[0] == 'P'):
                        pdups.add(id_1)
                        pdups.add(id_2)
                        pgens.update(mems1)
                    else:
                        odups.add(id_1)
                        odups.add(id_2)
                        ogens.update(mems1)
    fh.close()

    cnt1 = len(odups)
    cnt2 = len(pdups)
    tot  = cnt1 + cnt2
    cnt3 = len(ogens)
    cnt4 = len(pgens)
    totg = cnt3 + cnt4
    cnt5 = len(ogens.intersection(pgens))
    print(f"Found a total of {tot} duplicated orthogroups containing {totg} genes")
    print(f"{cnt1} duplicated orthogroups are orthologous containing {cnt3} genes")
    print(f"{cnt2} duplicated orthogroups are paralogous containing {cnt4} genes")
    print(f"{cnt5} genes grouped in both duplicated paralgous and duplicated orthologous groups")

    return 0

def compare_orthogroups(dist: int) -> int:
    """iterate over their orthogroups and report for each a comparison summary"""

    global otherGroups, synologGroups, nSpp

    N         = len(otherGroups)
    summaries = list()
    for i in range(N):
        summaries.append(Summary())

    for i in range(N):
        summarize_other_orthogroup(otherGroups[i], summaries[i], dist)

    # keep a couple counters
    equalCnt    = 0
    splitCnt    = 0
    splitSynCnt = 0
    subsetCnt   = 0
    supersetCnt = 0
    conflictCnt = 0
    confParaCnt = 0
    missingCnt  = 0
    missParaCnt = 0
    multiGrpRes = 0
    singleCnts  = 0
    paraCnt     = 0

    # construct the header for the tsv
    header = "OrthoGroupID\tStatus\tOrthology\tNumMultiGrouped\t" + \
             "Resolved\tNumSynologOrthoGroups\tSynologOrthoGroupIDS\t" + \
             "NumGenesInOrthoGroup\tNumSynologGenes\tNumMissingOtherGenes\t" + \
             "MissingOtherGenes\tNumSynologGenesMissing\tSynologMissingGenes\t" + \
             "OtherGeneNotes\tSynologGeneNotes\n"
    
    fh  = open("comparison.tsv", 'w')
    fh2 = open("MissingInSynolog.txt" ,'w')
    fh3 = open("Supersets.txt", 'w')
    fh4 = open("Subsets.txt", 'w')
    fh5 = open("Conflicts.txt", 'w')
    fh6 = open("SplitOrthogroups.txt", 'w')

    fh.write(header)

    for summary in summaries:
        fh.write(summary.outline())
        if (summary.has_missing()):
            fh2.write(summary.get_missing())

        if (summary.status == "Equal"):
            equalCnt += 1
        elif (summary.status == "Split"):
            splitCnt    += 1
            splitSynCnt += summary.groupCnt
            fh6.write(summary.orthoID + '\t' + summary.groupIDs + '\n')
        elif (summary.status == "Subset"):
            subsetCnt += 1
            process_subset(summary, fh4)
        elif (summary.status == "Superset"):
            supersetCnt += 1
            process_superset(summary, fh3)
        elif (summary.status == "Conflict"):
            conflictCnt += 1
            process_conflict(summary, fh5)
            if (summary.isParalog):
                confParaCnt += 1
        elif (summary.status == "MultiGroups"):
            multiGrpRes += 1
        else:
            missingCnt += 1
            if (summary.isParalog):
                missParaCnt += 1
        if (summary.singleCpy):
            singleCnts += 1
        if (summary.isParalog):
            paraCnt += 1

    fh.close()
    fh2.close()
    fh3.close()
    fh4.close()
    fh5.close()
    fh6.close()

    synSingle = 0
    # count the number of single copies for synolog
    for group in synologGroups:
        if (len(group.spp) == nSpp and len(group.members) == nSpp):
            synSingle += 1

    print_counts(equalCnt, splitCnt, splitSynCnt, subsetCnt, 
                 supersetCnt, conflictCnt, confParaCnt,
                 missingCnt, missParaCnt, multiGrpRes,
                 singleCnts, synSingle, paraCnt)

    find_duplicated_groups()

    return 0

def main() -> int:
    """Compare the orthologs.tsv file to either the orthofinder or orthomcl software"""

    # get arguments
    fltr, omtd, gdir, osyn, dist = get_arguments()

    # load the orthologs/groups
    load_synolog(osyn)
    load_other_method(omtd, fltr)

    # load the gene locations
    make_gene_map(gdir)

    # compare the methods
    compare_orthogroups(dist)

    return 0

if __name__ == "__main__":
    main()
