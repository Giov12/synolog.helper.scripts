#!/bin/env python3

import argparse
import os
import gzip
import glob
import sys
from   collections import defaultdict

genesMap      = dict()
synologGroups = list()
synologOGs    = list() # list of orthogroups from synolog
refOGs        = list() # list of tuples of (str, list) where str == file name, list == gene IDs

class Gene:
    def __init__(self, spp: str, id_: str, chrom: str, start: int, end: int):
        self.spp   = spp
        self.id    = id_
        self.chrom = chrom
        self.idx   = -1
        self.sotho = -1 # synolog orthogroup idx
        self.start = min(start, end)
        self.isRef = False
class OrthoGroup():
    def __init__(self, id_: str) -> None:
        self.id      = id_
        self.members = set()

class Comparison:

    """
    classification: Equal | Superset | Subset | Split | Absent
    missing: list of (gene_id, reason) for RefOG genes not recovered
    split_detail: list of (gene_id, reason) for genes in >1 synolog orthogroup
    extra: set of gene ids present in Synolog's orthogroups but absent in RefOGd (over-merge signal)
    majority_idx : the Synolog orthogroup index containing the most RefOG members (None if nothing was grouped)
    """

    def __init__(self):
        self.classification = ''
        self.missing        = list()
        self.split_reason   = list()
        self.extra_mems     = set()
        self.majority_idx   = None

def get_arguments() -> tuple[str, str, str, int]:
    """get the arguments"""

    d = "Some python code to compare the Synolog orthologs.tsv file to the recoded RefOGs from orthobench"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-r", "--refOGs",  help="path to the directory containing the recoded .txt RefOGs", required=True)
    parser.add_argument("-g", "--gtfs",    help="directory containing gtf files for focal species", required=True)
    parser.add_argument("-s", "--synolog", help="orthologs.tsv file from Synolog",                  required=True)
    parser.add_argument("-d", "--dist",    help="distance used for sliding window",                 default=100)

    args = parser.parse_args()
    rdir = args.refOGs
    gdir = args.gtfs
    osyn = args.synolog
    dist = args.dist

    assert os.path.isdir(rdir), f"Could not locate directory: {rdir}"
    assert os.path.isfile(osyn), f"Could not locate file: {osyn}"
    assert os.path.isdir(gdir),  f"Could not locate directory: {gdir}"
    assert dist > 0, "--dist must be greater than 0"
    
    return (rdir, gdir, osyn, dist)

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
            gffs = glob.glob(f"{gdir}/*.gff")
    
    gtfs.extend(gffs)

    assert len(gtfs) > 0, f"Could not locate any annotation files in {gdir}"

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

    global genesMap

    # get gtf files
    gtfs = get_gtfs(gdir)

    for gtf in gtfs:
        fh     = gzip.open(gtf, "rt") if gtf.endswith(".gz") else open(gtf, 'r')
        chroms = defaultdict(list)
        spp    = os.path.basename(gtf)
        idx    = spp.find('.')
        spp    = spp[:idx]

        for line in fh:
            if (line[0] == '#'): continue
            fields = line.split('\t')
            if (fields[2] != "gene"): continue
            Chr = fields[0]
            atr = fields[8] # attributes
            gID = get_gene_id(atr)
            pSt = int(fields[3])
            pEd = int(fields[4])
            chroms[Chr].append(Gene(spp, gID, Chr, pSt, pEd))
        fh.close()

        for genes in chroms.values():
            genes.sort(key = lambda g: g.start)
            for i, gene in enumerate(genes):
                gene.idx          = i
                genesMap[gene.id] = gene
    
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

    global synologGroups, genesMap

    fh = gzip.open(osyn, "rt") if osyn.endswith(".gz") else open(osyn, 'r')
    OG = None
    ct = 0

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.split('\t')
        grpID  = fields[0] # group id
        gene   = fields[5]
        if (gene not in genesMap):
            msg = f"Error: {gene} not found in any gtf files provided"
            sys.exit(msg)
        if (OG == None):
            OG  = OrthoGroup(grpID)
        if (grpID != OG.id):
            synologGroups.append(OG)
            OG = OrthoGroup(grpID)
        OG.members.add(gene)
        genesMap[gene].sotho = len(synologGroups)
        ct += 1
        
    # add last group
    synologGroups.append(OG)

    fh.close()

    msg = f"{ct} genes across {len(synologGroups)} orthogroups loaded " + \
           "from Synolog results\n"
    print(msg)

    return 0

def load_refOGs(rdir: str) -> int:
    """read in the genes/groups from the other method"""

    global refOGs, genesMap

    refog_files = glob.glob(f"{rdir}/*.txt") if rdir[-1] != '/' else glob.glob(f"{rdir}*.txt")
    if (len(refog_files) == 0):
        msg = f"No *.txt files found in {rdir}"
        sys.exit(msg)

    for refog_file in refog_files:
        fh = open(refog_file, 'r')
        og = list()
        bn = os.path.basename(refog_file) # base name
        for line in fh:
            if (len(line) == 0):
                continue
            mem = line.strip()
            og.append(mem)
            if (mem in genesMap):
                genesMap[mem].isRef = True
        fh.close()
        refOGs.append((bn, og))

    return 0


def compare_to_refOG(memGenes: list[Gene], dist: int) -> Comparison:
    """compare a specific refOG to the orthogroups in synolog"""
    global genesMap, synologGroups

    comparison = Comparison()

    # different classifications for this comparison
    equal    = "Equal"
    superset = "Superset"
    subset   = "Subset"
    split    = "Split"
    absent   = "Absent"

    memSet  = set()
    grouped = list() # list of Gene objs
    missAnn = defaultdict(list)

    for gene in memGenes:
        if (gene.spp == ''): # this is a dummy gene object
            comparison.missing.append((gene.id, "Not in Annotation"))
        elif (gene.sotho == -1):
            missAnn[gene.spp].append(gene)
        else:
            grouped.append(gene)
        memSet.add(gene.id)

    synologOGs  = defaultdict(list)
    synolog_set = set()
    for gene in grouped:
        synologOGs[gene.sotho].append(gene)
        synolog_set.add(gene.id)

    if (synolog_set == memSet and len(missAnn) == 0 and len(comparison.missing) == 0):
        comparison.classification = equal
        comparison.majority_idx   = list(synologOGs.keys())[0]
        return comparison

    # find out which synolog orthogroup has most of the members
    majority = -1
    best     = 0
    for idx, syn_mems in synologOGs.items():
        if (best < len(syn_mems)):
            majority = idx
    majority_members = synologOGs[majority]
    majority_ids     = set()
    for gene in majority_ids:
        set.add(gene.id)

    majority_syn_group = synologGroups[majority].members
    extra = majority_syn_group.difference(memSet)

    split_detail = list()
    if (len(synologOGs) > 1):
        # build species -> Gene lookup within the majority group, for distance checks
        majority_by_species = defaultdict(list)
        for gene in majority_members:
            majority_by_species[gene.spp].append(gene)

        for idx, members in synologOGs.items():
            if idx == majority:
                continue
            for gid, gene in members:
                ref_gene = majority_by_species.get(gene.spp)
                if ref_gene is None:
                    # majority group has no representative of this species at all
                    split_detail.append((gid, "Foreign Group"))
                    continue
                if ref_gene.chrom != gene.chrom:
                    split_detail.append((gid, "Different Chromosome"))
                    continue
                gap = abs(ref_gene.idx - gene.idx)
                if gap <= dist:
                    split_detail.append((gid, "Sequence Similarity"))
                else:
                    split_detail.append((gid, "Distance"))

    # classify
    if (len(synologOGs) > 1):
        classification = split
    elif (len(missing) > 0):
        classification = subset
    elif (len(extra) > 0):
        classification = subset
    else:
        classification = subset

    comparison.classification = classification

    return comparison

def process_refOGs(dist: int) -> int:
    """function to compare refOGs to synolog (unidirectional comparison)"""

    global genesMap, refOGs

    fh = open("refOG.discrepancies.txt", 'w')

    # discrepancies
    discreps   = list()
    totReasons = defaultdict(int)

    for entry in refOGs:
        fname      = entry[0]
        refMems    = entry[1]
        memGenes   = list()
        synologIds = set()
        for mem in refMems:
            gene = genesMap.get(mem, None)
            if (gene == None):
                sidx = -1 # make a dummy gene
                gene = Gene('', mem, '', -1, -1)
            else:
                sidx = gene.sotho
            memGenes.append(gene)
            if (sidx != -1):
                synologIds.add(sidx)
        if (len(synologIds) == 0):
            outline = f"{fname}: Completely missing\n"
            fh.write(outline)
            totReasons["Completely missing"] += 1
            continue
        comparison = compare_to_refOG(memGenes, dist)
        
    fh.close()

    return 0
  
def main() -> int:
    """Compare the orthologs.tsv file to the recoded RefOGs from the Orthobench dataset"""

    # get arguments
    rdir, gdir, osyn, dist = get_arguments()

    # create the gene objects
    make_gene_map(gdir)

    # load the synolog orthogroups
    load_synolog(osyn)

    # load ref OGs
    load_refOGs(rdir)

    # compare to the truth set
    process_refOGs(dist)

    return 0

if __name__ == "__main__":
    main()
