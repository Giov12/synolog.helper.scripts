#!/bin/env python3

import argparse
import os
import gzip
import glob
import sys
from   collections import defaultdict

# global variables
genesMap      = dict()
synologGroups = list()
synologOGs    = list() # list of orthogroups from synolog
refOGs        = list() # list of tuples of (str, list) where str == file name, list == gene IDs
notInSynolog  = set()
synologUniq   = set()

class Gene:
    def __init__(self, spp: str, id_: str, chrom: str, start: int, end: int):
        self.spp   = spp
        self.id    = id_
        self.chrom = chrom
        self.idx   = -1
        self.sotho = -1 # synolog orthogroup idx
        self.sID   = -1 # synolog orthogroup ID
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
    split_reason: list of (gene_id, reason) for genes in >1 synolog orthogroup
    extra_mems: set of gene ids present in Synolog's orthogroups but absent in RefOGs (over-merge signal)
    diff_RefOg: set of gene_ids that are in a seperate RefOG than the current one beign compared
    syngroupsIDs: IDs of synolog orthogroups found represented in the comparison

    """

    def __init__(self):
        self.classification = ''
        self.missing        = list()
        self.split_reason   = list()
        self.extra_mems     = list()
        self.diffRefOG      = list()
        self.syngroupsIDs   = list()

    def get_outputlines(self, fname: str, nMems: int) -> list[str]:
        """generate a summary of the comparison"""

        # fname == file name for the refOG
        # nMems == number of members in fname
        lines  = list()
        synIDs = ", ".join(self.syngroupsIDs)
        header = f"{fname} ({nMems}): {self.classification} (# of Synolog Groups: {len(self.syngroupsIDs)} [{synIDs}])\n"
        lines.append(header)

        missing_genes = defaultdict(list)
        split_genes   = defaultdict(list)

        # bucket the genes by reason
        for entry in self.missing:
            gene_id = entry[0]
            reason  = entry[1]
            missing_genes[reason].append(gene_id)

        for entry in self.split_reason:
            gene_id = entry[0]
            reason  = entry[1]
            split_genes[reason].append(gene_id)

        # is there missing genes
        if (len(missing_genes) > 0):
            parts = list()
            for reason, geneIDs in missing_genes.items():
                mems = ','.join(geneIDs)
                part = f"Missing: {reason} ({len(geneIDs)}): {mems}"
                parts.append(part)
            line = '\t' + "; ".join(parts) + '\n'
            lines.append(line)
        
        # if there are split genes (i.e., multiple synolog orthogroups)
        if (len(split_genes) > 0):
            parts = list()
            for reason, geneIDs in split_genes.items():
                mems = ','.join(geneIDs)
                part = f"Split: {reason} ({len(geneIDs)}): {mems}"
                parts.append(part)
            line = '\t' + "; ".join(parts) + '\n'
            lines.append(line)

        # if there are genes we added that they did not
        if (len(self.extra_mems) > 0):
            mems = ','.join(sorted(self.extra_mems))
            part = f"Extra ({len(self.extra_mems)}): {mems}"
            line = '\t' + part + '\n'
            lines.append(line)

        # if there are genes we added that they did not
        if (len(self.diffRefOG) > 0):
            mems = ','.join(sorted(self.diffRefOG))
            part = f"Misplaced ({len(self.diffRefOG)}): {mems}"
            line = '\t' + part + '\n'
            lines.append(line)
        
        return lines

def get_arguments() -> tuple[str, str, str, int]:
    """get the arguments"""

    d = "Some python code to compare the Synolog orthologs.tsv file to the recoded RefOGs from orthobench"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-r", "--refOGs",  help="path to the directory containing the recoded .txt RefOGs", required=True)
    parser.add_argument("-g", "--gtfs",    help="directory containing gtf files for focal species", required=True)
    parser.add_argument("-s", "--synolog", help="orthologs.tsv file from Synolog",                  required=True)
    parser.add_argument("-d", "--dist",    help="distance used for sliding window", type=int,       default=100)

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
        if (len(gffs) == 0):
            gffs = glob.glob(f"{gdir}*.gff")
    else:
        gffs = glob.glob(f"{gdir}/*.gff.gz")
        if (len(gffs) == 0):
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
            # gene at idx 1 would be 0 positions from index 0 if 
            # starting enumaration at 0
            for i, gene in enumerate(genes, start=1): 
                gene.idx          = i
                genesMap[gene.id] = gene

    print(f"Loaded {len(genesMap)} from {len(gtfs)} annotations")
    
    return 0

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
        genesMap[gene].sID   = grpID
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

    # orthobench data has a duplicated entry: FBpp0309618 / FBgn0003048
    all_genes = set()

    for refog_file in refog_files:
        fh = open(refog_file, 'r')
        og = list()
        bn = os.path.basename(refog_file) # base name
        for line in fh:
            if (len(line) == 0):
                continue
            mem = line.strip()
            og.append(mem)
            all_genes.add(mem)

            if (mem in genesMap):
                genesMap[mem].isRef = True
        fh.close()
        refOGs.append((bn, og))

    print(f"Total number of genes across {len(refOGs)} RefOGs: {len(all_genes)}")

    return 0

def get_reason(gene: Gene, synolog_members: list[Gene], dist: int) -> str:
    """determine why this gene was not grouped into this synolog orthogroup"""

    # reasons
    no_detection = "No Ortholog Detected"
    seq_similar  = "Sequence Similarity"
    diff_chrom   = "Different Chromosome"
    distance     = "Distance"

    # we just didn't find a reason to add an ortholog
    # for this species at all
    if (len(synolog_members) == 0):
        return no_detection

    chrom_found = False
    for mem in synolog_members:
        # shouldn't happen, but just being defensive
        if (mem.id == gene.id):
            continue
        # it could be distance
        if (mem.chrom == gene.chrom):
            chrom_found = True
            gap = abs(mem.idx - gene.idx)
            if (gap <= dist):
                return seq_similar
            
    if (chrom_found == False):
        return diff_chrom

    return distance

def compare_to_refOG(memGenes: list[Gene], dist: int) -> Comparison:
    """compare a specific refOG to the orthogroups in synolog"""

    global genesMap, synologGroups, notInSynolog, synologUniq

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
    missCnt = 0

    for gene in memGenes:
        if (gene.spp == ''): # this is a dummy gene object
            comparison.missing.append((gene.id, "Not in Annotation"))
            missCnt += 1
            notInSynolog.add(gene.id)
        elif (gene.sotho == -1):
            missAnn[gene.spp].append(gene)
            missCnt += 1
            notInSynolog.add(gene.id)
        else:
            grouped.append(gene)
            memSet.add(gene.id)

    if(len(grouped) == 0):
        comparison.classification = absent
        for gene in memGenes:
            comparison.missing.append((gene.id, "No Ortholog Detected"))
        return comparison

    synologOGs  = defaultdict(list)
    synologSpp  = defaultdict(list)
    synolog_set = set()
    synolog_ids = set()
    for gene in grouped:
        synologOGs[gene.sotho].append(gene)
        synologSpp[gene.spp].append(gene)
        synolog_ids.add(genesMap[gene.id].sID)

    # convert to a list
    synolog_ids = sorted(list(synolog_ids))
    comparison.syngroupsIDs.extend(synolog_ids)

    # collect all synolog members
    synolog_set = set()
    for idx in synologOGs.keys():
        synolog_set.update(synologGroups[idx].members)

    # if all of the refOG members are within the collected
    # synolog orthogroups & there are no extra members
    # between the two
    if (synolog_set == memSet and missCnt == 0):
        if (len(synologOGs) == 1):
            comparison.classification = equal
            return comparison
        elif (len(synologOGs) > 1):
            comparison.classification = split
            return comparison

    # now use the current members for this species
    # to see if this is due to synteny (i.e., distance)
    # or sequence similarity
    for spp, ungrouped in missAnn.items():
        grouped_spp = synologSpp.get(spp, list())
        for gene in ungrouped:
            reason = get_reason(gene, grouped_spp, dist)
            comparison.missing.append((gene.id, reason))

    # now to identify why the refOG is split across multiple
    # synolog orthogroups
    spp_by_synOG = defaultdict(lambda : defaultdict(list))
    for gene in grouped:
        spp_by_synOG[gene.spp][gene.sotho].append(gene)

    for spp, sppGrp in spp_by_synOG.items():
        if (len(sppGrp) < 2):
            continue # only one orthogroup with this spp

        # figure out which synolog orthogroup
        # has the most members 
        majority = -1
        best     = 0
        for sotho, mems in sppGrp.items():
            if (best < len(mems)):
                majority = sotho
                best     = len(mems)

        majorityMems = sppGrp[majority]
        for sotho, mems in sppGrp.items():
            if (sotho == majority):
                continue
            for gene in mems:
                reason = get_reason(gene, majorityMems, dist)
                comparison.split_reason.append((gene.id, reason))

    # get the over-merged members
    extra_genes = synolog_set.difference(memSet)

    for gene_id in extra_genes:
        if (genesMap[gene_id].isRef):
            comparison.diffRefOG.append(gene_id)
        else:
            comparison.extra_mems.append(gene_id)
            synologUniq.add(gene_id)

    # note the number of synolog orthogroups
    is_split             = len(synologOGs) > 1
    has_missing          = missCnt > 0
    has_extra            = len(comparison.extra_mems) > 0

    # now add the classification
    split_sub_sup = f"{split}+{subset}+{superset}"
    split_sub     = f"{split}+{subset}"
    split_sup     = f"{split}+{superset}"
    sub_sup       = f"{subset}+{superset}"
    
    if (is_split and has_missing and has_extra):
        # messies of all cases
        comparison.classification = split_sub_sup
    elif (is_split and has_missing):
        # we are split and missing some members
        comparison.classification = split_sub
    elif (is_split and has_extra):
        # we are split and have additional members
        comparison.classification = split_sup
    elif (is_split):
        # we are completely split
        comparison.classification = split
    elif (has_missing and has_extra):
        # single orthogroup with possibly messy assignments
        comparison.classification = sub_sup
    elif (has_missing):
        # we are missing some members from this single group
        comparison.classification = subset
    elif (has_extra):
        # single group with extra members
        comparison.classification = superset
    else:
        # these two orthogroups are equal
        comparison.classification = equal

    return comparison

def process_refOGs(dist: int) -> int:
    """function to compare refOGs to synolog (unidirectional comparison)"""

    global genesMap, refOGs

    fh       = open("refOG.discrepancies.txt", 'w')
    classMap = defaultdict(int)

    for entry in refOGs:
        fname      = entry[0]
        refMems    = entry[1]
        memGenes   = list()
        for mem in refMems:
            gene = genesMap.get(mem, None)
            if (gene == None):
                gene = Gene('', mem, '', -1, -1)
            memGenes.append(gene)
        # do the comparison & construct a summary of the comparison
        comparison = compare_to_refOG(memGenes, dist)
        outlines   = comparison.get_outputlines(fname, len(refMems))
        # write the summary
        for line in outlines:
            fh.write(line)
        classMap[comparison.classification] += 1

    fh.close()

    global synologUniq, notInSynolog
    print(f"Undetected number of genes from RefOGs: {len(notInSynolog)}")
    print(f"Number of introduced genes from Synolog not found in RefOGs: {len(synologUniq)}\n")

    # print the overall findings
    total = len(refOGs)
    for classification, class_count in classMap.items():
        if (class_count > 0):
            percent = round((class_count / total) * 100, 2)
        else:
            percent = 0
        print(f"{class_count} {classification} classifications ({percent}%) across {total} RefOGs")

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
