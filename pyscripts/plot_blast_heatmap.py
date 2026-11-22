#!/usr/bin/env python3
from   __future__ import annotations

import argparse
import sys
import os
import gzip
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from   glob import glob
from   typing import TextIO
from   itertools import combinations
from   collections import defaultdict
from   matplotlib.patches import Rectangle
from   concurrent.futures import ThreadPoolExecutor, as_completed

matplotlib.use("Agg") # no display needed

# 
# this code is specifically designed for working with the
# OrthoDB test case data in Synolog's manuscript
#

sppCache   = ''
orthoFile  = ''
orgOrder   = list() # order of orgs in heatmap
orgs       = dict()
groupIDs   = set()
numThreads = 1
label_mat  = False

class GENE:

    __slots__ = ("id", "transcripts", "orthologs", "paralogs", "group_id")

    #
    # orthologs -> Other Org -> Gene_ID -> Best Blast Hit (same orthogroup)
    # paralogs  -> Other Org -> Gene_ID -> Best Blast Hit (other orthogroup)
    #
    def __init__(self, id_: str, group_id: str) -> None:
        self.id          = id_
        self.group_id    = group_id
        self.orthologs   = defaultdict(lambda: defaultdict(float))
        self.paralogs    = defaultdict(lambda: defaultdict(float))
        self.transcripts = set()

    def add_similarity(self, org: str, gene: GENE, similarity: float) -> int:
        # helper function to determine which category this falls into
        if (gene.group_id == self.group_id):
            self.orthologs[org][gene.id] = max(self.orthologs[org][gene.id], similarity)
        else:
            self.paralogs[org][gene.id] = max(self.paralogs[org][gene.id], similarity)
        return 0

class Org:

    __slots__ = ("id", "annID", "name", "genes", "transcripts")

    def __init__(self, id_ : str) -> None:
        fields           = id_.split('.')
        self.id          = id_
        self.annID       = fields[0] + '.' + fields[2]
        self.name        = ''
        self.genes       = dict() # gene id -> Gene
        self.transcripts = dict() # transcript id -> gene id


    def add_transcript(self, gene_id: str, rec_id: str) -> None:
        self.genes[gene_id].transcripts.add(rec_id)
        self.transcripts[rec_id] = gene_id

    def missing_transcripts(self) -> bool:
        return len(self.transcripts) == 0

    def genes_in_group(self, group_id: str) -> list[GENE]:
        """helper function to get all genes in this orthogroup"""
        return [gene for gene in self.genes.values() if gene.group_id == group_id]

    
def parse_command_line() -> int:
    """helper function to get the user's arguments to ensure a proper start"""

    global sppCache, orthoFile, groupIDs, numThreads, label_mat
    
    desc   = "Code used to parse the species cache to generate a heatmap of blast hits used in the Synolog Manuscript"
    parser = argparse.ArgumentParser(description=desc)

    #
    # help messages to shorten argument lines
    #
    ihelp = "Comma separated list of orthogroup IDs to plot (eg., 1,4389,3) [required]"
    chelp = "path to species cache containing the species used in synolog analysis [required]"
    ohelp = "orthologs.tsv file containing the orthogroups to plot [required]"
    thelp = "number of threads [default: 1]"
    lhelp = "If set, label each cell with the average percent identity"
    rhelp = "Comma separated species codes giving the axis order (eg., hsap,mmus,gacu) " + \
            "[default: alphabetical by species name]"
    
    parser.add_argument("--cache",     help=chelp, type=str, required=True)
    parser.add_argument("--orthologs", help=ohelp, type=str, required=True)
    parser.add_argument("--ids",       help=ihelp, type=str, required=True)
    parser.add_argument("--order",     help=rhelp, type=str, default='')
    parser.add_argument("--threads",   help=thelp, type=int, default=1)
    parser.add_argument("--label",     help=lhelp, action="store_true")

    args       = parser.parse_args()
    sppCache   = args.cache
    orthoFile  = args.orthologs
    ids_list   = args.ids
    numThreads = args.threads
    label_mat  = args.label

    assert os.path.isdir(sppCache),   f"Could not locate directory {sppCache}"
    assert os.path.isfile(orthoFile), f"Could not locate file {orthoFile}"

    for id_ in ids_list.split(','):
        if (id_.isdigit() == False):
            sys.exit(f"Invalid orthogroup ID provided: {id_}")
        groupIDs.add(id_)

    assert len(groupIDs) > 1, f"Fewer than 2 orthogroup IDs were provided"
    assert numThreads > 0,    f"--threads must be at least 1"

    if (args.order != ''):
        global orgOrder
        for spp in args.order.split(','):
            spp = spp.strip()
            if (spp != ''):
                orgOrder.append(spp)

    return 0

def open_file(infile: str) -> TextIO:
    """helper function to shorten and in-line conditional"""

    return gzip.open(infile, "rt") if infile.endswith(".gz") else open(infile, 'r')

def run_threaded(func, jobs: list[tuple]) -> int:
    """helper function to run functions in parallel"""

    global numThreads

    with ThreadPoolExecutor(max_workers=numThreads) as executer:
        futures = [executer.submit(func, *job) for job in jobs]
        for future in as_completed(futures):
            future.result()

    return 0

def add_species_names() -> int:
    """parse the species map to give each org their species name"""

    global orgs, sppCache

    spmap = f"{sppCache}/species.map"

    assert os.path.isfile(spmap), f"Error: Unable to locate species.map in synolog database"

    fh     = open_file(spmap)
    orgMap = {org.split('.')[0]: org for org in orgs.keys()}

    for line in fh:
        if (len(line) == 0 or line[0] == '#' or line[0] == '\n'):
            continue
        fields = line.split('\t')
        org    = fields[0]
        if (org in orgMap):
            entry   = fields[1]
            entries = entry.split(' ')
            species = entries[0][0] + ". " + entries[1]
            orgs[orgMap[org]].name = species

    fh.close()

    missing = list()
    for org in orgs.values():
        if (org.name == ''):
            missing.append(org.id)

    if (len(missing) > 0):
        msg = f"Unable to find all orgs in species map\n" + \
              f"Missing: {', '.join(missing)}"
        sys.exit(msg)

    return 0

def get_orthogenes() -> int:
    """helper function to create each Org and populate it with the target GEs"""

    global groupIDs, orthoFile, orgs

    fh      = open_file(orthoFile)
    NColmns = (10, 12)
    count   = 0
    found   = set()

    #
    # if tandem duplicates are include vs not include
    #
    #    with tandems         without tandems
    #    0  -> group id       0 -> group id
    #    1  -> type           1 -> type
    #    2  -> ElementCount   2 -> ElementCount
    #    3  -> OrgCount       3 -> OrgCount
    #    4  -> Org            4 -> Org
    #    5  -> ElementId      5 -> ElementId
    #    6  -> Name           6 -> Name
    #    7  -> CopyNumState   7 -> Chr
    #    8  -> ParalogCount   8 -> StartBp
    #    9  -> Chr            9 -> EndBp
    #    10 -> StartBp
    #    11 -> EndBp
    #

    for line in fh:
        if (len(line) == 0 or line[0] == '#' or line[0] == '\n'):
            continue
        idx = line.find('\t')
        if (line[:idx] not in groupIDs):
            continue

        fields = line.split('\t')
        Ncol   = len(fields)

        if (Ncol not in NColmns):
            msg = f"Unexpected number of columns found in {os.path.basename(orthoFile)}\n" + \
                "Expected 12 columns if tandem duplicates included; else 10\n" + \
                f"Offending line: {line}"
            sys.exit(msg) 
        group_id = fields[0]
        org      = fields[4]
        gene_id  = fields[5]
        if (org not in orgs):
            orgs[org] = Org(org)
        orgs[org].genes[gene_id] = GENE(gene_id, group_id)
        count += 1
        found.add(group_id)

    fh.close()

    if (count == 0):
        msg = "Unable to find any orthogroups in orthologs.tsv matching target orthogroups"
        sys.exit(msg)

    for group_id in groupIDs:
        if (group_id not in found):
            msg = f"Orthogroup {group_id} not found in {os.path.basename(orthoFile)}"
            sys.exit(msg)

    return 0

def get_annotation(org: Org) -> str:
    """grab the annotation file for this org"""

    global sppCache

    gidir  = f"{sppCache}/genome_integrations"

    if (os.path.isdir(gidir) == False):
        msg = f"genome_integrations directory not found in {sppCache}/\n" + \
               "Was this removed?"
        sys.exit(msg)

    anns = glob(f"{gidir}/{org.annID}/*.g*.gz", recursive=True)

    if (len(anns) == 0):
        msg = f"Unable to find an annotation file for {org.id} in {sppCache}/"
        sys.exit(msg)

    ann = ''
    for an_ in anns:
        annFile = os.path.basename(an_)
        if (annFile.startswith(org.annID)):
            ann = an_
            break

    if (ann == ''):
        msg = f"Unable to find an anotation file matching species id {org.id} in {sppCache}"
        sys.exit(msg)

    return ann

def get_gene_id(attrb: str) -> str:
    """return the gene_id for this record"""

    fields  = attrb.split(';')
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

def make_attrb_map(attrb: str) -> dict[str, str]:
    """return the attributes as key value pairs"""

    map    = dict()
    attrb  = attrb.strip(' "\n') # remove new line char
    fields = attrb.split(';')

    for field in fields:
        field = field.strip(' "')
        if (field == ''):
            continue
        idx = field.find(' ') # find index of first space
        if (idx == -1 or idx == len(field) - 1):
            idx = field.find('=')
            if (idx == -1 or idx == len(field) - 1):
                continue
        key = field[:idx]
        key = key.strip(' "')
        val = field[idx + 1:]
        val = val.strip(' "')
        map[key] = val

    return map

def get_org_transcript_ids(org: Org) -> int:
    """grab the protein ids for the genes in org"""

    # get the org's annotation
    ann = get_annotation(org)

    # now to loop through and parse the annotation
    fh      = open_file(ann)
    lineNum = 0
    curGene = ''
    entries = ["cds", "exon"]

    for line in fh:
        lineNum += 1
        if (len(line) == 0 or line[0] == '#' or line[0] == '\n'):
            continue

        fields = line.split('\t')

        if (len(fields) < 9):
            continue # safety mechanism

        feat = fields[2].lower()

        # capture records such as ncgenes or pseudogenes
        if (feat.endswith("gene")):
            gene_id = get_gene_id(fields[8])
            if (gene_id == ''):
                msg = f"Failed to find a gene_id/ID at line {lineNum} in {ann} for {org.id}"
                print(msg, file=sys.stderr)
            if (gene_id in org.genes):
                curGene = gene_id
            else:
                curGene = ''
            continue

        elif (curGene == '' or feat not in entries):
            continue

        attrbMap = make_attrb_map(fields[8])

        if ("protein_id" in attrbMap):
            prot_id = attrbMap["protein_id"]
            if ("protein_version" in attrbMap):
                prot_vs = attrbMap['protein_version']
                prot_id = prot_id + '.' + prot_vs
            org.add_transcript(curGene, prot_id)
        elif ("transcript_id" in attrbMap):
            tran_id = attrbMap["transcript_id"]
            if ("transcript_version" in attrbMap):
                tran_vs = attrbMap["transcript_version"]
                tran_id = tran_id + '.' + tran_vs
            org.add_transcript(curGene, tran_id)
        elif ("Parent" in attrbMap):
            tran_id = attrbMap["Parent"]
            org.add_transcript(curGene, tran_id)

    fh.close()

    if (org.missing_transcripts()):
        msg = f"Unable to find transcript IDs the gene_ids for {org.id} in {ann}"
        sys.exit(msg)

    return 0

def get_transcript_ids() -> int:
    """helper function to loop through and get the transcript ids for each gene"""

    global orgs

    print("Parsing annotations..", end=' ', flush=True, file=sys.stderr)
    run_threaded(get_org_transcript_ids, [(org,) for org in orgs.values()])
    print("done", file=sys.stderr)
    
    return 0

def get_blast_files(org: Org, all_orgs: list[Org]) -> list[tuple[str, Org]]:
    """grab the blast files where this org is the query, paired with the other subject orgs"""

    global sppCache

    bstdir = f"{sppCache}/annotation_blasts"

    if (os.path.isdir(bstdir) == False):
        msg = f"annotation_blasts directory not found in {sppCache}/\n" + \
               "Was this removed?"
        sys.exit(msg)

    # construct orgID.AnnotationID
    bfiles = list()
    for other in all_orgs:
        if (org.id == other.id):
            continue # skip self
        bfile = f"{bstdir}/blast_{org.annID}_{other.annID}.tsv.gz"
        if (os.path.isfile(bfile) == False):
            msg = f"Unable to find a blast file for {org.id} and {other.id}"
            sys.exit(msg)
        bfiles.append((bfile, other))

    return bfiles

def add_blast_hits(org: Org, bfiles: list[tuple[str, Org]]) -> int:
    """collect all the blast hits for the pairwise files"""

    for bfile, orgB in bfiles:
        fh = open_file(bfile)

        for line in fh:
            if (len(line) == 0 or line[0] == '#' or line[0] == '\n'):
                continue
            fields = line.split('\t')
            transA = fields[0]
            transB = fields[1]
            if (transA not in org.transcripts or transB not in orgB.transcripts):
                continue
            geneA     = org.transcripts[transA]
            geneB     = orgB.transcripts[transB]
            percnt_id = float(fields[2])
            org.genes[geneA].add_similarity(orgB.id, orgB.genes[geneB], percnt_id)
        
        fh.close()

    return 0

def load_blast_hits() -> int:
    """orechestrate loading the blast hits"""

    global orgs

    #
    # get all the blast files. Each thread will be responsible for
    # writing to its assigned org
    #
    org_list = list(orgs.values())
    org_list.sort(key= lambda org: org.annID)

    jobs = list()
    for org in org_list:
        bfiles = get_blast_files(org, org_list)
        jobs.append((org, bfiles))

    print("Loading blast hits..", end=' ', flush=True, file=sys.stderr)
    run_threaded(add_blast_hits, jobs)
    print("done", file=sys.stderr)

    return 0

def ortholog_similarity(orgA: Org, geneA: GENE, orgB: Org, geneB: GENE) -> float:
    """similarity between two genes of the SAME orthogroup"""

    AtoB = geneA.orthologs[orgB.id][geneB.id]
    BtoA = geneB.orthologs[orgA.id][geneA.id]

    return (AtoB + BtoA) / 2.0

def paralog_similarity(orgA: Org, geneA: GENE, orgB: Org, geneB: GENE) -> float:
    """similarity between two genes of DIFFERENT orthogroups"""

    AtoB = geneA.paralogs[orgB.id][geneB.id]
    BtoA = geneB.paralogs[orgA.id][geneA.id]

    return (AtoB + BtoA) / 2.0

def mean_ortholog_similarity(orgA: Org, orgB: Org, group: str) -> float:
    """average ortholog similarity over all orgA x orgB gene pairs in one orthogroup"""

    genesA = orgA.genes_in_group(group)
    genesB = orgB.genes_in_group(group)

    # one of the orgs has no gene in this orthogroup
    if (len(genesA) == 0 or len(genesB) == 0):
        return np.nan

    total = 0.0
    count = 0
    for geneA in genesA:
        for geneB in genesB:
            total += ortholog_similarity(orgA, geneA, orgB, geneB)
            count += 1

    return total / count

def mean_paralog_similarity(orgA: Org, groupA: str, orgB: Org, groupB: str) -> float:
    """average paralog similarity over all (orgA genes in groupA) x (orgB genes in groupB)"""

    genesA = orgA.genes_in_group(groupA)
    genesB = orgB.genes_in_group(groupB)

    # one of the orgs has no gene in one of the orthogroups
    if (len(genesA) == 0 or len(genesB) == 0):
        return np.nan

    total = 0.0
    count = 0
    for geneA in genesA:
        for geneB in genesB:
            total += paralog_similarity(orgA, geneA, orgB, geneB)
            count += 1

    return total / count

def get_org_order() -> list[str]:
    """return the order the orgs will be layed out in the heatmap"""

    global orgs, orgOrder

    # default to alphabetical by species name
    if (len(orgOrder) == 0):
        order = list(orgs.keys())
        order.sort(key = lambda org: (orgs[org].name, org))
        return order

    #
    # ensure that the ordering matches what is actually present
    #
    orderMap = dict()
    for org_id in orgs.keys():
        orderMap[org_id.split('.')[0]] = org_id

    order = list()
    for org_id in orgOrder:
        if (org_id not in orderMap):
            msg = f"--order contains a species not in the selected orthogroups: {org_id}"
            sys.exit(msg)
        order.append(orderMap[org_id])

    for org_id in orderMap.keys():
        if (org_id not in orgOrder):
            msg = f"--order is missing species: {org_id}"
            sys.exit(msg)
            
    return order

def construct_matrices(ogroup1: str, ogroup2: str, order: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """construct the 2D matrix by organizing genes by their orthogroup ID"""

    global orgs

    # 
    # Matrix description
    # The matrix will be N by N where N = len(orgs)
    # In most cases, only 2 orthogroup IDs will be provided
    # So the left/bottom triangle will have genes from one orthogroup
    # while the right/top triangle will have genes from the other orthogroup
    # Diagonal lines will be black to illustrate that there was no percent identity
    # calculated to self, because that is always going to be 100% and is not
    # actually computed in Synolog's framework. For each square linking two orthologs,
    # we will take the average percent identity between the two
    #
    # To illustrate how this is advantagous, we will use this matrix to build a very similar
    # heatmap next to the first one, but we instead of the average sequence similarity between
    # orthologs, it will be the average sequence similarity between the paralogs
    # orgs will be used to sort the attributes. 
    #

    N    = len(order)
    orth = np.full((N, N), np.nan) # self comparisons are NAN
    para = np.full((N, N), np.nan)

    for i in range(N):
        orgI = orgs[order[i]]
        for j in range(N):
            if (i == j):
                continue
            orgJ = orgs[order[j]]
            if (i > j):
                orth[i][j] = mean_ortholog_similarity(orgI, orgJ, ogroup1)
            else:
                orth[i][j] = mean_ortholog_similarity(orgI, orgJ, ogroup2)
            para[i][j] = mean_paralog_similarity(orgI, ogroup1, orgJ, ogroup2)

    return (orth, para)

def draw_heatmap(ax, mat: np.ndarray, labels: list[str], cmap, vmin: float, vmax: float) -> object:
    """draw one heatmap panel; diagonal is black, missing data (NaN) is light grey"""

    global label_mat

    N  = len(labels)
    im = ax.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest", aspect="equal")

    # black diagonal (self comparisons not computed)
    for i in range(N):
        ax.add_patch(Rectangle((i - 0.5, i - 0.5), 1, 1, facecolor="black", edgecolor="none"))

    # thin white gaps between cells
    ax.set_xticks(np.arange(-0.5, N, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, N, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", length=0)
    ax.tick_params(which="major", length=0)
    
    # species names on both axes
    ax.set_xticks(np.arange(N))
    ax.set_yticks(np.arange(N))
    ax.set_xticklabels(labels, rotation=45, ha="right", rotation_mode="anchor", style="italic")
    ax.set_yticklabels(labels, style="italic")

    for spine in ax.spines.values():
        spine.set_visible(False)
    
    if (label_mat == False):
        return im
    
    # write the value in each cell; white text on dark cells
    mid = (vmin + vmax) / 2.0
    for i in range(N):
        for j in range(N):
            if (i == j or np.isnan(mat[i][j])):
                continue
            if (mat[i][j] > mid):
                color = "white"
            else:
                color = "#1a1a1a"
            ax.text(j, i, f"{mat[i][j]:.0f}", ha="center", va="center", fontsize=7, color=color)
    
    return im

def plot_pair(group1: str, group2: str, order: list[str], prefix: str) -> int:
    """make the side-by-side ortholog / paralog figure for one pair of orthogroups"""

    global orgs

    orth, para = construct_matrices(group1, group2, order)
    N          = len(order)

    labels = list()
    for org_id in order:
        labels.append(orgs[org_id].name)

    # one shared color scale so the two panels are directly comparable
    vmin = float("inf")
    vmax = -1
    for mat in (orth, para):
        for i in range(N):
            for j in range(N):
                if (np.isnan(mat[i][j]) == False):
                    vmin = min(mat[i][j], vmin)
                    vmax = max(mat[i][j], vmax)

    if (vmax == -1):
        msg = f"No data to plot for orthogroups {group1} and {group2}"
        sys.exit(msg)

    # edge case when they are equal
    if (vmin >= vmax):
        vmax = vmin + 1.0

    cmap = plt.get_cmap("Blues").copy()
    cmap.set_bad("#d9d9d9") # org has no gene in that orthogroup

    side      = max(3.0, 0.42 * N + 1.5)
    fig, axes = plt.subplots(1, 2, figsize=(2 * side + 1.2, side + 0.6), constrained_layout=True)

    draw_heatmap(axes[0], orth, labels, cmap, vmin, vmax)
    im = draw_heatmap(axes[1], para, labels, cmap, vmin, vmax)

    axes[0].set_title(f"Orthologs\nlower: OG {group1}   upper: OG {group2}", fontsize=10)
    axes[1].set_title(f"Paralogs\nrows: OG {group1}   columns: OG {group2}", fontsize=10)

    cbar = fig.colorbar(im, ax=axes, shrink=0.8, pad=0.02)
    cbar.set_label("Mean reciprocal best percent identity (%)")
    cbar.outline.set_visible(False)

    for ext in ("png", "pdf"):
        fig.savefig(f"{prefix}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote {prefix}.png and {prefix}.pdf", file=sys.stderr)

    return 0

def plot() -> int:
    """one figure per pair of orthogroups (just one figure when 2 IDs are given)"""

    global groupIDs

    order  = get_org_order()
    groups = list(groupIDs)
    groups.sort(key=int)
    pairs  = list(combinations(groups, 2))

    for group1, group2 in pairs:
        prefix = f"{group1}_{group2}_percent_identity_heatmap"
        plot_pair(group1, group2, order, prefix)

    return 0

def main() -> int:
    """entry point to this illustration-focused tool"""
    
    # retrieve the arguments
    parse_command_line()

    # determine where the gene_ids are coming from
    get_orthogenes()

    # add species names
    add_species_names()

    # get the target genes' transcripts ids
    get_transcript_ids()

    # load the blast hits
    load_blast_hits()

    # now to plot
    plot()

    return 0
        
if __name__ == '__main__':
    main()
