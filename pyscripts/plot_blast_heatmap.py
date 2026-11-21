#!/usr/bin/env python3

import argparse
import sys
import os
import gzip
from   glob import glob
from   typing import TextIO
from   collections import defaultdict
from   concurrent.futures import ThreadPoolExecutor, as_completed


# 
# this code is specifically designed for working with the
# OrthoDB test case data in Synolog's manuscript
#

sppCache   = ''
orthoFile  = ''
orgs       = dict()
groupIDs   = set()
numThreads = 1

class GENE:

    __slots__ = ("id", "transcripts", "others", "group_id")

    #
    # others -> Other Org -> Gene_ID -> Blast Hit
    #
    def __init__(self, id_: str, group_id: str) -> None:
        self.id          = id_
        self.group_id    = group_id
        self.orthologs   = defaultdict(lambda d: defaultdict(float))
        self.paralogs    = defaultdict(lambda d: defaultdict(float))
        self.transcripts = set()

class Org:

    __slots__ = ("id", "name", "genes", "transcripts")

    def __init__(self, id_ : str) -> None:
        self.id          = id_
        self.name        = ''
        self.genes       = dict() # gene id -> Gene
        self.transcripts = dict() # transcript id -> gene id


    def add_transcript(self, gene_id: str, rec_id: str) -> None:
        self.genes[gene_id].transcripts.add(rec_id)
        self.transcripts[rec_id] = gene_id

    def missing_transcripts(self) -> bool:
        return len(self.transcripts) == 0

    
def parse_command_line() -> int:
    """helper function to get the user's arguments to ensure a proper start"""

    global sppCache, orthoFile, groupIDs, numThreads
    
    desc   = "Code used to parse the species cache to generate a heatmap of blast hits used in the Synolog Manuscript"
    parser = argparse.ArgumentParser(description=desc)

    #
    # help messages to shorten argument lines
    #
    ihelp = "Comma separated list of orthogroup IDs to plot (eg., 1,4389,3) [required]"
    chelp = "path to species cache containing the species used in synolog analysis [required]"
    ohelp = "orthologs.tsv file containing the orthogroups to plot [required]"
    thelp = "number of threads [default: 1]"
    
    parser.add_argument("--cache",       help=chelp, type=str, required=True)
    parser.add_argument("--orthologs",   help=ohelp, type=str, required=True)
    parser.add_argument("--ids",         help=ihelp, type=str, required=True)
    parser.add_argument("--threads",     help=thelp, type=int, default=1)

    args       = parser.parse_args()
    sppCache   = args.cache
    orthoFile  = args.orthologs
    ids_list   = args.ids
    numThreads = args.threads

    assert os.path.isdir(sppCache),   f"Could not locate directory {sppCache}"
    assert os.path.isfile(orthoFile), f"Could not locate file {orthoFile}"

    for id_ in ids_list.split(','):
        if (id_.isdigit() == False):
            sys.exit(f"Invalid orthogroup ID provided: {id_}")
        groupIDs.add(id_)

    assert len(groupIDs) > 0, f"No orthogroup IDs were provided"
    assert numThreads > 0,    f"--threads must be at least 1"

    return 0

def open_file(infile: str) -> TextIO:
    """helper function to shorten and in-line conditional"""

    return gzip.open(infile, "rt") if infile.endswith(".gz") else open(infile, 'r')

def populate(orgs: dict[str, Org]) -> int:
    """helper function to initiate all defaults"""

    for orgA_id, orgA in orgs.items():
        genes = orgA.genes
        for orgB_id, orgB in orgs.items():
            if (orgA_id == orgB_id):
                continue
            for geneA_id, geneA in genes.items():
                for geneB_id, geneB in orgB.genes.items():
                    if (geneA.group_id == geneB.group_id):
                        geneA.orthologs[orgB_id][geneB_id] = 0.0
                    else:
                        geneA.paralogs[orgB_id][geneB_id] = 0.0
    return 0


def get_orthogenes() -> int:
    """helper function to create each Org and populate it with the target GEs"""

    global groupIDs, orthoFile, orgs

    fh      = open_file(orthoFile)
    lineNum = 0
    NColmns = (10, 12)
    count   = 0
    ogroup  = None

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
        lineNum += 1
        if (len(line) == 0 or line[0] == '#' or line[0] == '\n'):
            continue
        idx = line.find('\t')
        if (line[:idx] not in groupIDs):
            continue

        fields = line.split('\t')
        Ncol   = len(fields)

        if (Ncol not in NColmns):
            msg = f"Unexpected number of columns found in {os.path.basename(orthoFile)}\n" + \
                "Expected 12 columns if tandem duplicates included; else 11\n" + \
                f"Offending line: {line}"
            sys.exit(msg) 
        group_id = fields[0]
        org      = fields[4]
        gene_id  = fields[5]
        if (org not in orgs):
            orgs[org] = Org(org)
        orgs[org].genes[gene_id] = GENE(gene_id, group_id)
        count += 1


    fh.close()

    if (count == 0):
        msg = "Unable to find any orthogroups in orthologs.tsv matching target orthogroups"
        sys.exit(msg)

    # create all pairwise mappings
    populate(orgs)

    return 0

def get_annotation(org: str) -> str:
    """grab the annotation file for this org"""

    global sppCache

    gidir  = f"{sppCache}/genome_integrations"
    target = org.split('.')

    if (os.path.isdir(gidir) == False):
        msg = f"genome_integrations directory not found in {sppCache}/\n" + \
               "Was this removed?"
        sys.exit(msg)

    # construct orgID.IntegrationID
    target = target[0] + '.' + target[2]

    anns = glob(f"{gidir}/{target}/*.g*.gz", recursive=True)

    if (len(anns) == 0):
        msg = f"Unable to find an annotation file for {org} in {sppCache}/"
        sys.exit(msg)

    ann = ''
    for an_ in anns:
        annFile = os.path.basename(an_)
        if (annFile.startswith(target)):
            ann = an_
            break

    if (ann == ''):
        msg = f"Unable to find an anotation file matching species id {org} in {sppCache}"
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

def get_org_sequence_ids(org: Org) -> int:
    """grab the protein ids for the genes in org"""

    # get the org's annotation
    ann = get_annotation(org.id)

    # now to loop through and parse the annotation
    fh      = open_file(ann)
    lineNum = 0
    curGene = ''
    entries = ["cds", "exon"]

    for line in fh:
        lineNum += 1
        if (len(line) == '' or line[0] == '#'):
            continue

        fields = line.split('\t')
        feat   = fields[2].lower()

        if (feat == "gene"):
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

    global orgs, numThreads

    print("Parsing annotations..", end=' ', flush=True, file=sys.stderr)

    with ThreadPoolExecutor(max_workers=numThreads) as executor:
        futures = [executor.submit(get_org_sequence_ids, org) for org in orgs]
        for future in as_completed(futures):
            future.result()

    print("done", file=sys.stderr)
    
    return 0

def get_blast_files(idx: int, all_orgs: list[str]) -> list[str]:
    """grab the blast files where this org is the query"""

    global sppCache

    bstdir = f"{sppCache}/annotation_blasts"
    target = all_orgs[idx]

    if (os.path.isdir(bstdir) == False):
        msg = f"annotation_blasts directory not found in {sppCache}/\n" + \
               "Was this removed?"
        sys.exit(msg)

    # construct orgID.AnnotationID
    bfiles = list()
    for i in range(len(all_orgs)):
        if (idx == i):
            continue # skip self hit
        bfile = f"{bstdir}/blast_{target}-{all_orgs[i]}.tsv.gz"
        if (os.path.isfile(bfile) == False):
            msg = f"Unable to find a blast file for {target} and {all_orgs[i]}"
            sys.exit(msg)
        bfiles.append(bfile)

    if (len(bfiles) == 0):
        msg = f"Unable to find any blast files for for {target} in {sppCache}/annotation_blasts"
        sys.exit(msg)

    return bfiles

def add_blast_hits(org: Org, bfiles: list[str], transcript_map: dict[str, Org]):
    """collect all the blast hits for the pairwise files"""

    for bfile in bfiles:
        fname = os.path.basename(bfile)
        other = fname.split('_')[2]
        other = other[:other.find(".tsv")]
        fh    = open_file(bfile)

        for line in fh:
            if (len(line) == 0 or line[0] == '#' or line[0] == '\n'):
                continue
        

        fh.close()


def main() -> int:

    global allOGroups, Segmental
    
    # retrieve the arguments
    parse_command_line()

    # determine where the gene_ids are coming from
    if (Segmental):
        get_segmental_duplications()
    else:
        # retrieve the orthogroup ids if provided
        get_orthogroups()

        # get the target gene ids
        get_orthogenes()

    # get the target genes' transcripts ids
    get_transcript_ids()

    # now to get all the sequences
    get_sequences()

    return 0
        
if __name__ == '__main__':
    main()
