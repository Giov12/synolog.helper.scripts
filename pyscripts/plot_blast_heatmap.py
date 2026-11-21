#!/usr/bin/env python3

import argparse
import sys
import os
import gzip
import threading
from   glob import glob
from   collections import defaultdict
from   concurrent.futures import ThreadPoolExecutor, as_completed


# 
# this code is specifically designed for working with the
# OrthoDB test case data in Synolog's manuscript
#

sppCache   = ''
orthoFile  = ''
orgs       = list()
groupIDs   = set()
numThreads = 1

class GENE:

    __slots__ = ("id", "grp", "transcripts", "records")
    
    def __init__(self, id_: str, orthogroup) -> None:
        self.id          = id_
        self.grp         = orthogroup
        self.transcripts = set()
        self.records     = list()

    def process(self, transcript: str) -> None:
        self.transcripts.discard(transcript)

    def done(self) -> bool:
        return len(self.transcripts) == 0
    
    def has_records(self) -> bool:
        return len(self.records) > 0
    
    def submit(self, org: str, incr: int) -> None:
        # add this gene to the orthogroup
        self.grp.add_mem(org, self, incr)

 
class Org:

    __slots__ = ("id", "genes", "transcripts", "_left", "_pairs")

    def __init__(self, id_ : str) -> None:
        self.id          = id_
        self.genes       = dict() # gene id -> Gene
        self.transcripts = dict() # transcript id -> gene id
        self._left       = set()  # genes needing prot ids

    def add_pair(self, gene_id: str, dup_id: str, blockID: str) -> None:

        if (gene_id not in self.genes):
            self.genes[gene_id] = GENE(gene_id, None)
            self._left.add(gene_id)
        
        if (dup_id not in self.genes):
            self.genes[dup_id] = GENE(dup_id, None)
            self._left.add(dup_id)

        self._pairs.append((gene_id, dup_id, int(blockID)))

    def add_transcript(self, gene_id: str, rec_id: str) -> None:
        self.genes[gene_id].transcripts.add(rec_id)
        self.transcripts[rec_id] = gene_id
        self._left.discard(gene_id)

    def has_remaining(self) -> bool:
        return len(self._left) > 0

    def rem_count(self) -> int:
        return len(self._left)
    
    def get_remaining(self) -> str:
        return ','.join(self._left)
    
    def add_record(self, id_: str, header: str, seq: str, longest: bool) -> None:
        gene_id = self.transcripts[id_]  # id_ is a protein/transcript id
        rec     = Seq(header, seq)
        gene    = self.genes[gene_id]

        # we have processed this record
        gene.process(id_)

        if (longest == False or len(gene.records) == 0):
            gene.records.append(rec)
        else:
            if (len(gene.records[0]) < len(rec)):
                # replace the current id
                gene.records[0] = rec
        
        #
        # we have looked at every sequence id
        # and we can now add it to the orthogroup
        #
        if (gene.grp != None and gene.done()):
            gene.submit(self.id, 1)

    def write_remaining(self) -> None:
        for gene in self.genes.values():
            if (gene.done() == False and gene.has_records()):
                gene.submit(self.id, 0)

    def add_remaining(self) -> None:
        for gene in self.genes.values():
            if (gene.done() == False and gene.has_records()):
                gene.submit(self.id, 0)

    def write_pairs(self) -> None:

        global outPath

        #
        # prev will indicate when to
        # restart the pair num
        #
        prev    = -1
        pairNum = 0

        for pair in self._pairs:
            blockID = pair[2]
            if (blockID != prev):
                pairNum = 0
                prev    = blockID
            out = f"{outPath}/{self.id}-seg_block_{blockID}_pair{pairNum}.fasta"
            ofh = open(out, 'w')
            for i in range(2):
                geneID = pair[i]
                gene   = self.genes[geneID]
                for rec in gene.records:
                    header = rec.header + f" ({geneID})\n"
                    seq    = textwrap.fill(rec.seq, width=60) + '\n'
                    ofh.write(header)
                    ofh.write(seq)
            ofh.close()
            pairNum += 1

#############################################################################################################################
#                                                  Command Line Options                                                     #
#############################################################################################################################

def set_arguments(args) -> int:
    """helper function to check the required parameters shared across operations"""

    global path, sppCache, outPath, groupsFile, singleCopy, \
           allOGroups, Longest, numThreads, lock, Segmental

    #
    # path check
    #
    if (args.path == ''):
        msg = "--path is required"
        sys.exit(msg)
    elif (os.path.isdir(args.path) == False):
        msg = f"Could not locate path {args.path}/"
        sys.exit(msg)
    elif (args.path[-1] == '/'):
        args.path = args.path[:-1]

    #
    # cache check
    #
    if (args.cache == ''):
        msg = "--cache is required"
        sys.exit(msg)
    elif (os.path.isdir(args.cache) == False):
        msg = f"Could not locate path {args.cache}/"
        sys.exit(msg)
    elif (args.cache[-1] == '/'):
        args.cache = args.cache[:-1]
    
    #
    # outdir check
    #
    if (args.out_path == ''):
        args.out_path = '.'
    elif (os.path.isdir(args.out_path) == False):
        msg = f"Could not locate path {args.out_path}/"
        sys.exit(msg)
    elif (args.out_path[-1] == '/'):
        args.out_path = args.out_path[:-1]

    #
    # orthogroups file check
    #
    if (args.orthogroups != ''):
        if (os.path.isfile(args.orthogroups) == False):
            msg = f"Could not locate file {args.orthogroups}/"
            sys.exit(msg)
        if (args.single_copy):
            msg = f"--orthogroups will take precedence over --single_copy"
            print(msg)
            args.single_copy = False
    if (type(args.threads) != int):
        msg = "--threads requires an integer as an argument"
        sys.exit(msg)
    else:
        if (args.threads < 1):
            msg = "--threads must greater than 0"
            sys.exit(msg)
        else:
            args.threads = min(args.threads, os.cpu_count())

    # set the global variables
    path       = args.path
    sppCache   = args.cache
    outPath    = args.out_path
    groupsFile = args.orthogroups
    singleCopy = args.single_copy
    allOGroups = args.all
    Longest    = args.longest
    numThreads = args.threads
    Segmental  = args.segmental

    if (groupsFile == '' and allOGroups == False and singleCopy == False and Segmental == False):
        msg = "Neither --orthogroups, --singleCopy, --segmental, or --all were selected. Defaulting to all.."
        print(msg, file=sys.stderr)
        allOGroups = True

    if (numThreads > 1):
        lock = threading.Lock()

    return 0
    
def parse_command_line() -> int:
    """helper function to get the user's arguments to ensure a proper start"""
    
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

    args = parser.parse_args()

    if (len(vars(args)) == 0):
        parser.print_help()
        sys.exit()
    
    # now set the global variables
    set_arguments(args)

    return 0

#############################################################################################################################
#                                            Helper Functions (Segmentals)                                                  #
#############################################################################################################################

def get_segmental_files(segFiles: list[str]) -> int:
    """retrieve all the segmental tsv files in the provided path"""

    global path
    
    suffix = "*_segmentalduplication_clusters.tsv"
    segFiles.extend(glob(f"{path}/{suffix}"))

    if (len(segFiles) == 0):
        msg = f"Could not locate any *{suffix} files in {path}/"
        sys.exit(msg) 

    return 0

#############################################################################################################################
#                                    Parsing Functions (Orthologs/groups/Segmentals)                                        #
#############################################################################################################################

def read_orthogroup_file(group_ids: set[str]) -> int:
    """read in the group ids found in the orthogroup ids"""

    global groupsFile

    fh      = gzip.open(groupsFile, "rt") if groupsFile.endswith(".gz") else open(groupsFile, 'r')
    lineNum = 0

    for line in fh:
        lineNum += 1
        if (len(line) == 0):
            continue
        line = line.strip()
        if (line.isdigit() == False):
            msg = f"Error at line {lineNum} in {groupsFile}\n" + \
                  f"Expected an integer for group ID. " + \
                  f"Found {line}"
            sys.exit(msg)
        groupID = line
        group_ids.add(groupID)

    fh.close()

    return 0


def get_orthogenes() -> int:
    """helper function to create each Org and populate it with the target GEs"""

    global orgs, groupIDs, path, allOGroups, orthogroups, completed, \
           numThreads, outPath, singleCopy

    orthoFile = f"{path}/orthologs.tsv"
    if (os.path.isfile(orthoFile) == False):
        msg = f"Could not locate orthologs.tsv in {path}/"
        sys.exit(msg)

    fh      = open(orthoFile, 'r')
    lineNum = 0
    NColmns = [10, 12]
    count   = 0
    prev    = "-1"
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
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.split('\t')
        Ncol   = len(fields)

        if (Ncol not in NColmns):
            msg = f"Unexpected number of columns found at line {lineNum} in orthologs.tsv\n" + \
                "Expected 12 columns if tandem duplicates included; else 11\n" + \
                f"Offending line: {line}"
            sys.exit(msg) 
    
        groupID = fields[0]

        if (allOGroups == False and groupID not in groupIDs):
            continue
        org     = fields[4]
        gene_id = fields[5]
        added   = False
        retro   = False
        count  += 1

        if (groupID != prev):
            if (ogroup != None):
                if (numThreads == 1):
                    ogroup.outpath = outPath # will be used in single-threaded mode
                orthogroups.append(ogroup)
            # store a reference to this instance
            ogroup = OrthoGroup(groupID, completed)
            prev   = groupID

        for i in range(len(orgs)):
            if (orgs[i].id == org):
                if (singleCopy and Ncol == 12 and fields[7] == 'R'):
                    retro = True
                    break # exclude retro copies
                orgs[i].add_gene(gene_id, ogroup)
                added = True
                break
        if (retro):
            continue
        if (added == False):
            orgs.append(Org(org))
            orgs[-1].add_gene(gene_id, ogroup)

        #
        # internal counter to let know when all genes
        # have been processed
        #
        ogroup.increment()

    fh.close()

    if (ogroup != None):
        if (numThreads == 1):
            ogroup.outpath = outPath
        orthogroups.append(ogroup)

    if (count == 0):
        msg = "Unable to find any orthogroups in orthologs.tsv matching target orthogroups"
        sys.exit(msg)

    global totalCount
    totalCount = len(orthogroups)
    print(f"Total orthogroups to retrieve: {totalCount}", file=sys.stderr)
    print(f"Total number of genetic elements across orthogroups: {count}", file=sys.stderr)
    return 0


#############################################################################################################################
#                                         Parsing Functions (Annotations)                                                   #
#############################################################################################################################

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

    global lock

    # get the org's annotation
    ann = get_annotation(org.id)

    # now to loop through and parse the annotation
    fh      = gzip.open(ann, "rt") if ann.endswith(".gz") else open(ann, 'r')
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
                if (lock != None):
                    with lock:
                        print(msg, file=sys.stderr)
                else:
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

    if (org.has_remaining()):
        rem = org.get_remaining()
        msg = f"Unable to find elements IDs for {org.id} in {ann} for " + \
               "the following genes: " + rem
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

#############################################################################################################################
#                                         Parsing Functions (Fastas)                                                       #
#############################################################################################################################

def get_protein_file(org: str) -> str:
    """grab the peptide fasta for this org"""

    global sppCache

    gidir  = f"{sppCache}/genome_annotations"
    target = org.split('.')

    if (os.path.isdir(gidir) == False):
        msg = f"genome_annotations directory not found in {sppCache}/\n" + \
               "Was this removed?"
        sys.exit(msg)

    # construct orgID.AnnotationID
    target = target[0] + '.' + target[1]
    peps   = glob(f"{gidir}/{target}/*.peptide.fa.gz", recursive=True)

    if (len(peps) == 0):
        msg = f"Unable to find an protein fasta file for {org} in {sppCache}/"
        sys.exit(msg)

    pep = ''
    for pep_ in peps:
        pepFile = os.path.basename(pep_)
        if (pepFile.startswith(target)):
            pep = pep_
            break

    if (pep == ''):
        msg = f"Unable to find a protein fasta file matching species id {org} in {sppCache}"
        sys.exit(msg)

    return pep


def get_seq_id(header: str) -> str:
    """return the first characters before any spacing"""

    id_ = ''
    idx = header.find(' ')

    if (idx == -1):
        id_ = header[1:]
    else:
        id_ = header[1:idx]

    return id_

def flush_orthogroups() -> int:
    """flush out all the incomplete orthogroups"""
    
    global orthogroups, curCount, totalCount, completed

    if (curCount == totalCount):
        return 0

    unfinished = list()

    for orthogroup in orthogroups:
        if (orthogroup.has_remaining()):
            orthogroup.flush()
            unfinished.append(orthogroup.id)

    count = len(unfinished)
    if (count > 0):
        print(f"Failed to find all sequences/transcripts in fasta files for members in {count} orthogroups", file=sys.stderr)
        ids_ = ','.join(unfinished)
        print(f"Orthogroups missing sequences/transcripts: {ids_}", file=sys.stderr)

    return 0

def rename_unfinished() -> int:
    """rename the temporary files that are incomplete"""

    global orthogroups, curCount, totalCount

    if (curCount == totalCount):
        return 0
    
    unfinished = list()

    for orthogroup in orthogroups:
        if (orthogroup.has_remaining()):
            orthogroup.rename()
            unfinished.append(orthogroup.id)

    count = len(unfinished)
    if (count > 0):
        print("Failed to find all sequences/transcripts for members in", count, "orthologroups")
        ids_ = ','.join(unfinished)
        print("Orthogroups missing sequences/transcripts:", ids_)

    return 0



#############################################################################################################################
#                                                    MAIN                                                                   #
#############################################################################################################################

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
