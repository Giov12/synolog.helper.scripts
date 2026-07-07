#!/usr/bin/env python3
#
# Copyright 2025, Gio Madrigal <gm33@illinois.edu> & Julian Catchen <jcatchen@illinois.edu>
#
# This file is part of Synolog.
#
# Synolog is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Synolog is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Synolog.  If not, see <http://www.gnu.org/licenses/>.
#

import argparse
import sys
import os
import gzip
import textwrap
from   glob import glob
from   collections import defaultdict

#############################################################################################################################
#                                                  Global Variables                                                         #
#############################################################################################################################

path       = ''
orgID      = ''
annID      = ''
ref_fasta  = ''
outName    = ''
chroms     = defaultdict(list) # list of Transcript objects

#############################################################################################################################
#                                                      Classes                                                              #
#############################################################################################################################

class Seq:
    
    __slots__ = ("header", "seq")

    def __init__(self, header: str, seq: str) -> None:
        self.header = header
        self.seq    = seq

    def __len__(self) -> int:
        return len(self.seq)

class Transcript:

    __slots__ = ("id", "exons")

    def __init__(self, id_: str) -> None:
        self.id    = id_
        self.exons = list()

    def add_exons(self, start: int, end: int) -> None:
        self.exons.append((start - 1, end)) # set the indices to be 0-based

    def get_exons(self) -> list[tuple[int, int]]:
        self.exons.sort(key = lambda x: x[0])
        return self.exons
    
class Gene:

    __slots__ = ("id", "chr", "transcripts", "_coding")

    def __init__(self, id_: str, chrom: str) -> None:
        self.id          = id_
        self.chr         = chrom
        self.transcripts = list()
        self._coding     = False

    def is_coding(self) -> bool:
        return self._coding
    
    def set_as_coding(self) -> None:
        self._coding = True

    def add_transcript(self, transcript_id: str, start: int, end: int) -> int:
        if (self._coding):
            return 0 # this gene will not be used
        
        for i in range(len(self.transcripts)):
            if (self.transcripts[i].id == transcript_id):
                self.transcripts[i].add_exons(start, end)
                return 0
            
        #
        # if not executed, add the transcript
        #
        transcript = Transcript(transcript_id)
        transcript.add_exons(start, end)
        self.transcripts.append(transcript)
        return
    
    def get_transcripts(self) -> list[Transcript]:
        return self.transcripts
    
#############################################################################################################################
#                                                  Command Line Options                                                     #
#############################################################################################################################

def validate_query(path_: str, org_: str, annID_: str) -> int:
    """validate that this org and the annotation exist in the species database"""

    annMap = f"{path_}/annotation_integration.map"

    if (os.path.isfile(annMap) == False):
        msg = f"Unable to find the annotation map in the species database provided"
        sys.exit(msg)

    fh    = open(annMap, 'r')
    ln    = 0
    found = False

    for line in fh:
        ln += 1
        if (len(line) == 0 or line[0] == '#'):
            continue
        record = line.split('\t')

        if (len(record) != 6):
            msg = f"Expected 6 columns in annotation map. Found {len(record)} at line {ln}"
            sys.exit(msg)

        if (org_ == record[0] and annID_ == record[2]):
            found = True
            break 

    fh.close()

    if (found == False):
        msg = f"Unable to find integration {annID_} for {org_} in {path_}"
        sys.exit(msg)

    return 0

def set_arguments(args) -> int:
    """helper function to check the required parameters for this feature"""

    global path, orgID, annID, ref_fasta

    #
    # species db check
    #
    if (args.path == ''):
        msg = "--path is required"
        sys.exit(msg)
    elif (os.path.isdir(args.path) == False):
        msg = f"Could not locate path {args.path}/"
        sys.exit(msg)
    elif (args.path[-1] == '/'):
        args.path = args.path[:-1]

    if (args.id == ''):
        msg = "--id is required"
        sys.exit(msg)

    if (args.ann_id == ''):
        msg = "--ann-id is required" 
        sys.exit(msg)

    if (args.ref == ''):
        msg = "--ref is required"
        sys.exit(msg)
    elif (os.path.isfile(args.ref) == False):
        msg = f"Could not locate reference genome: {args.ref}"
        sys.exit(msg)

    #
    # check the annotation_integration.map
    # to verify that this query exists
    #
    validate_query(args.path, args.id, args.ann_id)

    # set the global variables
    path      = args.path
    orgID     = args.id
    annID     = args.ann_id
    ref_fasta = args.ref
   
    return 0
    
def parse_command_line() -> int:
    """helper function to get the user's arguments to ensure a proper start"""
    
    desc   = "Generate orthogroup-specific fasta files from a Synolog analysis"
    parser = argparse.ArgumentParser(description=desc)

    #
    # help messages to shorten argument lines
    #
    phelp = "Specify the base path to the Synolog data cache."
    ihelp = "A short alphanumeric label/ID to represent a species within Synolog."
    ahelp = "A short label to represent this genome structure"
    rhelp = "Reference genome in fasta format to extract non-coding regions from"

    parser.add_argument("--path",   help=phelp, type=str, default='')
    parser.add_argument("--id",     help=ihelp, type=str, default='')
    parser.add_argument("--ann-id", help=ahelp, type=str, default='')
    parser.add_argument("--ref",    help=rhelp, type=str, default='')

    args = parser.parse_args()
    
    # now set the global variables
    set_arguments(args)

    return 0

#############################################################################################################################
#                                         Helerp Functions (Annotations)                                                    #
#############################################################################################################################

def get_annotation() -> str:
    """grab the annotation file for this org"""

    global path, orgID, annID

    gidir  = f"{path}/genome_integrations"

    if (os.path.isdir(gidir) == False):
        msg = f"genome_integrations directory not found in {path}/\n" + \
               "Was this removed?"
        sys.exit(msg)

    # construct orgID.IntegrationID
    target = orgID + '.' + annID

    anns = glob(f"{gidir}/{target}/*.g*.gz", recursive=True)

    if (len(anns) == 0):
        msg = f"Unable to find an annotation file for {orgID} in {annID}/"
        sys.exit(msg)

    ann = ''
    for an_ in anns:
        annFile = os.path.basename(an_)
        if (annFile.startswith(target)):
            ann = an_
            break

    if (ann == ''):
        msg = f"Unable to find an anotation file matching species id {orgID} in {annID}"
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

def get_nc_ids() -> int:
    """grab the transcripts ids for the genes in target org integration"""

    global orgID

    # get the org's annotation
    ann = get_annotation()

    # now to loop through and parse the annotation
    fh      = gzip.open(ann, "rt") if ann.endswith(".gz") else open(ann, 'r')
    lineNum = 0
    feats   = ["exon", "cds"] # features
    genes   = dict()

    for line in fh:
        lineNum += 1
        if (len(line) == '' or line[0] == '#'):
            continue
        fields = line.split('\t')
        if (fields[2].lower() == "gene"):
            gene_id = get_gene_id(fields[8])
            if (gene_id == ''):
                msg = f"Failed to find a gene_id/ID at line {lineNum} in {ann} for {orgID}"
                sys.exit(msg)
            if (gene_id not in genes):
                genes[gene_id] = Gene(gene_id, fields[0])
            continue
        
        feat = fields[2].lower()
        if (feat not in feats):
            continue

        attrbMap = make_attrb_map(fields[8])

        if ("gene_id" not in attrbMap):
            if ("ID" not in attrbMap):
                msg = f"Failed to find a gene_id/ID at line {lineNum} in {ann} for {orgID}"
                sys.exit(msg)
            else:
                gene_id = attrbMap["ID"]
        else:
            gene_id = attrbMap["gene_id"]
        
        if (feat == "cds"):
            # this gene is a protein coding gene
            genes[gene_id].set_as_coding()
            continue
        elif (feat == "exon"):
            start   = int(fields[3])
            end     = int(fields[4])
            tran_id = ''
            if ("transcript_id" in attrbMap):
                tran_id = attrbMap["transcript_id"]
                if ("transcript_version" in attrbMap):
                    tran_vs = attrbMap["transcript_version"]
                    tran_id = tran_id + '.' + tran_vs
            if (tran_id == ''):
                msg = f"Failed to find a transcript ID at line {lineNum} in {ann} for {orgID}"
                sys.exit(msg)
            genes[gene_id].add_transcript(tran_id, start, end)
            
    fh.close()

    total = len(genes)
    count = 0
    global chroms

    # now to populate the chroms dictionary with the transcripts
    for gene in genes.values():
        if gene.is_coding():
            continue
        count += 1
        transcripts = gene.get_transcripts()
        chroms[gene.chr].extend(transcripts)

    per = round((count / total) * 100, 2)
    msg = f"Processed {total} genes. {count} ({per}%) were found to be non-coding"
    print(msg, file=sys.stderr)

    if (count == 0):
        sys.exit(0)

    return 0

#############################################################################################################################
#                                               Helper Functions (Fastas)                                                   #
#############################################################################################################################

def make_outname() -> str:
    """helper function to generate the output name in the species database"""

    global path, ref_fasta, orgID, annID

    bname = os.path.basename(ref_fasta)

    #
    # remove suffixes
    #
    if (bname.endswith(".gz")):
        bname = bname[:-1]

    for suffix in [".fna", ".fa", ".fasta"]:
        if (bname.endswith(suffix)):
            bname = bname.replace(f"{suffix}", '')
            break

    # if this cache was created without this component
    outDir = f"{path}/genome_nc_genes"

    if (os.path.isdir(outDir) == False):
        os.mkdir(outDir, mode=0o700)

    # query specific directory
    outDir = outDir + '/' + orgID + '.' + annID
    if (os.path.isdir(outDir) == False):
        os.mkdir(outDir, mode=0o700)

    out = f"{outDir}/{bname}_rna_from_genomic.fna.gz"
    
    return out

def make_symlink() -> int:
    """helper function to create the symlink to the new fasta file"""

    global path, orgID, annID, outName

    outDir = f"{path}/genome_nc_genes/{orgID}.{annID}"
    os.chdir(outDir)

    lnName = f"{orgID}.{annID}.nc.fa.gz"
    os.symlink(os.path.basename(outName), lnName)

    return 0

def get_header_id(header: str) -> str:
    """return the first characters before any spacing"""

    id_ = ''
    idx = header.find(' ')

    if (idx == -1):
        id_ = header[1:]
    else:
        id_ = header[1:idx]

    return id_

def get_org_sequences() -> int:
    """extract the exonic sequences"""

    global chroms, ref_fasta, outName

    # parse and process
    outName = make_outname()
    fh      = gzip.open(ref_fasta, "rt") if ref_fasta.endswith(".gz") else open(ref_fasta, 'r')
    outfh   = gzip.open(outName, "wt")
    lineNum = 0
    curRec  = ''
    curSeq  = list()
    count   = 0

    for line in fh:

        lineNum += 1

        if (len(line) == 0 or line[0] == '#'):
            continue
        line = line.strip()

        if (line[0] == '>'):
            if (curRec == ''):
                curRec = get_header_id(line)
                if (curRec not in chroms):
                    curRec = ''
            else:
                seq = ''.join(curSeq)
                for transcript in chroms[curRec]:
                    header = '>' + transcript.id + '\n'
                    exons  = transcript.get_exons()
                    rnas   = list()
                    # exons are tuples of start:end positions
                    for exon in exons:
                        rna = seq[exon[0]:exon[1]]
                        rnas.append(rna)
                    rna = ''.join(rna)
                    rna = textwrap.fill(rna, width=60) + '\n'
                    outfh.write(header)
                    outfh.write(rna)
                    count += 1
                chroms[curRec].clear()

                # update
                curSeq.clear()
                seq    = ''
                curRec = get_header_id(line)
                if (curRec not in chroms):
                    curRec = ''

        elif (curRec != ''):
            curSeq.append(line)

    fh.close()

    if (curRec != ''):
        seq = ''.join(curSeq)
        for transcript in chroms[curRec]:
            header = '>' + transcript.id + '\n'
            exons  = transcript.get_exons()
            rnas   = list()
            # exons are tuples of start:end positions
            for exon in exons:
                rna = seq[exon[0]:exon[1]]
                rnas.append(rna)
            rna = ''.join(rna)
            rna = textwrap.fill(rna, width=60) + '\n'
            outfh.write(header)
            outfh.write(rna)
            count += 1
        chroms[curRec].clear()

    outfh.close()

    missing = list()
    for chrom, transcripts in chroms.items():
        if (len(transcripts) > 0):
            missing.append(chrom)

    print(f"Wrote {count} transcript entries", file=sys.stderr)

    if (len(missing) > 0):
        msg = f"Unable to find genes on the following sequences: " + \
              ','.join(missing)
        print(msg, file=sys.stderr)

    # now to symlink the generated file
    make_symlink()

    return 0

#############################################################################################################################
#                                               Helper Functions (BLAST)                                                    #
#############################################################################################################################

def make_blast_db() -> int:
    """construct a blast db for the org elements"""

    global outName, orgID, annID

    bDir = os.path.dirname(outName)

    os.chdir(bDir)



#############################################################################################################################
#                                                    MAIN                                                                   #
#############################################################################################################################

def main() -> int:
    
    # retrieve the arguments
    parse_command_line()

    # get the non-coding gene locations
    get_nc_ids()

    # parse the reference genome for the rnas
    get_org_sequences()
    
    return 0
        
if __name__ == '__main__':
    main()
