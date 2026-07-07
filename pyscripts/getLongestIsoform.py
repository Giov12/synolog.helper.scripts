#!/bin/env python3

import argparse
import os
import textwrap
import gzip
import sys

"""
filter a protein fasta file for the longest transcript
"""

class Seq:
    __slots__ = ("header", "seq")
    def __init__(self, header: str, seq: str) -> None:
        self.header = header
        self.seq    = seq

    def __len__(self) -> int:
        return len(self.seq)

class GENE:

    __slots__ = ("id", "records")
    
    def __init__(self, id_: str) -> None:
        self.id       = id_
        self.records  = list()

class Org:

    __slots__ = ("id", "genes", "isoforms")

    def __init__(self, id_: str) -> None:
        self.id       = id_
        self.genes    = dict()
        self.isoforms = dict() # prot id -> gene id

    def add_gene(self, gene_id: str) -> None:
        self.genes[gene_id] = GENE(gene_id)

    def add_isoform(self, gene_id: str, prot_id: str) -> None:
        self.isoforms[prot_id] = gene_id

    def add_record(self, id_: str, header: str, seq: str) -> None:
        gene_id = self.isoforms[id_]  # id_ is a protein id
        rec     = Seq(header, seq)
        gene    = self.genes[gene_id]

        if (len(gene.records) == 0):
            gene.records.append(rec)
        elif(len(gene.records[0]) < len(rec)):
            # replace the current id
            gene.records[0] = rec
        
def get_arguments() -> tuple[str, str, str]:
    """get the arguments"""

    d = "Some python code to filter protein sequences to the longest transcript"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-g", "--gtf", help="gtf file", required=True, type=str)
    parser.add_argument("-f", "--fasta", help="fasta file", required=True, type=str)
    parser.add_argument("-o", "--org", help="name of org", required=True, type=str)

    args = parser.parse_args()
    faa  = args.fasta
    gtf  = args.gtf
    id_  = args.org

    assert os.path.isfile(faa), f"Could not locate file: {faa}"
    assert os.path.isfile(gtf), f"Could not locate file: {gtf}"
    
    return (faa, gtf, id_)

def get_prot_id(header: str) -> str:
    """return the first characters before any spacing"""

    id_ = ''
    idx = header.find(' ')

    if (idx == -1):
        id_ = header[1:]
    else:
        id_ = header[1:idx]

    return id_

def get_gene_id(attrb: str) -> str:
    """return the gene_id for this record"""

    fields  = attrb.split(';')
    gene_id = ''
    id_     = ''

    if (len(fields) == 1):
        if ("gene_id" in fields[0]):
            gene_id = fields[0].replace("gene_id", '')
            gene_id = gene_id.strip(' "\n')
        elif ("ID=" in fields[0]):
            gene_id = fields[0].replace("ID=", '')
            gene_id = gene_id.strip(' "\n')
        else:
            gene_id = fields[0].strip(' "\n')
        return gene_id
    
    for field in fields:
        field = field.strip(' "')
        if ((field.startswith("gene_id") == False) and (field.startswith("ID") == False)):
            continue
        subfields = field.strip(", \"=").split(' ')
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

def get_org_proteins(org: Org, faa: str) -> int:
    """grab the isoforms from the fasta file for the target genes"""

    # parse and process
    fh      = gzip.open(faa, "rt") if faa.endswith(".gz") else open(faa, 'r')
    lineNum = 0
    curIso  = ''
    curHdr  = '' # preserve the header
    curSeq  = list()

    for line in fh:

        lineNum += 1
        if (len(line) == 0 or line[0] == '#'):
            continue
        line = line.strip()
        if (line[0] == '>'):
            if (curIso == ''):
                isoID = get_prot_id(line)
                if (isoID in org.isoforms):
                    curIso = isoID
                    curHdr = line
            else:
                seq = ''.join(curSeq)
                org.add_record(curIso, curHdr, seq)
                # update
                curSeq.clear()
                seq   = ''
                isoID = get_prot_id(line)
                if (isoID in org.isoforms):
                    curIso = isoID
                    curHdr = line
                else:
                    curIso = ''
                    curHdr = ''
        elif (curIso != ''):
            curSeq.append(line)

    fh.close()

    if (curIso != ''):
        seq = ''.join(curSeq)
        org.add_record(curIso, curHdr, seq)

    return 0

def get_org_protein_ids(org: Org, gtf: str) -> int:
    """grab the protein ids for the genes in org"""

    # now to loop through and parse the annotation
    fh      = gzip.open(gtf, "rt") if gtf.endswith(".gz") else open(gtf, 'r')
    lineNum = 0
    curGene = ''

    for line in fh:
        lineNum += 1
        if (len(line) == '' or line[0] == '#'):
            continue
        fields = line.split('\t')
        if (fields[2].lower() == "gene"):
            gene_id = get_gene_id(fields[8])
            if (gene_id == ''):
                msg = f"Failed to find a gene_id/ID at line {lineNum} in {gtf} for {org.id}"
                print(msg, file=sys.stderr)
            else:
                curGene = gene_id
                org.add_gene(curGene)
            continue
        elif (curGene == '' or fields[2].lower() != "cds"):
            continue
        attrbMap = make_attrb_map(fields[8])

        if ("protein_id" in attrbMap):
            prot_id = attrbMap["protein_id"]
            if ("protein_version" in attrbMap):
                prot_vs = attrbMap['protein_version']
                prot_id = prot_id + '.' + prot_vs
            org.add_isoform(curGene, prot_id)
        elif ("transcript_id" in attrbMap):
            tran_id = attrbMap["transcript_id"]
            if ("transcript_version" in attrbMap):
                tran_vs = attrbMap["transcript_version"]
                tran_id = tran_id + '.' + tran_vs
            org.add_isoform(curGene, tran_id)

    fh.close()

    return 0

def write_proteins(org: Org) -> int:

    out  = f"{org.id}.def.faa"
    err  = f"{org.id}.missing.genes.txt"
    fh   = open(out, 'w')
    eh   = open(err, 'w')
    cnt  = 0
    tot  = 0

    for gene_id, gene in org.genes.items():
        tot += 1
        if (len(gene.records) == 0):
            cnt += 1
            line = f"No isoforms found for {gene_id}\n"
            eh.write(line)
            continue
        rec = gene.records[0]
        fh.write(f">{gene_id}\n")
        fh.write(textwrap.fill(rec.seq, width=60))
        fh.write('\n')

    fh.close() 
    eh.close()

    per = round((cnt / tot) * 100, 2)
    print(f"Missing {cnt} out of {tot} ({per}%) genes")
    return 0

def main() -> int:
    """entry point to this application"""

    # get arguments
    faa, gtf, id_ = get_arguments()

    # create an org object & get the genes
    org = Org(id_)
    get_org_protein_ids(org, gtf)

    # now add the prot seqs to each gene
    get_org_proteins(org, faa)

    # write the protein seqs using the gene IDs
    write_proteins(org)

    return 0

if __name__ == "__main__":
    main()
