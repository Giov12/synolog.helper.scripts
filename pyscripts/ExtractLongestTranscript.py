#!/bin/env python3

import argparse
import os
import sys
import textwrap
import gzip

ann         = ''
fna         = ''
genes       = dict()
transcripts = dict()

class Transcript:

    __slots__ = ("id", "seq")

    def __init__(self, id_: str) -> None:
        self.id    = id_
        self.seq   = ''

class Gene:

    __slots__ = ("id", "transcripts", "longest_trans", "longest_size")

    def __init__(self, id_: str) -> None:
        self.id            = id_
        self.transcripts   = dict()
        self.longest_trans = ''
        self.longest_size  = 0

    def add_transcript(self, transcript_id: str) -> int:

        self.transcripts[transcript_id] = Transcript(transcript_id)
  
        return 0
    
    def update(self, transcript_id: str, seq: str) -> int:

        seq_len = len(seq)
        if (self.longest_size < seq_len):
            self.longest_trans = transcript_id
            self.longest_size  = seq_len
            self.transcripts[transcript_id].seq = seq

        return 0
    
    def get_longest(self) -> Transcript:
        if (self.longest_trans == ''):
            print("No transcript found for", self.id)
            return None
        return self.transcripts[self.longest_trans]

def set_arguments() -> int:
    """get & set the arguments"""

    d = "Filter for the longest transcript given an annotations file"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-a", "--ann", help="a gtf or gff3 file", required=True, type=str)
    parser.add_argument("-f", "--fasta", help="a fasta file containing the transcripts/isoforms", required=True, type=str)

    args = parser.parse_args()

    assert os.path.isfile(args.ann), f"Could not locate file: {args.ann}"
    assert os.path.isfile(args.fasta), f"Could not locate file: {args.fasta}"

    global ann, fna
    ann = args.ann
    fna = args.fasta
    
    return 0

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

    amap   = dict()
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
        amap[key] = val

    return amap

def get_genes() -> int:
    """grab the transcripts ids for the genes in target org integration"""

    global ann, genes, transcripts

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
        feat   = fields[2].lower()
        
        if (feat == "gene"):
            gene_id = get_gene_id(fields[8])
            if (gene_id == ''):
                msg = f"Failed to find a gene_id/ID at line {lineNum} in {ann}"
                sys.exit(msg)
            if (gene_id not in genes):
                genes[gene_id] = Gene(gene_id)
            continue
        
        elif (feat not in feats):
            continue

        attrbMap = make_attrb_map(fields[8])

        if ("gene_id" not in attrbMap):
            if ("ID" not in attrbMap):
                msg = f"Failed to find a gene_id/ID for the following line: {line}"
                print(msg, end='')
                continue
            else:
                gene_id = attrbMap["ID"]
        else:
            gene_id = attrbMap["gene_id"]
        
        tran_id = ''
        if ("protein_id" in attrbMap):
            tran_id = attrbMap["protein_id"]
            if ("protein_version" in attrbMap):
                tran_vs = attrbMap["protein_version"]
                tran_id = tran_id + '.' + tran_vs
        elif ("transcript_id" in attrbMap):
            tran_id = attrbMap["transcript_id"]
            if ("transcript_version" in attrbMap):
                tran_vs = attrbMap["transcript_version"]
                tran_id = tran_id + '.' + tran_vs
        elif ("Parent" in attrbMap):
            tran_id = attrbMap["Parent"]

        if (tran_id == ''):
            msg = f"Failed to find a transcript ID at line {lineNum} in {ann}"
            sys.exit(msg)
        if (gene_id not in genes):
            tmp_id = gene_id.split('.')[0]
            if (tmp_id in genes):
                gene_id = tmp_id
        if (gene_id not in genes and gene_id.startswith("unassigned_gene_")):
            print("Skipping entries for the following line:", line, end='')
            continue
        genes[gene_id].add_transcript(tran_id)
        transcripts[tran_id] = gene_id
            
    fh.close()

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

def get_sequences() -> int:
    """extract the exonic sequences"""

    global genes, transcripts, fna

    # parse and process
    fh      = gzip.open(fna, "rt") if fna.endswith(".gz") else open(fna, 'r')
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
                if (curRec not in transcripts):
                    print("Not in annotation:", curRec)
                    curRec = ''
            else:
                seq = ''.join(curSeq)
                genes[transcripts[curRec]].update(curRec, seq)
                # update
                curSeq.clear()
                seq    = ''
                curRec = get_header_id(line)
                if (curRec not in transcripts):
                    print("Not in annotation:", curRec)
                    curRec = ''
        elif (curRec != ''):
            curSeq.append(line)

    fh.close()

    if (curRec != ''):
        seq = ''.join(curSeq)
        genes[transcripts[curRec]].update(curRec, seq)

    return 0

def write_fasta() -> None:
    """write the longest transcripts to an out file"""

    global genes, fna

    fname = os.path.basename(fna)
    fname = f"Filtered.New.{fname}"

    if (fname.endswith(".gz")):
        fh = gzip.open(fname, "wt")
    else:
        fname += ".gz"
        fh     = gzip.open(fname, 'wt')

    for gene in genes.values():
        trans = gene.get_longest()
        if (trans == None):
            continue
        fh.write('>' + trans.id + '\n')
        fh.write(textwrap.fill(trans.seq, width=60))
        fh.write('\n')
    fh.close()

def main() -> int:
    """Filter for the longest transcripts/isoforms"""

    # get arguments
    set_arguments()

    # get the genes
    get_genes()

    # get the sequences
    get_sequences()

    # now write out the longest sequences
    write_fasta()

    return 0

if __name__ == "__main__":
    main()
