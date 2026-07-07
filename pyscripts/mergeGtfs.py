#!/bin/env python3

import argparse
import os
import sys
import gzip
from collections import defaultdict

genes = dict()
gtfs  = list()

class Gene:

    __slots__ = ("chrom", "start", "end", "lines")

    def __init__(self, chrom: str, start: int, end: int) -> None:
        self.chrom = chrom
        self.start = start
        self.end   = end
        self.lines = list()

    def add_line(self, line: str) -> int:
        self.lines.append(line)
    
    def get_lines(self) -> list[str]:
        return self.lines
    
def set_arguments() -> int:
    """get & set the arguments"""

    global gtfs

    d = "Merge a set of gtf files"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-g", "--gtfs", help="list of gtf files", required=True, nargs='*')

    args = parser.parse_args()

    assert len(args.gtfs) > 1, "Program requires at least 2 gtf files"
    for g in args.gtfs:
        assert os.path.isfile(g), f"Could not locate {g}"

    gtfs = args.gtfs
    
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


def parse_gtf(gtf: str) -> int:
    """grab all the genes in this file"""

    global genes

    # now to loop through and parse the annotation
    fh      = gzip.open(gtf, "rt") if gtf.endswith(".gz") else open(gtf, 'r')
    lineNum = 0

    for line in fh:
        lineNum += 1
        if (len(line) == '' or line[0] == '#'):
            continue
        fields = line.split('\t')
        gene_id = get_gene_id(fields[8])
        if (gene_id == ''):
            msg = f"Failed to find a gene_id/ID at line {lineNum} in {gtf}"
            sys.exit(msg)
        if (fields[2] == "gene"):
            genes[gene_id] = Gene(fields[0], int(fields[3]), int(fields[4]))
        genes[gene_id].add_line(line)    
            
    fh.close()

    return 0

def merge() -> int:
    """write a combined file for all the genes parsed"""

    global genes

    chroms = defaultdict(list)
    outfh  = open("Combined.gtf", 'w')

    for gene in genes.values():
        chroms[gene.chrom].append(gene)

    keys = sorted(chroms.keys())

    for chrom in keys:
        genes_list = chroms[chrom]
        genes_list.sort(key= lambda ge: (ge.start, ge.end))
        for gene in genes_list:
            outlines = gene.get_lines()
            for line in outlines:
                outfh.write(line)
    
    outfh.close()

    return 0

def main() -> int:
    """Filter for the longest transcripts/isoforms"""

    # get arguments
    set_arguments()

    # get the genes
    global gtfs
    for gtf in gtfs:
        parse_gtf(gtf)

    # combine
    merge()

    return 0

if __name__ == "__main__":
    main()
