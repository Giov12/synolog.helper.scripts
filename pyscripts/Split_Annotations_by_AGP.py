#!/bin/env python3

import argparse
import os
import sys
import gzip
from collections import defaultdict

ann    = ''
agp    = ''
genes  = dict()
chroms = defaultdict(list)

class Contig:

    __slots__ = ("parent", "size", "start", "end", "name", "_reversed")
    def __init__(self, parent: str, size: int, start: int, end: int, name: str, 
                 orientaiton: str) -> None:
        self.parent    = parent
        self.size      = size
        self.start     = start
        self.end       = end
        self.name      = name
        self._reversed = (orientaiton == '-') 

    def is_reversed(self) -> bool:
        return self._reversed

def set_arguments() -> int:
    """get & set the arguments"""

    d = "Split a gff3 file based on an agp file"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-g", "--gff", help="a gtf or gff3 file", required=True, type=str)
    parser.add_argument("-a", "--agp", help="an AGP describing the structure of the genome", required=True, type=str)

    args = parser.parse_args()

    assert os.path.isfile(args.agp), f"Could not locate file: {args.agp}"
    assert os.path.isfile(args.gff), f"Could not locate file: {args.gff}"

    global ann, agp
    agp = args.agp
    ann = args.gff
    
    return 0

def parse_agp() -> int:
    """add the contigs to the chromosomes list"""

    global agp, chroms

    fh = gzip.open(agp, "rt") if agp.endswith(".gz") else open(agp, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.split('\t')

        # skip gaps
        if (fields[4] == 'N'):
            continue
        chrom  = fields[0]
        start  = int(fields[1])
        end    = int(fields[2])
        name   = fields[5]
        size   = int(fields[7])
        orien  = fields[8]
        chroms[chrom].append(Contig(chrom, size, start, end, name, orien))

    fh.close()

    return 0

def make_outname() -> str:
    """create a output name for this program"""

    global ann

    fname = os.path.basename(ann)
    if (fname.endswith(".gz")):
        fname = fname[:-3]
    fname = "contiged-" + fname
    return fname

def rewrite() -> int:
    """parse the annotations gff/gtf file"""

    global ann, chroms

    # now to loop through and parse the annotation
    fh      = gzip.open(ann, "rt") if ann.endswith(".gz") else open(ann, 'r')
    ofh     = open(make_outname(), 'w')
    lineNum = 0
    skips   = 0
    skip    = False
    
    for line in fh:
        lineNum += 1
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields  = line.split('\t')
        chrom   = fields[0]
        start   = int(fields[3])
        end     = int(fields[4])
        recType = fields[2].lower()
        contigs = chroms[chrom]
        contig  = None
        left    = 0
        right   = len(contigs) - 1

        if (recType == "gene"):
            skip = False
        elif (skip):
            continue

        while (left <= right):   
            # candidate contig
            mid   = (right + left) // 2
            candi = contigs[mid]
            if (candi.end < start):
                left = mid + 1
            elif (candi.start > end):
                right = mid - 1
            elif (candi.start <= start <= end <= candi.end):
                contig = candi
                break
            elif (recType == "gene"):
                skip   = True
                skips += 1
                break
            
        if (contig == None):
            continue

        newStart = start - contig.start
        newEnd   = end   - contig.start

        if (contig.is_reversed()):
            newStart = contig.size - newEnd
            newEnd   = contig.size - newStart
            Strand   = '+' if fields[6] == '-' else '-'
        else:
            newStart += 1
            newEnd   += 1
            Strand    = fields[6]

        fields[0] = contig.name
        fields[3] = str(newStart)
        fields[4] = str(newEnd)
        fields[6] = Strand

        ofh.write('\t'.join(fields))
            
    fh.close()
    ofh.close()

    print("Skipped", skips, "gene records not contained in a single contig")

    return 0

def main() -> int:
    """Entry point to this application"""

    # get arguments
    set_arguments()

    # get the contigs
    parse_agp()

    # create the new annotations
    rewrite()

    return 0

if __name__ == "__main__":
    main()
