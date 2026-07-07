#!/bin/env python3

import argparse
import os
import gzip
from collections import defaultdict

tbl = ''

def get_arguments() -> int:
    """get the input file"""

    global tbl

    parser = argparse.ArgumentParser(
        description="Convert an cmScan table format 2 file to a gtf file")
    parser.add_argument(
        "-t", "--table", help="table file [required]", required=True)

    args = parser.parse_args()

    assert os.path.isfile(args.table), f"Could not find {args.table}"

    tbl = args.table

    return 0

def make_gtf():
    """create a gtf file from the table"""

    global tbl

    bname = os.path.basename(tbl)
    ofh   = open(bname + ".gtf", 'w')
    fh    = gzip.open(tbl, "rt") if bname.endswith(".gz") else open(tbl, 'r')
    cnts  = defaultdict(int)

    # now parse the gff
    for lineNum, line in enumerate(fh):
        if (len(line) == 0) or (line[0] == '\t') or (line[0] == '#'):
            continue
        fields = line.strip('\n').replace('\t', ' ').split(' ')
        idx    = 1 # skip first column
        jdx    = 0
        type_  = ''
        feat_  = ''
        desc   = ''
        chrom  = ''
        start  = -1
        end    = -1
        strand = ''
        nums   = 0

        while (idx < len(fields)):
            while (fields[idx] == ''):
                idx += 1
            if (type_ == ''):
                type_ = fields[idx]
            elif (feat_ == ''):
                feat_ = fields[idx]
            elif (chrom == ''):
                chrom = fields[idx]
            elif (fields[idx] == '-'):
                jdx = idx
                if (end > -1 and strand == ''):
                    strand = fields[idx]
            elif (fields[idx] == '+'):
                strand = fields[idx]
            if (fields[idx].isnumeric()):
                nums += 1
                if (start == -1 and nums == 3):
                    start = int(fields[idx])
                elif (end == -1 and nums == 4):
                    end = int(fields[idx])
            idx += 1
        desc = ' '.join(fields[jdx + 1:])
        cnt  = cnts[type_]
        cnts[type_] += 1
        gene_id = f"ncRNA-{type_}-{cnt}"
        desc = f"gene_id \"{gene_id}\"; desc \"{desc}\""
        outline = [chrom, "cmscan", "gene", str(start), str(end), '.', strand, '.', desc + '\n']
        trans   = outline.copy()
        trans[2] = "transcript"
        trans[8] = desc + f"; transcript_id {gene_id}.t1;\n"
        exon     = trans.copy()
        exon[2]  = "exon"
        ofh.write('\t'.join(outline))
        ofh.write('\t'.join(trans))
        ofh.write('\t'.join(exon))
        

    fh.close()
    ofh.close()

def main():
    """start the program"""

    # get input
    get_arguments()

    # write output
    make_gtf()
    
if __name__ == "__main__":
    main()
