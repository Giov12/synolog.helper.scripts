#!/bin/env python3

import argparse
import os

ofile = ''
        
def get_arguments() -> int:
    """get the argument"""

    global ofile

    d = "Some python code to reformat the synolog orthologs.tsv to a format compatible for benchmarking"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-o", "--orthologs", help="orthologs.tsv file generated from synolog", required=True, type=str)

    args  = parser.parse_args()
    ofile = args.orthologs

    assert os.path.isfile(ofile), f"Could not locate file: {ofile}"
    
    return 0

def reformat() -> int:
    """single function that will parse and write as it goes"""

    global ofile

    fh   = open(ofile, 'r')
    ofh  = open("benchmarking.orthologs.txt", 'w')
    cur  = ''
    mems = list()

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.split('\t')
        grp    = fields[0] # this is safe for orthologs.tsv regardless
        geneID = fields[5] # of tandem duplications or not

        if (cur == ''): # first value
            cur = grp
        if (cur == grp):
            mems.append(geneID)
        else:
            outline = ", ".join(mems)
            outline = cur + ": " + outline + '\n'
            ofh.write(outline)
            mems.clear()
            cur = grp
            mems.append(geneID)

    if (len(mems) > 0):
        outline = ", ".join(mems)
        outline = cur + ": " + outline + '\n'
        ofh.write(outline)

    fh.close()
    ofh.close()

    return 0

def main() -> int:
    """entry point to this application"""

    # get the single input argument
    get_arguments()

    # reformat as we go
    reformat()

    return 0

if __name__ == "__main__":
    main()
