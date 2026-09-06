#!/bin/env python3

import argparse
import os
import glob
import sys

homologs = dict()

class Pair:
    __slots__ = ("chicken", "turtle")

    def __init__(self):
        self.chicken = ''
        self.turtle  = ''
  
def get_arguments() -> tuple[str, str]:
    """get the arguments"""

    d = "Some python code to find the chicken gene symbols for the turtle orthogroups of interest"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-H", "--homologs",  help="directory containing homolog files of with model organisms", required=True)
    parser.add_argument("-o", "--orthologs", help="orthologs.tsv file from synolog_info.py", required=True)

    args  = parser.parse_args()
    hDir  = args.homologs
    ofile = args.orthologs

    assert os.path.isdir(hDir),   f"Could not locate directory: {hDir}"
    assert os.path.isfile(ofile), f"Could not locate file: {ofile}"

    return (hDir, ofile)

def parse_homologs(hFile: str) -> int:
    """parse a homologs file"""

    global homologs

    fh      = open(hFile, 'r')
    chicken = "ggal"
    turtle  = "Gevg"

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.split('\t')
        orgA   = fields[2]
        orgB   = fields[3]
        idA    = fields[4]
        nameA  = fields[5]
        idB    = fields[9]
        nameB  = fields[10]
        useA   = False
        if (orgA.startswith(chicken) or orgA.startswith(turtle)):
            key  = f"{orgB}:{idB}"
            useA = True
        else:
            key = f"{orgA}:{idA}"

        if (key not in homologs):
            homologs[key] = Pair()
        
        if (useA):
            if (orgA.startswith(chicken)):
                homologs[key].chicken = nameA
            else:
                homologs[key].turtle  = nameA
        else:
            if (orgB.startswith(chicken)):
                homologs[key].chicken = nameB
            else:
                homologs[key].turtle  = nameB
    
    fh.close()

    return 0

def get_homologs(hDir: str) -> int:
    """get the homolog files"""

    if (hDir[-1] == '/'):
        hFiles = glob.glob(f"{hDir}*homologs.tsv")
    else:
        hFiles = glob.glob(f"{hDir}/*homologs.tsv")

    if (len(hFiles) == 0):
        msg = f"No *homologs.tsv files found in {hDir}"
        sys.exit(msg)

    for hFile in hFiles:
        parse_homologs(hFile)

    return 0

def add_gene_names(ofile: str) -> int:
    """parse the text file for orthogroups of interest"""

    global homologs

    fh   = open(ofile, 'r')
    out  = "Named.Helpers-" + os.path.basename(ofile)
    out  = out.replace(".txt", ".tsv")
    ofh  = open(out, 'w')
    out  = out.replace("Named", "Group")
    ofh2 = open(out, 'w')
    cur  = "-1"
    nams = set()

    for line in fh:
        if (len(line) == 0 or line[0].isdigit() == False):
            continue
        fields = line.split('\t')
        id_    = fields[0]
        org    = fields[4]
        gid    = fields[5]
        name   = fields[6]

        if (cur == "-1"):
            cur = id_
        
        if (id_ != cur):
            outs = ','.join(nams)
            ofh2.write(f"{cur}\t{outs}\n")
            nams.clear()
            cur = id_

        key = f"{org}:{gid}"
        if (name == '' or name.startswith("LOC")):
            if (key in homologs):
                if (homologs[key].turtle != ''):
                    name = homologs[key].turtle
                elif (homologs[key].chicken != ''):
                    name = homologs[key].chicken

        outline = f"{id_}\t{org}\t{gid}\t{name}\n"
        ofh.write(outline)
        if (name != '' and name.startswith("LOC") == False):
            nams.add(name)

    fh.close()
    ofh.close()

    outs = ','.join(nams)
    ofh2.write(f"{cur}\t{outs}\n")
    nams.clear()
    ofh2.close()

    return 0

def main() -> int:
    """entry point to this little helper code"""

    # get the input file & directory
    hDir, ofile = get_arguments()

    # load all the homologs
    get_homologs(hDir)

    # add names
    add_gene_names(ofile)

    return 0

if __name__ == "__main__":
    main()
