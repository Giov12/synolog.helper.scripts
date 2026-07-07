#!/bin/env python3

import argparse
import os

total = 0

class entry:
    def __init__(self, fields: list[str]) -> int:
        self.id    = fields[0]
        self.org   = fields[4]
        self.type  = fields[7]
        self.gid   = fields[5]
        self.count = int(fields[8]) 
    def __str__(self):
        return f"{self.id}\t{self.org}\t{self.gid}\t{self.type}\t{self.count}"

        
def get_arguments() -> str:
    """get the arguments"""

    d = "Some python code to parse the orthologs.tsv from the turtles synteny run"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-o", "--orthologs",  help="orthologs.tsv from the turtles synteny analysis", required=True)

    args = parser.parse_args()
    ofile = args.orthologs

    assert os.path.isfile(ofile), f"Could not locate file: {ofile}"

    return ofile

def test(tracker: dict[str, entry], res: list[list[entry]]) -> None:
    """helper function to print the orthogroup meeting the conditions"""

    global total

    sppKey = ["Cmyd.def.def", "Mter.def.def", "Tscr.def.def", "Agig.def.def", "Gfla.def.def"]

    c = tracker[sppKey[0]]
    m = tracker[sppKey[1]]
    t = tracker[sppKey[2]]
    a = tracker[sppKey[3]]
    g = tracker[sppKey[4]]

    # options
    A  = ['A']
    C  = ['C']
    E  = ['E']
    AC = ['A', 'C']
    AE = ['A', 'E']

    if (c.type in A and m.type in E and t.type in AC and a.type in AC and g.type in AC):
        if (a.count > g.count):
            return
        spp = [c, m, t, a, g]
        spp.sort(key=lambda x: x.count, reverse=True)
        res.append(spp)
        # for s in spp:
        #     print(s)
        # print()

        total += 1

def parse_orthologs(ofile: str) -> None:
    """find the orthogroups of interst"""
    
    cur     = ''
    fh      = open(ofile, 'r')
    tracker = dict()
    res     = list()

    for line in fh:
        if (line[0] == '#'):
            continue
        fields = line.split('\t')
        if (fields[7] == 'D' or fields[3] != "5"):
            continue
        if (cur == ''):
            cur     = fields[0]
            tracker = dict() 
        if (cur == fields[0]):
            tracker[fields[4]] = entry(fields)
        else:
            test(tracker, res)
            cur                = fields[0]
            tracker            = dict()
            tracker[fields[4]] = entry(fields)
        
    fh.close()

    test(tracker, res)

    global total
    print("Total orthogroups:", total)

    res.sort(key=lambda x: x[0].count, reverse=True)

    for enList in res:
        for en in enList:
            print(en)
        print()

def main() -> int:
    """entry point to this little helper code"""

    # get the orthologs.tsv file
    ofile = get_arguments()

    # parse
    parse_orthologs(ofile)

    return 0

if __name__ == "__main__":
    main()
