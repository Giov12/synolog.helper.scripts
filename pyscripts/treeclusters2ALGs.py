#!/bin/env python3

import argparse
import os
from  collections import defaultdict

clusters = defaultdict(list)

class Gene:

    __slots__ = ("id", "name")
    def __init__(self, gid: str, name: str) -> None:
        self.id   = gid
        self.name = name

class Cluster:
    def __init__(self, id_: int) -> None:
        self.id      = id_
        self.members = defaultdict(list)

    def add_member(self, org: str, genes: list[Gene]) -> None:
        self.members[org].extend(genes)

    def empty(self) -> bool:
        return len(self.members) == 0

def get_arguments() -> str:
    """get the single argument"""

    d = "script to take a tree_clusters_membership.tsv file and count the number of ALGs"

    parser  = argparse.ArgumentParser(description = d)
    parser.add_argument("-t", "--tree-clusters", help="tree_clusters_membership.tsv file", required=True)
    args    = parser.parse_args()
    tcfile  = args.tree_clusters

    assert os.path.isfile(tcfile), f"Could not locate directory {tcfile}"

    return tcfile

def parse_clusters(tcfile: str) -> int:
    """parse the tree clusters for the algs"""

    global clusters

    fh    = open(tcfile, 'r')
    N     = -1
    cur   = "-1"
    clstr = defaultdict(list)
    mems  = set()
    total = 0

    for line in fh:
        if (len(line) == 0):
            continue
        elif (line[0] == '#'):
            if (line.startswith("#Orgs")):
                j = -1
                for i in range(len(line)):
                    if (line[i] == ' '):
                        j = i + 1
                    elif (line[i] == '\n'):
                        N = int(line[j:i])
            continue
        fields = line.split('\t')
        count  = int(fields[2])
        if (count != N):
            if (len(mems) > 0):
                keys   = sorted(list(mems))
                keys   = '@'.join(keys)
                total += 1
                nclstr = Cluster(total)
                for mem, genes in clstr.items():
                    nclstr.add_member(mem, genes)
                clusters[keys].append(nclstr)
                mems  = set()
                clstr = defaultdict(list)
            continue
        clID  = fields[0]
        org   = fields[3]
        chrom = fields[4]
        key   = org + ':' + chrom
        gid   = fields[5]
        name  = fields[6]
        gene  = Gene(gid, name)

        if (cur == "-1"):
            cur = clID
        if (cur == clID):
            mems.add(key)
            clstr[org].append(gene)
        else:
            if (len(mems) > 0):
                keys = sorted(list(mems))
                keys = '@'.join(keys)
                total += 1
                nclstr = Cluster(total)
                for mem, genes in clstr.items():
                    nclstr.add_member(mem, genes)
                clusters[keys].append(nclstr)
            mems  = set()
            clstr = defaultdict(list)
            mems.add(key)
            clstr[org].append(gene)
            cur = clID        

    fh.close()

    if (len(mems) > 0):
        keys = sorted(list(mems))
        keys = '@'.join(keys)
        total += 1
        nclstr = Cluster(total)
        for mem, genes in clstr.items():
            nclstr.add_member(mem, genes)
        clusters[keys].append(nclstr)

    return 0

def make_algs() -> int:
    """take the clusters and quantify the members"""

    global clusters

    print("Number of ALGs:", len(clusters))
    fh = open("ALGs.infered.txt", 'w')

    for key, clstrs in clusters.items():
        mems = key.split('@')
        if (len(mems) < 5):
            continue
        entry = dict()
        for mem in mems:
            org, chrom = mem.split(':')
            entry[org] = [chrom]
        members = defaultdict(list)
        total   = 0
        for clstr in clstrs:
            for org, genes in clstr.members.items():
                members[org].extend(genes)
        for org, genes in members.items():
            entry[org].append(str(len(genes)))
            total += len(genes)
        out = list()
        for org, entries in entry.items():
            chrom = entries[0]
            count = entries[1]
            rec   = f"{org} ({chrom}, {count})"
            out.append(rec)
        outline = ", ".join(out) + ' ' + str(total) + '\n'
        fh.write(outline)
        
    fh.close()

    return 0

def main() -> int:
    """entry point to this little subprogram"""

    # get the input file
    tcfile = get_arguments()

    # parse the clusters
    parse_clusters(tcfile)

    # now make the clusters
    make_algs()

    return 0

if __name__ == "__main__":
    main()
