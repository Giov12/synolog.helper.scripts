#!/bin/env python3

import argparse
import os
import sys
from collections import defaultdict
from glob import glob

path          = ''
assignments   = dict() # tuple (org, gene_id) -> group ID
consyn_assign = defaultdict(dict)

def set_arguments() -> int:
    """get & set the arguments"""

    global path

    d = "Check if genes marked as conserved syntenic pairs are within the same orthogroup"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-p", "--path", help="Path to results generated from synolog", required=True, type=str)

    args = parser.parse_args()
    assert os.path.isdir(args.path), f"Could not verify that {path} is a directory"
    path = args.path

    if (path[-1] == '/'):
        path = path[:-1]
    
    return 0

def parse_orthogroups() -> int:
    """parse the orthogroups file and populate the orthogroups & assignments data structs"""

    global path, assignments

    ofile = f"{path}/orthologs.tsv"
    assert os.path.isfile(ofile), f"Could not locate orthologs.tsv in {path}/"
    fh    = open(ofile, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.split('\t')
        group  = fields[0]
        org    = fields[4]
        gid    = fields[5]
        entry  = (org, gid)
        assignments[entry] = group

    fh.close()

    return 0

def get_consyn_files() -> list[str]:
    """get the conserved synteny block files from the path"""

    global path

    consyn_files = glob(f"{path}/*_consyn_cluster_membership.tsv")

    assert len(consyn_files) > 0, f"Could not locate any conserved synteny files at {path}/"

    return consyn_files

def report_contradictions() -> int:
    """write out the number of contradicting results"""

    global assignments, consyn_assign

    ofh = open("Contradictions.tsv", 'w')
    cnt = 0
    mis = 0

    consyn_files = get_consyn_files()

    for cfile in consyn_files:
        cCnt  = 0
        mCnt  = 0
        fh    = open(cfile, 'r')
        bname = os.path.basename(cfile)
        for line in fh:
            if (len(line) == 0 or line[0] == '#'):
                continue
            fields = line.split('\t')
            assert len(fields) == 16, f"Unexpected number of columns found in {bname}"
            orgA   = fields[2]
            orgB   = fields[3]
            gidA   = fields[4]
            gidB   = fields[10]
            entryA = (orgA, gidA)
            entryB = (orgB, gidB)
            consyn_assign[entryA][orgB] = gidB
            consyn_assign[entryB][orgA] = gidA
            if (entryA not in assignments or entryB not in assignments):
                mCnt += 1
                if (entryA not in assignments):
                    msg = f"Missing {orgA} {gidA} fron {bname}\n"
                    ofh.write(msg)
                if (entryB not in assignments):
                    msg = f"Missing {orgB} {gidB} fron {bname}\n"
                    ofh.write(msg)
            else:
                grpA = assignments[entryA]
                grpB = assignments[entryB]
                if (grpA != grpB):
                    cCnt += 1
                    msg = f"Members {orgA} {gidA} and {orgB} {gidB} assigned to different orthogroups in {bname}\n"
                    ofh.write(msg)

        fh.close()
        print(f"{bname} contradictions: {cCnt}, missing {mCnt}")
        cnt += cCnt
        mis += mCnt

    ofh.close()

    print(f"Total number of contradictions: {cnt}; Total number of missing: {mis}")

    return 0

def parse_tree_members(entries: str) -> list[tuple[str, str]]:
    """parse the last columns in a tree clusters file"""

    members = list()
    fields  = entries.split(", ")

    for field in fields:
        subfields = field.split(' ')
        org       = subfields[0].strip()
        gid       = ''
        j         = 0
        for i in range(len(subfields[1])):
            if (subfields[1][i] == '('):
                j = i + 1
            elif (subfields[1][i] == ')'):
                gid   = subfields[1][j:i]
                entry = (org, gid)
                members.append(entry)
                break

    return members

def check_tree_clusters() -> int:
    """helper function to verify that the tree clusters match what is reported in the conserved synteny files"""

    global path, consyn_assign

    tfile = f"{path}/TreeClustersMembership.tsv"

    if (os.path.isfile(tfile) == False):
        msg = f"Unable to find TreeClustersMembership.tsv at {path}/"
        sys.exit(msg)

    fh   = open(tfile, 'r')
    ofh  = open("TreeClusterContradictions.tsv", 'w')
    mcnt = 0 # missing
    dcnt = 0 # different

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.split('\t')
        org    = fields[3]
        gid    = fields[5]
        entry  = (org, gid)
        if (entry not in consyn_assign):
            mcnt += 1
            msg = f"Missing {org} {gid} in consyn files\n"
            ofh.write(msg)
            continue
        diff = list()
        mems = parse_tree_members(fields[-1])
        for mem in mems:
            mem_org = mem[0]
            mem_gid = mem[1]
            if (mem_org not in consyn_assign[entry]):
                diff.append(mem)
                continue
            assigned = consyn_assign[entry][mem_org]
            if (assigned != mem_gid):
                diff.append(mem)
        if (len(diff) > 0):
            dcnt += 1
            out   = list()
            for mem in diff:
                dorg  = mem[0]
                dgid  = mem[1]
                entry = f"({dorg} {dgid})"
                out.append(entry)
            out = ','.join(out)
            msg = f"Issues with {org} {gid} and {out}\n"
            ofh.write(msg)

    fh.close()
    ofh.close()

    print(f"Total number of missing syntenic genes on tree: {mcnt}; Total number of contradicting genes {dcnt}")

    return 0

def main() -> int:
    """Get the path to the synolog output files"""

    # get the path
    set_arguments()

    # get the orthogroups
    parse_orthogroups()

    # check for contradictions
    report_contradictions()

    # now do tree clusters
    check_tree_clusters()

    return 0

if __name__ == "__main__":
    main()
