#!/bin/env python3

import argparse
import os
import gzip

"""
    todos:
        is subset
        is split in othergroups
        is not in any (check other softwares)
        is equal
    
        compare to other softwares

"""

class OrthoRes:
    def __init__(self, methodName: str) -> None:
        self.name      = methodName
        self.geneMap   = dict()
        self.orthogrps = list()

    def get_name(self) -> str:
        return self.name
    
    def get_geneMap(self) -> dict:
        return self.geneMap
    
    def get_orthogrps(self) -> list:
        return self.orthogrps
    
    def parse_reformat_tsv(self, tsv: str) -> None:
        """fill the geneMap dict and create a list of the ortholog groups"""

        sppList = list()
        fh      = gzip.open(tsv, "rt") if tsv.endswith(".gz") else open(tsv, 'r')

        for line in fh:
            fields = line.strip().split('\t')
            if (line[0] == '#'):
                sppList = fields[1:]
                continue
            orthogroup = set()
            idx        = len(self.orthogrps)
            for i in range(1, len(fields)):
                if (fields[i] == ''): continue
                spp      = sppList[i - 1]
                subfield = fields[i].split(", ")
                for gene in subfield:
                    gene = f"{spp}:{gene}"
                    orthogroup.add(gene)
                    self.geneMap[gene] = idx
            
            self.orthogrps.append(orthogroup)
        fh.close()

    def load_synolog(self, osyn: str) -> None:
        """"this function does the actual reading and loading"""

        fh = gzip.open(osyn, "rt") if osyn.endswith(".gz") else open(osyn, 'r')

        prev = None
        og   = set()
        idx  = 0 # ortholog group id for class
        for line in fh:
            if (line[0] == '#'):
                continue
            fields = line.split('\t')
            gi     = fields[0] # group id
            spp    = fields[4].split('.')[0]
            ge     = fields[5]
            gene   = f"{spp}:{ge}"
            if (prev == None):
                prev = gi
            if (gi == prev):
                og.add(gene)
                self.geneMap[gene] = idx
            else:
                self.orthogrps.append(og)
                og   = set()
                idx += 1
                prev = gi
                og.add(gene)
                self.geneMap[gene] = idx

        
        # add last group
        self.orthogrps.append(og)

        fh.close()
        
def get_arguments() -> tuple:
    """get the arguments"""

    d = "Some python code to unify the results of OrthoFinder, InParanoid, & OrthoMCL"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-F", "--orthofinder", help="reformated tsv of orthofinder",       required=True)
    parser.add_argument("-I", "--inparanoid",  help="table of merged inparnoid tables",    required=True)
    parser.add_argument("-M", "--orthoMCL",    help="reformated tsv of OrthoMCL pipeline", required=True)
    parser.add_argument("-S", "--synolog",     help="orthologs.tsv file from Synolog",     required=True)

    args = parser.parse_args()
    ofdr = args.orthofinder
    iprd = args.inparanoid
    omcl = args.orthoMCL
    osyn = args.synolog

    assert os.path.isfile(ofdr), f"Could not locate file: {ofdr}"
    assert os.path.isfile(omcl), f"Could not locate file: {omcl}"
    assert os.path.isfile(iprd), f"Could not locate file: {iprd}"
    assert os.path.isfile(osyn), f"Could not locate file: {osyn}"
    
    return (ofdr, iprd, omcl, osyn)

def parse_methods(ofdr: str, iprd: str, omcl: str) -> list:
    """return a list of each methods results in the OrthoRes data structure"""

    res     = []
    mkey    = {0 : "OrthoFinder", 1 : "Inparanoid", 2: "OrthoMCL"}
    methods = [ofdr, iprd, omcl]

    for i, m in enumerate(methods):
        or_obj = OrthoRes(mkey[i])
        or_obj.parse_reformat_tsv(m)
        res.append(or_obj)

    return res

def compare_methods(res: list) -> None:
    """compare the ortholo groups across the three different methods"""

    visited = set() # avoid repeated work

    #
    # orthofinder = 0, inparanoid = 1, orthomcl = 2
    #
    
    newOrthoGroups = list()

    for m, or_obj in enumerate(res):
        gmap   = or_obj.get_geneMap()
        ogrps1 = or_obj.get_orthogrps()
        for gene, i in gmap.items():
            if (gene in visited):
                continue
            og1 = ogrps1[i]
            visited.update(og1) # add all genes from this orthogrp to visited
            nog = og1.copy() # will be pruned to create a new orthogroup
            for j, or_obj2 in enumerate(res):
                if (j == m): continue
                ogrps2 = or_obj2.get_orthogrps()
                gmap2  = or_obj2.get_geneMap()
                for g in og1:
                    if (g not in gmap2):
                        nog.discard(g)
                        continue
                    k   = gmap2[g]
                    og2 = ogrps2[k]
                    visited.update(og2)
                    nog = nog.intersection(og2)
                    break

            if (len(nog) > 1):
                newOrthoGroups.append(nog)
    print("Number of new orthogroups:", len(newOrthoGroups))
    return newOrthoGroups

def write_merged(mlist: list, omcl: str) -> None:
    """write the merged list to a new tsv file"""

    if (os.path.isfile("Merged.OrthoGroups.tsv")):
        return

    # take the header line from one of the input files
    ofh = open("Merged.OrthoGroups.tsv", 'w')

    fh     = gzip.open(omcl, 'r') if omcl.endswith(".gz") else open(omcl, 'r')
    header = fh.readline()
    fh.close()
    ofh.write(header)

    sppList = header.strip().split('\t')[1:]

    for i, ogrp in enumerate(mlist):
        spp_dict = {s : [] for s in sppList}
        row      = [f"{i}"]
        for gene in ogrp:
            field = gene.split(':')
            spp   = field[0]
            ge    = field[1]
            spp_dict[spp].append(ge)
        for spp in sppList:
            row.append(", ".join(spp_dict[spp]) if spp_dict[spp] else '')
        line = '\t'.join(row) + '\n'
        ofh.write(line)
 
    ofh.close()

def read_synolog(osyn: str) -> OrthoRes:
    """load in the orthologs.tsv file"""

    or_obj = OrthoRes("synolog")

    or_obj.load_synolog(osyn)

    return or_obj

def create_merge_map(mlist: list) -> dict:
    """create a gene map of the merged orthologs"""

    geneMap = {}

    for i, ogrp in enumerate(mlist):
        for g in ogrp:  
            geneMap[g] = i
    return geneMap

def filter_by_missing(miss: set, gmap: dict, revisit: set) -> None:
    """filter the missing genes out that have their orthogroup already revisit"""

    for g in miss:
        idx = gmap[g]
        if (idx not in revisit):
            revisit.add(idx)


def compare_syn_with_merged(mlist: list, orSyn: OrthoRes) -> set:
    """"compare the synolog results to the merged coherent groups"""

    ogrps = orSyn.get_orthogrps()
    gmap  = orSyn.get_geneMap()
    mmap  = create_merge_map(mlist)

    # create remainers file?
    splitcnt = 0
    notinmer = 0 # not in merged
    eqlCnter = set()
    subcnter = set()
    supcnter = set()
    missCnt  = set()
    visited  = set()
    revisit  = set()

    # collection of file handlers
    fh1 = open("Not.in.merged.txt", 'w')
    fh2 = open("Subset.txt", 'w')
    fh3 = open("Superset.txt", 'w')
    fh4 = open("Equal.txt", 'w')
    fh5 = open("Diff.txt", 'w')
    fh6 = open("Split.txt", 'w')

    for gene, i in gmap.items():
        if (i in visited): continue
        visited.add(i)
        ogrp = ogrps[i]
        oids = set() # other ids
        subsetTrue = False
        supsetTrue = False
        equalTrue  = False
        for g in ogrp:
            if (g not in mmap):
                fh1.write(g + '\n')
                notinmer += 1
                revisit.add(i)
                missCnt.add(i)
                continue
            ### continued logic
            j     = mmap[g]
            ogrp2 = mlist[j]
            if (len(ogrp2) > len(ogrp)):
                if (ogrp.issubset(ogrp2)):
                    subcnter.add(i)
                    fh2.write(g + '\n')
                else:
                    revisit.add(i)
                # check subset approach
            elif (len(ogrp2) < len(ogrp)):
                if (ogrp2.issubset(ogrp)):
                    supcnter.add(i)
                    fh3.write(g + '\n')
                else:
                    revisit.add(i)
                # check superset & get left overs
            elif (ogrp2 == ogrp):
                eqlCnter.add(i)
                fh4.write(g + '\n')
            else:
                revisit.add(i)
                fh5.write(g + '\n')
            oids.add(j)
                # check for left overs / remainders
            # should be call cases
        if (len(oids) > 1):
            revisit.add(i)
            splitcnt += 1
            for g in ogrp:
                fh6.write(g + '\n')
    
    print("Equal count:", len(eqlCnter))
    print("Subset count:", len(subcnter))
    print("Superset count:", len(supcnter))
    print("Number of splits:", splitcnt)
    print("Number not in merged:", notinmer, "and", len(missCnt))
    print("Number of orthogroups to revisit:", len(revisit))

    for f in [fh1, fh2, fh3, fh4, fh5, fh6]:
        f.close()

    return revisit

def write_missing(miss: set) -> None:
    """write out the list of genes that no method was able to find"""

    misSorted = sorted(miss)
    fh        = open("MissingGenes.txt", 'w')

    for g in misSorted:
        fh.write(f"{g}\n")

    fh.close()

def check_other_software(res: list, orSyn : OrthoRes, revisit: set) -> None:
    """compare the orthogroups that did not match the merged results to all other software"""

    rvts  = sorted(revisit)
    ogrps = orSyn.get_orthogrps()
    miss  = set()
    ct1 = ct2 = ct3 = 0

    print("Total count:", len(ogrps))
    print("Before:", len(revisit))
    for i in rvts:
        ogrp = ogrps[i]
        fnd  = False # if found in one of the prev methods
        for or_obj in res:
            gmap2  = or_obj.get_geneMap()
            ogrps2 = or_obj.get_orthogrps()
            for g in ogrp:
                if (g not in gmap2):
                    miss.add(g)
                    continue
                miss.discard(g) # if another method ended up finding it
                j     = gmap2[g]
                ogrp2 = ogrps2[j]
                if (len(ogrps2) > len(ogrp)):
                    if (ogrp.issubset(ogrp2)):
                        revisit.discard(i)
                        fnd = True
                        ct1 += 1
                        break
                elif (len(ogrp2) < len(ogrp)):
                    if (ogrp2.issubset(ogrp)):
                        revisit.discard(i)
                        fnd = True
                        ct2 += 1
                        break
                elif (ogrp2 == ogrp):
                    revisit.discard(i)
                    fnd = True
                    ct3 += 1
                    break
            if fnd:
                break

    print("After:", len(revisit))
    print("Missing count:", len(miss))

    print("Is subset:", ct1)
    print("Is superset:", ct2)
    print("Is equal:", ct3)

    # gmap = orSyn.get_geneMap()
    # filter_by_missing(miss, gmap, revisit)
    write_missing(miss)

def compare_to_all_methods(res: list, orSyn: OrthoRes) -> None:
    """Just print out how much is the same with the other three methods"""

    ogrps    = orSyn.get_orthogrps()
    syngenes = set(orSyn.get_geneMap().keys())

    for or_obj in res:
        subcnt  = 0
        supcnt  = 0
        equal   = 0
        split   = 0
        gmap    = or_obj.get_geneMap()
        ogrp2Gs = set(gmap.keys())
        ogrps2  = or_obj.get_orthogrps()
        for ogrp in ogrps:
            oids = set()
            for g in ogrp:
                if (g not in gmap):
                    continue
                j     = gmap[g] # index of other group
                ogrp2 = ogrps2[j]
                if (len(ogrp2) > len(ogrp)):
                    if (ogrp.issubset(ogrp2)):
                        subcnt += 1
                        break
                elif (len(ogrp2) < len(ogrp)):
                    if (ogrp2.issubset(ogrp)):
                        supcnt += 1
                        break
                elif (ogrp == ogrp2):
                    equal += 1
                    break
                oids.add(j)
            split += 1 if (len(oids) > 1) else 0
        # 
        # 
        nfnd = len(ogrp2Gs.difference(syngenes))
        miss = len(syngenes.difference(ogrp2Gs))
        
        # 
        # 
        print(f"For method {or_obj.name}:")
        print("\t\tMissing:", miss)
        print("\t\tNotFound:", nfnd)
        print("\t\tEqual:", equal)
        print("\t\tSubset:", subcnt)
        print("\t\tSupercnt:", supcnt)
        print("\t\tSplit:", split)


def write_revisit(revisit: set, orSyn: OrthoRes) -> None:
    """write the orthogroups that need to be revisited"""

    orgrps = orSyn.get_orthogrps()
    rvist  = sorted(revisit)
    fh     = open("RevisitList.txt", 'w')

    tot    = 0

    for i in rvist:
        ogrp  = orgrps[i]
        tot  += len(ogrp)
        ogrp = sorted(ogrp) 
        genes = ", ".join(ogrp)
        fh.write(genes + '\n')

    fh.close()

    print("Total genes to revisit:", tot)

def main() -> int:
    """Parse and merge results across ortholog-calling software"""

    # get arguments
    ofdr, iprd, omcl, osyn = get_arguments()

    # load the results
    res = parse_methods(ofdr, iprd, omcl)

    # create the merged list
    mlist = compare_methods(res)

    # save the merged table
    write_merged(mlist, omcl)

    # read in synolog orthologs
    orSyn = read_synolog(osyn)

    # compare the merged to synolog
    revisit = compare_syn_with_merged(mlist, orSyn)

    # see if we can resolve orthogroups that do not match the other groups
    check_other_software(res, orSyn, revisit)

    write_revisit(revisit, orSyn)

    # now a general comparison to all three softwares
    compare_to_all_methods(res, orSyn)

    return 0

if __name__ == "__main__":
    main()
