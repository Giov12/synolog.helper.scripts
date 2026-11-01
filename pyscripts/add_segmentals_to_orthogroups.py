#!/bin/env python3

import argparse
import os
import glob
import sys
from  collections import defaultdict

class Candidate:

    #
    # class to hold the gene duplication information
    # needed to record it in orthologs.tsv
    # 
    
    __slots__ = ("Name", "Chr", "StartBp", "EndBp", "syntenic_pairs")
    
    def __init__(self, name: str, chrom: str, start_bp: str, end_bp: str):
        self.syntenic_pairs = list()
        self.StartBp        = start_bp
        self.EndBp          = end_bp
        self.Name           = name
        self.Chr            = chrom

    def spp_support(self) -> int:
        """count the number of species supporting this duplication"""

        global orthologs

        # in practice, ~10 species will be present.
        # fast to iteratve over
        species = list()

        for _, other_spp in self.syntenic_pairs:
            if (other_spp in orthologs):
                spp = other_spp.split('_', 1)[0]
                if (spp not in species):
                    species.append(spp)
        return len(species)


# globals
orthologs     = {}               # str -> str [str will be a string int | group ID]
orthogroups   = defaultdict(set) # str -> set[str]
candidates    = dict()           # str -> Candidate object
selected      = {}               # str -> str [str will be a string int | group ID]
path          = ''
retroset      = set()
selected_grps = defaultdict(list)
min_spp       = 2                # num of species supporting the segmental duplication

def get_arguments() -> int:
    """get the arguments"""

    global path, min_spp

    desc  = "Add unplaced segmental duplications to an orthogroup"
    phelp = "Path to the directory containing the orthogroups.tsv and *segmentalduplications.tsv files from synolog"
    mhelp = "Minimum number of other species to support the segmental duplication"

    parser = argparse.ArgumentParser(description = desc)
    parser.add_argument("-p", "--path",    help=phelp, required=True)
    parser.add_argument("-m", "--min-spp", help=mhelp, default=min_spp, type=int)
    args = parser.parse_args()
    path = args.path

    assert os.path.isdir(path), f"Could not locate directory {path}"

    if (path[-1] == '/'):
        path = path[:-1] # remove path char

    assert args.min_spp >= 0, "--min-spp cannot be less than 0"
    min_spp = args.min_spp

    return 0

def get_segmental_duplications() -> list[str]:
    """get the segmental duplication files in the provided path"""

    global path

    suffix   = "_segmentalduplications.tsv"
    segfiles = glob.glob(f"{path}/*{suffix}")

    assert len(segfiles) > 0, f"Could not locate any *{suffix} files in {path}/"

    return segfiles

def load_orthologs() -> int:
    """load the current members into memory"""

    global path, orthologs, orthogroups, retroset

    infile = f"{path}/orthologs.tsv"

    assert os.path.isfile(infile), f"Could not locate orthologs.tsv in {path}/"

    fh = open(infile, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#' or line[0] == '\n'):
            continue
        fields = line.split('\t')
        assert len(fields) == 12, f"Error: Expected 12 columns in orthologs.tsv. Offending line:\n{line}"
        ogrp     = fields[0]
        spp      = fields[4]
        gene_id  = fields[5]
        eState   = fields[7]
        member   = f"{spp}_{gene_id}"
        if (eState != 'R'):
            orthogroups[ogrp].add(member)
            orthologs[member] = ogrp
        else:
            retroset.add(member)

    fh.close()

    print(f"Loaded {len(orthologs)} orthologs across {len(orthogroups)} orthogroups")
    print(f"Loaded {len(retroset)} genes flagged as retrogenes")
    
    return 0

def load_candidates() -> int:
    """find segmental duplications not in orthogroups.tsv"""

    global orthologs, retroset, candidates

    segfiles = get_segmental_duplications()

    for segfile in segfiles:
        fh = open(segfile, 'r')

        for line in fh:
            if (len(line) == 0 or line[0] == '#'):
                continue
            fields = line.split('\t')
            assert len(fields) == 22, \
                f"Error: Expected 22 columns in {os.path.basename(segfile)}. Offending line:\n{line}"
            spp    = fields[2]
            dup_id = fields[3]
            cand   = f"{spp}_{dup_id}"
            if (cand in orthologs or cand in retroset):
                continue # was placed into orthologs.tsv
            gene_id = f"{spp}_{fields[9]}"         # duplicate of this gene
            other   = f"{fields[15]}_{fields[16]}" # species B
            if (cand not in candidates):
                candidates[cand] = Candidate(fields[4], fields[5], fields[7], fields[8])
            candidates[cand].syntenic_pairs.append((gene_id, other))

        fh.close()

    print(f"Loaded {len(candidates)} candidate segmental duplications")
    return 0

def select_candidates() -> int:
    """helper function to go through the candidates and see if they can be added into the orthogroup"""

    global orthologs, candidates, selected, selected_grps, min_spp

    ambigous = 0
    low_supp = 0

    for mem, cand in candidates.items():
        syn_pairs = cand.syntenic_pairs
        ogroups   = set()
        for syn_pair in syn_pairs:
            same_spp      = syn_pair[0]
            other_spp     = syn_pair[1]
            if (same_spp in orthologs):
                ogroups.add(orthologs[same_spp])
            if (other_spp in orthologs):
                ogroups.add(orthologs[other_spp])
        if (len(ogroups) != 1):
            ambigous += 1
        elif (cand.spp_support() < min_spp):
            low_supp += 1
        else:
            grp           = next(iter(ogroups))
            selected[mem] = grp
            selected_grps[grp].append(mem)

    if (len(selected) == 0):
        msg = "No candidates past any of the filters. Nothing to write."
        sys.exit(msg)
    print(f"{ambigous} segmental duplications were mapped to multiple orthogroups and ignored")
    print(f"{low_supp} segmental duplications did not meet the minimum number of species support required")
    print(f"Found {len(selected)} candidate segmental duplications to add to orthogroups.tsv")

    return 0

def write_candidates() -> int:
    """helper function to create a new orthologs.tsv"""

    global path, selected_grps, selected, candidates

    infile  = f"{path}/orthologs.tsv"
    outfile = f"{path}/orthologs.with.segdups.tsv"
    fh      = open(infile, 'r')
    ofh     = open(outfile, 'w')
    cur_grp = defaultdict(list)
    prev    = ''
    tot     = 0

    #
    # header
    # [GroupID] [Type] [ElementCount] [OrgCount] [Org] [ElementId] [Name]
    # [CopyNumState] [ParalogCount] [Chr] [StartBp] [EndBp]
    #

    for line in fh:
        if (len(line) == 0 or line[0] == '\n'):
            continue
        if (line[0] == '#'): # header
            ofh.write(line)
            continue
        fields = line.split('\t')
        grp    = fields[0]
        if (prev == ''):
            prev = grp
        if (prev != grp):
            if (prev in selected_grps):
                # get the segmentdal duplications for this group
                cands = selected_grps[prev]
                tot  += len(cands)
                for mem in cands:
                    cand         = candidates[mem]
                    spp, elem_id = mem.split('_', 1)
                    newline      = '\t'.join([prev, 'gene', str(tot), str(len(cur_grp)), spp, elem_id, 
                                              cand.Name, 'S', '0', cand.Chr, cand.StartBp, cand.EndBp + '\n'])
                    cur_grp[spp].append(newline)
            for spp, newlines in cur_grp.items():
                for newline in newlines:
                    line_fields    = newline.split('\t')
                    line_fields[2] = str(tot)
                    line_fields[3] = str(len(cur_grp))
                    ofh.write('\t'.join(line_fields))
            cur_grp.clear()
            prev = grp
            tot  = 0
        spp = fields[4]
        cur_grp[spp].append(line)
        tot += 1

    fh.close()

    # write the last group
    if (prev in selected_grps):
        cands = selected_grps[prev]
        tot  += len(cands)
        for mem in cands:
            cand         = candidates[mem]
            spp, elem_id = mem.split('_', 1)
            newline      = '\t'.join([prev, 'gene', str(tot), str(len(cur_grp)), spp, elem_id, 
                                      cand.Name, 'S', '0', cand.Chr, cand.StartBp, cand.EndBp + '\n'])
            cur_grp[spp].append(newline)
    for spp, newlines in cur_grp.items():
        for newline in newlines:
            line_fields    = newline.split('\t')
            line_fields[2] = str(tot)
            line_fields[3] = str(len(cur_grp))
            ofh.write('\t'.join(line_fields))

    ofh.close()

    return 0

def main() -> int:
    """entry point to this helper script"""

    # get the input path
    get_arguments()

    # load the core members
    load_orthologs()

    # filter the segmental duplications for candidate genes
    load_candidates()

    # now to parse through and select candidates that pass the screening
    select_candidates()

    # write a new orthologs.tsv file
    write_candidates()

    return 0

if __name__ == "__main__":
    main()
