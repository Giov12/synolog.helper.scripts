#!/bin/env python3

import argparse, cairo, math, os, sys, gzip
from collections import defaultdict
import numpy as np


def get_arguments():
    """
        Function to get membership files for Gene Anchors, 
        GFF TE locations, and chromosome.tsv from Synolog 
        Output
    """

    arg_description = "A program to detect putative orthologous TEs across Species by using neighboring orthologous genes"

    parser = argparse.ArgumentParser(description=arg_description)

    parser.add_argument("-m", "--memberships_tsv", help="Membership file between the two species", type=str, required=True)
    parser.add_argument("-g1", "--species1_gff", help="GFF file of the first species containing TE locations",type=str, required=True)
    parser.add_argument("-g2", "--species2_gff", help="GFF file of the second species containing TE locations",type=str, required=True)
    parser.add_argument("-c", "--chromosomes_tsv", help="chromosomes tsv file generated from Synolog analysis",type=str, required=True)
    parser.add_argument("-M", "--min_len", help = "Minimum length to retain a scaffold [default = 1e6]", type = int, default = 1e6)
    parser.add_argument("-d", "--max_distance", help = "Maximum distance from anchor to consider a TE [default = 1e6]", type = int, default = 1e6)

    args = parser.parse_args()
    mem_tsv = args.memberships_tsv
    spA_gff = args.species1_gff
    spB_gff = args.species2_gff
    chroms_tsv = args.chromosomes_tsv
    min_len = args.min_len
    max_dist = args.max_distance

    # check if they all exist
    for f in [mem_tsv, spA_gff, spB_gff, chroms_tsv]:
        assert os.path.isfile(f), f"Could not locate {f}"

    # assert limits/thresholds
    assert min_len >= 0, f"{min_len} is not a legal minimum length"
    assert max_dist >= 0, f"{max_dist} is not a legal maximum TE count"

    return mem_tsv, spA_gff, spB_gff, chroms_tsv, min_len, max_dist


class TE:
    def __init__(self, family: str, subfamily: str, start: int, end: int, orientation: str):
        self.family = family
        self.subfamily = subfamily
        self.start = start
        self.end = end
        self.orientation = orientation
        self.assigned = False # set true when/if assigned
        self.ortho_indices = [] # indices of the (co-) orthologs

    def assign(self, ortho_index: int):
        self.assigned = True # constant update
        self.ortho_indices.append(ortho_index)

# loop through and make each row an object to store in a list
class Ortho_Pair:
    def __init__(self, spA: str, spA_chrom: str, spB: str, spB_chrom: str, 
                 startA: int, endA: int, startB: int, endB: int, geneID_A: str,
                 geneID_B: str, cluster_id: str):
        self.spA = spA
        self.spB = spB
        self.spA_chrom = spA_chrom
        self.spB_chrom = spB_chrom
        self.start_posA = startA
        self.end_posA = endA
        self.start_posB = startB
        self.end_posB = endB
        self.geneID_A = geneID_A
        self.geneID_B = geneID_B
        self.cluster_id = cluster_id

def get_species_ids(mem_tsv: str):
    """parse membership file name to get species IDs"""
    
    bname = os.path.basename(mem_tsv)

    # will remove this extension
    ex = "_consyn_cluster_membership.tsv"
    if bname.endswith(".gz"):
        ex = ex + ".gz"
    
    bname = bname.replace(ex, '')
    
    spp = bname.split('-')
    
    return spp[0], spp[1]

def parse_chromosomes_tsv(chroms_tsv: str, mem_tsv: str, min_len: int):
    """parse the chromosomes file for chrom lengths"""

    # get species IDs from memberships file

    fh = gzip.open(chroms_tsv, "rt") if chroms_tsv.endswith(".gz") else open(chroms_tsv, 'r')
    spA_lengths = {}
    spB_lengths = {}

    # start with first species
    spA, spB = get_species_ids(mem_tsv)

    print(f"Using {spA} as the first species\nUsing {spB} as the second species")

    for line in fh:
        # skip comments
        if line[0] == '#':
            continue
        fields = line.strip().split('\t')
        if fields[0] == spA and int(fields[2]) >= min_len:
            spA_lengths[fields[1]] = int(fields[2])
        elif fields[0] == spB and int(fields[2]) >= min_len:
            spB_lengths[fields[1]] = int(fields[2])

    fh.close()

    return spA_lengths, spB_lengths, spA, spB


def parse_memberships_tsv(mem_tsv: str, spA: str, spB: str):
    """parse the membership files to get gene orthologous pairs"""

    fh = gzip.open(mem_tsv, "rt") if mem_tsv.endswith(".gz") else open(mem_tsv, 'r')

    # create a list to store all the Single_Ortho_Row objects
    ortho_pairs = []

    for line in fh:
        if line[0] == '#':
            continue
        fields = line.strip().split('\t')
        # get the components of each row
        cluster_id = fields[0]
        spA_id = fields[1]
        spB_id = fields[2]
        geneID_A = fields[3]
        chrom_A = fields[5]
        start_A = int(fields[7])
        end_A = int(fields[8])
        geneID_B = fields[9]
        chrom_B = fields[11]
        start_B = int(fields[13])
        end_B = int(fields[14])
        # do we need to switch?
        if spA_id == spB or spB_id == spA:
            spA_id, spB_id = spB_id, spA_id
            geneID_A, geneID_B = geneID_B, geneID_A
            chrom_A, chrom_B = chrom_B, chrom_A
            start_A, start_B = start_B, start_A
            end_A, end_B = end_B, end_A
        ortho_pair = Ortho_Pair(spA_id, chrom_A, spB_id, chrom_B,
                                start_A, end_A, start_B, end_B, 
                                geneID_A, geneID_B, cluster_id)
        ortho_pairs.append(ortho_pair)

    fh.close()

    return ortho_pairs

def get_TE_familyID(fields: str):
    """take the last field in the gff & get the family ID"""

    subfields = fields.split(';')
    family = '' # default to empty string

    for subfield in subfields:
        subfield_split = subfield.split('=')
        if subfield_split[0].upper() == "ID":
            family = subfield_split[1]
            break
    
    return family
    

def parse_TEs_gff(gff_file: str, sp_lengths: dict):
    """get the TEs for the chromosomes that were retained after filtering for length"""

    sp_TEs = defaultdict(list)

    fh = gzip.open(gff_file, "rt") if gff_file.endswith(".gz") else open(gff_file, 'r')

    line_cnt = 0

    for line in fh:
        line_cnt += 1
        if line[0] == '#':
            continue
        fields = line.strip().split('\t')
        if fields[0] in sp_lengths:
            chrom = fields[0]
            start_pos = int(fields[4])
            end_pos = int(fields[5])
            subfamily = fields[2]
            family = get_TE_familyID(fields[8])
            if family == '': # did not find a family
                msg = f"Error: Could not identify family ID at line {line_cnt}\n{line}"
                sys.exit(msg)
            orientation = fields[6]
            te = TE(family, subfamily, start_pos, end_pos, orientation)
            sp_TEs[chrom].append(te)

    fh.close()

    return sp_TEs 

def same_family(familyA: str, familyB: str):
    """return whether the two families are the same between the TEs"""

    # remove the trailing interger
    te_famA = '-'.join(familyA.split('-')[:-1])
    te_famB = '-'.join(familyB.split('-')[:-1])

    return te_famA == te_famB

def get_window_bounds(te_list: list, gene_start: int, gene_end: int, max_dist: int):
    """return the indexes of TEs that fall within the window along with the index closest to the end of the gene"""

    # using 0 is not a safe default
    left, right = -1, -1

    leftbound = gene_start - max_dist
    rightbound = gene_end + max_dist

    # to keep track of the closest TE from the end of the gene [current anchor]
    closest_index = None
    cur_min_dist = float("inf")

    for i, te in enumerate(te_list):
        if (leftbound <= te.end <= rightbound) or (leftbound <= te.start <= rightbound):
            if left == -1: # assign first one
                left = i
                right = i # default just in case there's only 1
            # otherwise, we keep updating our right pointer
            else:
                right = i
            start_dist = abs(gene_end - te.start)
            end_dist =  abs(gene_end - te.end)
            if min(start_dist, end_dist) <= cur_min_dist: # equal to shift closer to the end
                closest_index = i
                cur_min_dist = min(start_dist, end_dist)
            
    return left, right, closest_index


def scan_window(start_a: int, end_a, start_b : int, end_b, chromA: list, chromB: list, closestA: int, closestB: int):
    """scan through the sections of the window to assign orthologous TEs"""

    # first, let's assume that the left most family member is the first
    # i.e., tandem duplicates are right flanking

    cur_familyA = chromA[closestA].family

    while closestA >= 1:
        if chromA[closestA - 1].family == cur_familyA:
            closestA -= 1
        else:
            break
    
    cur_familyB = chromB[closestB].family

    while closestB >= 1:
        if chromB[closestB - 1].family == cur_familyB:
            closestB -= 1
        else:
            break

    # now that the left most family member is located
    # search from A to B if A has more TEs, else
    # search from B to A

    if (end_a - closestA) >= (end_b - closestB):
        first_start = closestA
        first_end = end_a
        subject1 = chromA
        second_start = closestB
        second_end = end_b
        subject2 = chromB
        # for out of bound checks
        numTEs_sub1 = len(chromA)
        numTEs_sub2 = len(chromB)

    else:
        first_start = closestB
        first_end = end_b
        subject1 = chromB
        second_start = closestA
        second_end = end_a
        subject2 = chromA
        numTEs_sub1 = len(chromB)
        numTEs_sub2 = len(chromA)

    for i in range(first_start, first_end + 1):
        te_sub1 = subject1[i]
        if te_sub1.assigned:
            continue
        for j in range(second_start, second_end + 1):
            te_sub2 = subject2[j]
            if te_sub2.assigned:
                continue
            # assign left most
            if same_family(te_sub1.family, te_sub2.family):
                # check for co-orthlogous
                te_sub1.assign(j)
                te_sub2.assign(i)
                if (i < numTEs_sub1 - 1) and (j < numTEs_sub2 - 1): # out of window legal
                    if te_sub1.family != subject1[i + 1].family:
                        k = j + 1
                        while (k < numTEs_sub2 - 1) and (te_sub2.family == subject2[k].family):
                            te_sub1.assign(k)
                            subject2[k].assign(i)
                            k += 1
                    # check the other way around
                    elif te_sub2.family != subject2[j + 1].family:
                        k = i + 1
                        while (k < numTEs_sub1 - 1) and (te_sub1.family == subject1[k].family):
                            te_sub2.assign(k)
                            subject1[k].assign(j)
                            k += 1
            # else:
            #     print(te_sub1.family, te_sub2.family)
    

def identify_ortho_TEs(ortho_pairs: list, spA_Tes: dict, spB_Tes: dict, max_dist: int):
    """assume orthology based on neighboring genes"""

    TE_ortho_pairs = []

    # go through each gene to get the chromosome IDs
    for gene_ortho_pair in ortho_pairs:
        chromA = gene_ortho_pair.spA_chrom
        chromB = gene_ortho_pair.spB_chrom
        chromA_TEs = spA_Tes[chromA]
        chromB_TEs = spB_Tes[chromB]
        if (len(chromA_TEs) == 0) or (len(chromB_TEs) == 0):
            continue
        # okay, there are TEs here
        geneA_start = gene_ortho_pair.start_posA
        geneA_end = gene_ortho_pair.end_posA
        geneB_start = gene_ortho_pair.start_posB
        geneB_end = gene_ortho_pair.end_posB
        # identify window within the range for both species
        leftA, rightA, closetA = get_window_bounds(chromA_TEs, geneA_start, geneA_end, max_dist)
        leftB, rightB, closetB = get_window_bounds(chromB_TEs, geneB_start, geneB_end, max_dist)
        if (closetA == None) or (closetB == None):
            continue # no TEs within the window
        # scan and do the assignments
        scan_window(leftA, rightA, leftB, rightB, chromA_TEs, chromB_TEs, closetA, closetB)

    for a in chromA_TEs:
        if a.assigned:
            print(a.ortho_indices)


    return TE_ortho_pairs

def main():
    """Entry point to the program to detect putative orthologous TEs"""

    # get args
    mem_tsv, spA_gff, spB_gff, chroms_tsv, min_len, max_dist \
          = get_arguments()

    # start with the membership file to get the species IDs
    # and by getting the desired scaffolds/chromosomes
    spA_lengths, spB_lengths, spA, spB = \
        parse_chromosomes_tsv(chroms_tsv, mem_tsv, min_len)
    
    # get the orthologous genes
    ortho_pairs = parse_memberships_tsv(mem_tsv, spA, spB)

    # now get the TEs
    spA_TEs = parse_TEs_gff(spA_gff, spA_lengths)
    spB_TEs = parse_TEs_gff(spB_gff, spB_lengths)

    # now to do the ortholog detection
    TE_ortho_pairs = identify_ortho_TEs(ortho_pairs, spA_TEs, spB_TEs, max_dist)

if __name__ == "__main__":
    main()
