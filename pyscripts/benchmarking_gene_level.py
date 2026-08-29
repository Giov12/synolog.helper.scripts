#!/usr/bin/env python3

############################################################################################
#    NOTE: This version compares orthogroups at the GENE level rather than the protein     #
#    level used by the original benchmark.py. It was written because Synolog               #
#    considers all annotated isoforms per gene rather than the longest isoform available.  #
#    So the orthogroups of Synolog are reported as gene IDs and not protein IDs. As such   #
#    the RefOGs fed to this script are assumed to contain gene IDs when possible.          #
############################################################################################

"""
This script calculates the benchmarks for an input set of orthogroups

Instructions:
1. Predict the complete set of orthogroups for the genes in the "gpMaps/" directory
2. Write Orthogroups to a file, one orthogroup per line. First line is a header and is
   ignored. Genes can be separated by commas, spaces or tabs. Lines starting with '#' 
   are ignored:
   E.g. "OG0000001: Gene1, Gene2"
   Optionally, each line can start with the name of the orthogroup followed by a colon.
   I.e. any text before the first colon on each line will be ignored.
   This aim of this format is to be flexible, e.g. OrthoFinder .tsv and OrthoMCL output
   files are both valid under this format.
3. Call the script with the the orthogroup filename as the only argument
"""

import os
import re
import sys
import glob
import gzip
import argparse

delims = " |,|\t"

############################################################################################
#    NOTE: these key-value pairs have been changed to now use GENE level IDs               #
#    instead of protein coding IDs                                                         #
#                                                                                          #
#    ALSO: n_genes_total is overwritten in get_expected_genes() so that it                 #
#    can count the number of genes used with the newer annotations used in this            #
#    study                                                                                 #
############################################################################################

expected_gene_names_base = {"ENSG": "Homo sapiens", 
    "ENSRNOG":"Rattus norvegicus", 
    "ENSCAFG":"Canis familiaris", 
    "ENSMUSG":"Mus musculus", 
    "ENSMODG":"Monodelphis domestica", 
    "ENSGALG":"Gallus gallus", 
    "ENSCING":"Ciona intestinalis",
    "ENSTNIG":"Tetraodon nigroviridis",
    "ENSPTRG":"Pan troglodytes",
    "ENSDARG":"Danio rerio",
    "FBgn":"Drosophila melanogaster",
    "WBGene":"C elegans"}

n_genes_total = 251378

def exit():
    sys.exit()

def read_hierarchical_orthogroup(infile):
    ogs = []
    iSkip = 3 
    for line in infile:
        genes = [g for s in line.rstrip().split("\t")[iSkip:] for g in s.split(", ") if g != ""]
        ogs.append(set(genes))
    return ogs


def read_orthogroups(fn, exp_genes, n_col_skip=1):
    """
    Read the orthogroups from a file formatted as specified above
    """
    ogs = []
    q_past_header = False
    with open(fn, 'r') as infile:
        for l in infile:
            if l.startswith("#"):
                continue
            t = l.rstrip().split(None, n_col_skip)[-1]
            genes = re.split(delims, t)
            # remove white spaces
            genes = [g.strip() for g in genes]
            genes = set([g for g in genes if g != ""])
            if q_past_header or (exp_genes & genes):
                q_past_header = True
                ogs.append(genes)
    return ogs


def get_n_col_skip(fn):
    with open(fn, 'r') as infile:
        header = next(infile)
        if header.startswith("HOG\tOG\tGene Tree Parent Clade"):
            return 3
    return 1


def get_expected_genes():

    global n_genes_total

    ############################################################################################
    #    NOTE: This function has been modified to look at all the protein IDs available.       #
    #    So since synolog uses the annotations to associate isoforms with genomic positions,   #
    #    we will use some pre-parsed gene-protein maps                                         #
    ############################################################################################

    d_gp     = os.path.dirname(__file__) + os.sep + "gpMaps" + os.sep
    gp_files = list(glob.glob(d_gp + "*.gpMap.tsv"))
    assert 12 == len(gp_files), f"Number of gene-protein maps not equal to 12 orthobench orgs"

    all_genes = set()
    for fn in gp_files:
        with open(fn, 'r') as infile:
            for l in infile:
                if (len(l) == 0 or l[0] == '#'):
                    continue
                idx = l.find('\t')
                all_genes.add(l[:idx])

    # provide a score for all the genes
    n_genes_total = len(all_genes)

    return all_genes

def check_orthogroups(ogs, exp_genes):
    all_pred_genes = set([g for og in ogs for g in og])
    x = all_pred_genes.difference(exp_genes)
    if len(x) != 0:
        print("ERROR: found extra genes in input file, check its formatting is correct and there are no incorrect genes")
        print("Examples:")
        for g in list(x)[:10]:
            print(g)
    x = exp_genes.difference(all_pred_genes)
    if len(x) != 0:
        print("Examples of genes not in file:")
        for g in list(x)[:3]:
            print(g)

    n_genes = sum([len(og) for og in ogs])
    if n_genes < 0.5 * n_genes_total:
        print("ERROR: Too many missing genes in predicted orthogroups.")
        print("Orthogroups should contain at least 50% of all genes but")
        print("orthogroups file only contained %d genes" % n_genes)
        exit()
    all_genes = [g for og in ogs for g in og]
    n_genes_no_dups = len(set(all_genes))
    if n_genes_no_dups != n_genes:
        print("ERROR: Some genes appear in multiple orthogroups, benchmark are meaningless with such data.")
        print("with such data")
        print((n_genes_no_dups, n_genes))
        from collections import Counter
        c = Counter(all_genes)
        for g, n in c.most_common(10):
            print("%d: %s" % (n,g))
        raise Exception()
    # checked genes from each of the expected species are present
    for g_pat, sp in expected_gene_names_base.items():
        if not any(g.startswith(g_pat) for g in all_genes):
            print("ERROR: No genes found from %s" % sp)
    p = 100.*n_genes/float(n_genes_total)
    print("%d genes found in predicted orthogroups, this is %0.1f%% of all genes" % (n_genes, p))

def read_refogs(d_refogs):

    ############################################################################################
    #    NOTE: This function has been modified to used the Recoded RefOGs files.               #
    #    The recoding is the swapping of protein sequence IDs with gene IDs                    #
    ############################################################################################

    refogs = []
    for i in range(1, 71):
        fn = d_refogs + ("Recoded.RefOG%03d.txt" % i)
        if not os.path.exists(fn):
            print("ERROR: RefOG file not found: %s" % fn)
            exit()
        with open(fn, 'r') as infile:
            refogs.append(set([g.rstrip() for g in infile.readlines()]))
    n = sum([len(r) for r in refogs])
    n_expected = 1945
    if not n_expected == n:
        print("ERROR: There are genes missing from the RefOG files. Found %d, expected %d" % (n, n_expected))
    return refogs


def read_uncertain_refogs(d_refogs):

    ############################################################################################
    #    NOTE: This function has been modified to used the Recoded RefOGs files.               #
    #    The recoding is the swapping of protein sequence IDs with gene IDs                    #
    ############################################################################################

    refogs = []
    for i in range(1, 71):
        fn = d_refogs + ("Recoded.RefOG%03d.txt" % i)
        if not os.path.exists(fn): 
            refogs.append(set())
        else:
            with open(fn, 'r') as infile:
                refogs.append(set([g.rstrip() for g in infile.readlines()]))
    return refogs


def calculate_benchmarks_pairwise(ref_ogs, uncert_genes, pred_ogs, q_even=True, q_remove_uncertain=True):
    referenceOGs = ref_ogs
    predictedOGs = pred_ogs
    totalFP = 0.
    totalFN = 0.
    totalTP = 0.
    totalGroundTruth = 0.
    n_exact = 0
    # so as not to count uncertain genes either way remove them from the 
    # expected and remove them from any predicted OG (as though they never existed!)
    for refOg, uncert in zip(referenceOGs, uncert_genes):
        thisFP = 0.
        thisFN = 0.
        thisTP = 0.
        if q_remove_uncertain:
            refOg = refOg.difference(uncert)
        nRefOG = len(refOg)
        not_present = set(refOg)
        for predOg in predictedOGs:
            overlap = len(refOg.intersection(predOg))
            if overlap > 0:
                if q_remove_uncertain:
                    predOg = predOg.difference(uncert)   # I.e. only discount genes that are uncertain w.r.t. this RefOG
                overlap = len(refOg.intersection(predOg))
            if overlap > 0:
                not_present = not_present.difference(predOg)
                thisTP += overlap * (overlap - 1)/2    # n-Ch-2
                thisFP += overlap * (len(predOg) - overlap)
                thisFN += (nRefOG - overlap) * overlap
        # finally, count all the FN pairs from those not in any predicted OG
        thisFN += len(not_present)*(nRefOG-1)
        # All FN have been counted twice
        assert(thisFN % 2 == 0)
        thisFN /= 2 
        # sanity check
        nPairs1 = thisTP + thisFN
        nPairs2 = nRefOG * (nRefOG - 1) / 2
        if nPairs1 != nPairs2:
            print("ERROR: %d != %d" % (nPairs1,nPairs2))    
            # print(refOg)
            # print(predOg)
            # print((nRefOG, thisTP, thisFP, thisFN))
            # print("ERROR: Sanity check failed")
#            assert(nPairs1 == nPairs2)
#            raise Exception
        totalGroundTruth += nPairs1
        if thisFN == 0 and thisFP == 0:
            n_exact += 1
        if q_even:
            N = float(len(refOg)-1)
            totalFN += thisFN/N
            totalFP += thisFP/N
            totalTP += thisTP/N
        else:
            totalFN += thisFN
            totalFP += thisFP
            totalTP += thisTP
        # print("%d\t%d\t%d" % (thisTP, thisFN, thisFP))  
    TP, FP, FN = (totalTP, totalFP, totalFN)
    # print("%d Correct gene pairs" % TP)
    # print("%d False Positives gene pairs" % FP)
    # print("%d False Negatives gene pairs\n" % FN)
    pres = TP/(TP+FP)
    recall = TP/(TP+FN)
    f = 2*pres*recall/(pres+recall)
    print("%0.1f%% F-score" % (100.*f))
    print("%0.1f%% Precision" % (100.*pres))
    print("%0.1f%% Recall\n" % (100.*recall))
    print("%d Orthogroups exactly correct" % n_exact)
    return [100.*f, 100.*pres, 100.*recall]
   
def benchmark(ogs_filename, d_refogs):
    print("\nReading RefOGs from: %s" % d_refogs)
    ref_ogs = read_refogs(d_refogs)
    ref_ogs_uncertain = read_uncertain_refogs(d_refogs + "low_certainty_assignments/")
    exp_genes = get_expected_genes()
    print("\nReading predicted orthogroups from: %s" % ogs_filename)
    n_col_skip = get_n_col_skip(ogs_filename)
    pred_ogs = read_orthogroups(ogs_filename, exp_genes, n_col_skip)
    check_orthogroups(pred_ogs, exp_genes)
    print("\nCalculating benchmarks:")
    # print(os.path.basename(ogs_filename))
    x = calculate_benchmarks_pairwise(ref_ogs, ref_ogs_uncertain, pred_ogs)
    # print("\t".join(map(str, x)))

if __name__ == "__main__":

    ############################################################################################
    #    NOTE: The main function has been changed to now use the recoded RefOGs paths          #
    ############################################################################################

    parser = argparse.ArgumentParser()
    parser.add_argument("ogs_filename", help="File containing orthogroups.")
    args = parser.parse_args()
    d_refogs = os.path.dirname(__file__) + os.sep + "RefOGs.Recoded" + os.sep
    benchmark(args.ogs_filename, d_refogs)
