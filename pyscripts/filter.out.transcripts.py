#!/bin/env python3

import gzip
import argparse
import os
from collections import defaultdict

def get_arguments():
    """function to return the path to the orthoMCL out file file"""

    parser = argparse.ArgumentParser(description="Filter genes and groups for alternative transcripts")
    parser.add_argument("-O", "--orthologs", required=True, type=str)
    parser.add_argument("-o", "--output", required=False, type=str, default="Filtered_ortho_groups.out")

    args       = parser.parse_args()
    ortho_file = args.orthologs
    output     = args.output

    assert (os.path.isfile(ortho_file)), f"Could not locate {ortho_file}"

    return ortho_file, output

def parse_gene(gene: str):
    """remove transcript tag"""

    if (gene.count('.') == 0):
        return gene

    sections = gene.split('.')
    last_sec = sections[-1]

    if (last_sec[0] == 't'):
        last_sec = last_sec[1:]

    if (last_sec.isdigit()):
        sections.pop()

    return '.'.join(sections)



def filter_orthoFile(ortho_file: str, output: str):
    """one function to go line by line & filter genes/groups"""

    fh  = gzip.open(ortho_file, "rt") if (ortho_file.endswith(".gz")) else open(ortho_file, 'r')
    ofh = open(output, 'w')

    for line in fh:
        sections = line.strip('\n').split(' ')
        newline  = [sections[0]] # get group id
        gene_map = defaultdict(set)
        for entry in sections[1:]:
            spp, gene = entry.split('|')
            gene      = parse_gene(gene)
            gene_map[spp].add(gene)
        if (len(gene_map) == 1):
            genes = list(gene_map.values())[0]
            if (len(genes) == 1):
                continue
        for spp, genes in gene_map.items():
            out   = spp + '|'
            genes = ','.join(list(genes))
            newline.append(out + genes)
        newline = ' '.join(newline) + '\n'
        ofh.write(newline)

    fh.close()
    ofh.close()

def main():
    """entry point to the pipeline"""

    # get the inputs
    ortho_file, output = get_arguments()

    # parse the file
    filter_orthoFile(ortho_file, output)


if __name__ == "__main__":
    main()