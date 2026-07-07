#!/bin/env python3

import argparse
import os
import glob
from collections import defaultdict

class UnionFind:
    """Union-Find (Disjoint Set) class to manage gene groupings."""
    def __init__(self):
        self.parent    = {}
        self.rank      = {}
        self.parents   = defaultdict(set)
        self.contra    = defaultdict(set)
        self.groups    = defaultdict(set)
        self.orthogrp  = defaultdict(set)
        self.hasgroup  = set()
        self.newgrps   = defaultdict(set)

    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])  # Path compression
        return self.parent[x]

    def union(self, x, y):
        rootX = self.find(x)
        rootY = self.find(y)

        if rootX != rootY:
            # Union by rank
            if self.rank[rootX] > self.rank[rootY]:
                self.parent[rootY] = rootX
                self.parents[rootX].update(self.parents[rootY])
                # self.rank[rootX] += self.rank[rootY]
            elif self.rank[rootX] < self.rank[rootY]:
                self.parent[rootX] = rootY
                # self.rank[rootY] += self.rank[rootX]
                self.parents[rootY].update(self.parents[rootX])
            else:
                self.parent[rootY] = rootX
                self.parents[rootX].update(self.parents[rootY])
                self.rank[rootX] += 1 # self.rank[rootY]
            
            if (rootX in self.contra):
                del self.contra[rootX]
            if (rootY in self.contra):
                del self.contra[rootY]

            # self.parents[rootY].clear()
    def build_spp_map(self, genes):
        spMap = defaultdict(set)
        for g in genes:
            sp = g.split(':')[0]
            spMap[sp].add(g)
        return spMap

    def read_pair(self, gene_a: str, gene_b: str):
        self.groups[gene_a].add(gene_b)
        self.groups[gene_a].add(gene_a)
        self.groups[gene_b].add(gene_a)
        self.groups[gene_b].add(gene_b)


    def add(self, x: str, spp: str):
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x]   = 0
            self.parents[x].add(spp)
        else:
            self.parents[self.find(x)].add(spp)

    def make_groups(self):
        for gene in self.parent:
            root = self.find(gene)
            if (gene != root):
                self.hasgroup.add(gene)
            self.orthogrp[root].add(gene)
        
        for gene in self.orthogrp:
            if (len(self.orthogrp[gene]) > 1 and gene not in self.hasgroup):
                self.hasgroup.add(gene)

        for gene in self.hasgroup:
            if (gene in self.contra):
                del self.contra[gene]


    def check_contradictions(self, gene_a: str, gene_b: str):
        """Check if merging two genes creates a contradiction."""
        root_a = self.find(gene_a)
        root_b = self.find(gene_b)

        # if (gene_a == "cgun:g13961" or gene_b == "cgun:g13961"):
        #     print(root_a, gene_a, root_b, gene_b)

        if root_a != root_b:
            # Check if merging would cause a species contradiction
            groupa = self.groups[root_a]
            groupb = self.groups[root_b]

            if (len(groupa) > len(groupb)):
                if (groupb.issubset(groupa)):
                    return False
            elif (len(groupa) < len(groupb)):
                if (groupa.issubset(groupb)):
                    return False
            elif (groupa == groupb):
                return False
            else:
                self.contra[gene_a].add(gene_b)
                self.contra[gene_b].add(gene_a)
                return True
            # if not self.parents[root_a].isdisjoint(self.parents[root_b]):
            #     self.contra[gene_a].add(gene_b)
            #     self.contra[gene_b].add(gene_a)
            #     return True
        return False
    
    def merge_contra(self):
        a = 0
        remove = set()
        for gene, cgenes in self.contra.items():
            if (gene in self.hasgroup):
                remove.add(gene)
                continue
            valid = False
            for g in cgenes:
                if (g in self.hasgroup):
                    root = self.find(g)
                    othergenes = self.orthogrp[root]
                    if cgenes.issubset(othergenes):
                        othergenes.add(gene)
                        valid = True
                        break
            if valid:
                a+=1
                remove.add(gene)
        
        for gene in remove:
            del self.contra[gene]

    def write_contradictions(self):
        fh = open("Contradictions.txt", 'w')
        for gene_a, cgenes in self.contra.items():
            cgs = ','.join([g for g in cgenes])
            line = f"Contradiction between {gene_a} and {cgs}\n"
            fh.write(line)
        fh.close()

    def write_orthogroups(self, species_list: list):
        fh = open("Inparanoid.merged.tables.tsv", 'w')

        header = ["#GroupID"] + species_list
        header = '\t'.join(header) + '\n'
        fh.write(header)

        gi = 1 # group id
        for group in self.orthogrp.values(  ):
            row = [f"{gi}"]
            species_dict = {s : [] for s in species_list}
            found = set()
            for gene in group:
                species, gene_id = gene.split(':', 1)
                species_dict[species].append(gene_id)
                found.add(species)
            if (len(found) == 1):
                continue
            for spp in species_list:
                row.append(', '.join(species_dict[spp]) if species_dict[spp] else '')
            line = '\t'.join(row) + '\n'
            fh.write(line)
            gi += 1
        fh.close()



def get_arguments() -> tuple:
    """get the arguments"""

    d = "merge all the table*.faa files from inparanoid into a tsv where all the genes agree with one another"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-d", "--dir", help="directory containing all the table*.faa files", required=True)
    parser.add_argument("-m", "--maps", help="directory containing any map gene->protein id mappping", default='')
    args = parser.parse_args()
    idir = args.dir
    mdir = args.maps

    assert os.path.isdir(idir), f"Could not locate directory {idir}"

    if (mdir != ''):
        assert os.path.isdir(mdir), f"Could not locate directory {mdir}"

    return (idir, mdir)

def get_tables(idir: str) -> list:
    """return the list of input files"""

    ifiles = list()

    if (idir.endswith('/')):
        ifiles = glob.glob(f"{idir}table.*")
    else:
        ifiles = glob.glob(f"{idir}/table.*")

    assert len(ifiles) > 0, f"Could not locate any table.* files at {idir}"

    return ifiles

def get_maps(mdir: str) -> dict:
    """get any gene->protein mapping for easier comparison"""

    if (mdir == ''): return {}

    mfiles = []

    if (mdir.endswith('/')):
        mfiles = glob.glob(f"{mdir}*gpMap.tsv")
    else:
        mfiles = glob.glob(f"{mdir}/*gpMap.tsv")

    assert len(mfiles) > 0, f"Could not locate any *gpMap.tsv files at {mdir}"

    gpMap = {}

    for f in mfiles:
        fh         = open(f, 'r')
        spp        = os.path.basename(f).split('.')[0]
        gpMap[spp] = {}

        for line in fh:
            fields = line.strip().split('\t')
            gene   = fields[0]
            prot   = fields[1]
            gpMap[spp][prot] = gene

        fh.close()

    return gpMap

def get_spp_names(fname: str) -> tuple:
    """get the spp names from the file name"""

    bname = os.path.basename(fname)
    bname = bname.replace("table.", '')
    bname = bname.replace(".faa", '')
    spp   = bname.split('-')

    return (spp[0], spp[1])


def create_imap(ifiles: list, gpMap: dict) -> tuple:
    """go file by file, building upon the ortholog group"""

    # sort so we encounter the first spp >1 time
    ifiles.sort()

    imap = {} # iparanoid mapping
    taxa = []

    # next, read in all the data before we go in & merge
    for f in ifiles:
        spp  = get_spp_names(f)
        
        if (spp[0] not in imap):
            imap[spp[0]] = {}
            taxa.append(spp[0])
        if (spp[1] not in imap):
            imap[spp[1]] = {}
            taxa.append(spp[1])

        fh = open(f, 'r')
        for line in fh:
            if (line[0] == 'O'):
                continue
            fields = line.split('\t')
            genes1 = fields[2].strip()
            genes2 = fields[3].strip()
            fields = [genes1, genes2]
            add    = [[], []]

            for i, genes in enumerate(fields):
                genes = genes.split(' ')
                org  = spp[i]
                for g in genes:
                    if (g[0].isnumeric()): continue
                    g = g.split('.')[0]
                    if (org in gpMap):
                        g = gpMap[org][g]
                    if (g not in imap[org]):
                        imap[org][g] = defaultdict(set)
                    add[i].append(g)

            for i, genes in enumerate(add):
                j = 1 if (i == 0) else 0
                a = spp[i]
                b = spp[j]
                for g1 in genes:
                    for g2 in genes:
                        if (g1 != g2):
                            imap[a][g1][a].add(g2)
                    for g2 in add[j]:
                        imap[a][g1][b].add(g2)

        fh.close()

    return (imap, taxa)

# Step 1: Parsing Functions
def extract_species_from_filenames(pair_files: list) -> list:
    """Extract unique species names from filenames."""
    species = set()
    for file in pair_files:
        filename = os.path.basename(file)
        species_a, species_b = filename.split('.')[1:3]
        species_b = species_b.split('-')[1]
        species.add(species_a)
        species.add(species_b)
    return sorted(species)  # Sort for consistent column ordering

def parse_genes_with_species(column: str, species: str, gpMap: dict) -> list:
    """Parse a column with multiple genes and prepend species names."""
    column = column.split()
    genes  = []
    for i in range(0, len(column), 2):
        gene = column[i].split('.')[0]
        if (species in gpMap):
            gene = gpMap[species][gene]
        genes.append(f"{species}:{gene}")

    return genes

def parse_input(file_path: str, gpMap: dict) -> list:
    """Parse a file and include species names as prefixes to gene IDs."""
    filename = os.path.basename(file_path)
    species_a, species_b = filename.split('.')[1:3]
    species_b = species_b.split('-')[1]

    pairs = []
    with open(file_path, 'r') as f:
        for line in f:
            if (line[0] == 'O'):
                continue
            fields = line.split('\t')
            genes1 = fields[2].strip()
            genes1 = parse_genes_with_species(genes1, species_a, gpMap)
            genes2 = fields[3].strip()
            genes2 = parse_genes_with_species(genes2, species_b, gpMap)
            for gene_a in genes1:
                for gene_b in genes2:
                    pairs.append((gene_a, gene_b))
            if (len(genes1) > 1):
                for g in genes1:
                    for g2 in genes1:
                        if (g != g2):
                            pairs.append((g, g2))
            if (len(genes2) > 1):
                for g in genes2:
                    for g2 in genes2:
                        if (g != g2):
                            pairs.append((g, g2))
                

    return pairs

def is_same_spp(gene_a: str, gene_b: str) -> bool:
    """if the two genes are from the same species"""
        
    species_a, gene_a_id = gene_a.split(':', 1)
    species_b, gene_b_id = gene_b.split(':', 1)

    return species_a == species_b
    

# Step 2: Graph Construction
def build_UF(pair_files: list, gpMap: dict, species_list: list) -> dict:
    """Build an unweighted graph from pairwise ortholog files, with species-specific gene IDs."""
    
    uf = UnionFind()

    # fh = open("Contradictions.txt", 'w')

    all_pairs = []
    try_again = []
    pair_files.sort()
    for file_path in pair_files:
        pairs = parse_input(file_path, gpMap)
        new_pairs = []
        for a, b in pairs:
            uf.read_pair(a, b)
            if not is_same_spp(a, b):
                new_pairs.append((a, b))
        all_pairs.append(new_pairs)
        
    for pairs in all_pairs:
        for gene_a, gene_b in pairs:
            uf.read_pair(gene_a, gene_b)
            
            species_a, gene_a_id = gene_a.split(':', 1)
            species_b, gene_b_id = gene_b.split(':', 1)
            
            uf.add(gene_a, species_a)
            uf.add(gene_b, species_b)

            if uf.check_contradictions(gene_a, gene_b):
                try_again.append((gene_a, gene_b))
                # line = f"Contradiction between {gene_a} and {gene_b}\n"
                # fh.write(line)
            else:
                uf.union(gene_a, gene_b)

    # fh.close()

    for gene_a, gene_b in try_again:
        species_a, gene_a_id = gene_a.split(':', 1)
        species_b, gene_b_id = gene_b.split(':', 1)
            
        uf.add(gene_a, species_a)
        uf.add(gene_b, species_b)

        if not uf.check_contradictions(gene_a, gene_b):
            uf.union(gene_a, gene_b)
        # else:
        #     line = f"Contradiction between {gene_a} and {gene_b}\n"
        #     fh.write(line)

    uf.make_groups()
    uf.merge_contra()
    uf.write_contradictions()
    uf.make_groups()
    uf.write_orthogroups(species_list)

    # fh.close()

    return uf


# Step 4: Connected Components
def find_connected_components(uf: UnionFind) -> dict:
    """Find connected components in UF structure."""
    components = defaultdict(list)

    for gene in uf.parent:
        root = uf.find(gene)
        components[root].append(gene)

    return list(components.values())

# Step 5: Format Output
def write_orthogroups(groups: list, species_list: list) -> None:
    """Format ortholog groups into a table with species columns."""
    
    fh     = open("Merged.tables.tsv", 'w')
    header = ["#GroupID"] + species_list
    header = '\t'.join(header) + '\n'
    fh.write(header)

    for group_id, group in enumerate(groups, start=1):
        row = [f"{group_id}"]
        species_dict = {species: [] for species in species_list}
        for gene in group:
            species, gene_id = gene.split(':', 1)
            species_dict[species].append(gene_id)
        for species in species_list:
            row.append(', '.join(species_dict[species]) if species_dict[species] else '')
        line = '\t'.join(row) + '\n'
        fh.write(line)
        
    fh.close()

def main() -> int:
    """get the gtf file and create a tsv file to serve as a map"""

    # get arguments & do everything
    idir, mdir = get_arguments()
    gpMap      = get_maps(mdir)
    ifiles     = get_tables(idir)
    spList     = extract_species_from_filenames(ifiles)
    build_UF(ifiles, gpMap, spList)


    return 0

if __name__ == "__main__":
    main()
