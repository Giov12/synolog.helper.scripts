#!/bin/env python3

import argparse
import os
import gzip


### NOTE:
###      this is going to be a bit ugly
### TODO: 
#        come back in the future and have a class do the computation

def get_arguments():
    """get the gff file"""

    parser = argparse.ArgumentParser(
        description="Convert a gff3 file to a gtf file")
    parser.add_argument(
        "-g", "--gff", help="gff file [required]", required=True)

    args = parser.parse_args()
    gff  = args.gff 

    assert os.path.isfile(gff), f"Could not find {gff}"

    return gff

def make_gtf(gff: str):
    """create a gtf file for the gff file"""

    bname = os.path.basename(gff)
    ext   = ".gff3" if (".gff3" in bname) else ".gff"
    ext   = ext + ".gz" if bname.endswith(".gz") else ext
    ofh   = open(bname.replace(ext, ".gtf"), 'w')
    fh    = gzip.open(gff, "rt") if bname.endswith(".gz") else open(gff, 'r')
    gi    = '' # current gene id
    giNum = -1
    prev  = ''
    tid   = '' # transcript id

    # write some additional info
    ofh.write(f"##gff3gtf.py -g {gff}\n")

    # now parse the gff
    for lineNum, line in enumerate(fh):
        if (len(line) == 0) or (line[0] == '\t'):
            continue
        if line[0] == '#':
            ofh.write(line)
            continue
        fields = line.strip('\n').split('\t')
        ### gene
        if (fields[2] == "gene"):
            subfields = fields[-1].split(';')
            newfields = []
            found = ("gene_id" in fields[-1])
            for subfield in subfields:
                info = subfield.split('=')
                if (len(info) == 1):
                    info = subfield.split(' ')
                if (info[0] == "ID"):
                    if (found == False):      
                        info[0] = "gene_id"
                        gi      = info[1].strip(' "')
                        giNum   = lineNum
                    else:
                        continue
                elif (info[0] == "gene_id"):
                    if (info[1].strip('"').startswith("ABC")):
                        continue
                    gi    = info[1]
                    giNum = lineNum
                elif (info[0] == "Name"):
                    info[0] = "gene_name"
                elif (info[0] == "transcript_id" and info[1] != tid):
                    tid = info[1]
                info[1] = f"\"{info[1]}\""
                info = ' '.join(info)
                prev = info
                newfields.append(info)
            fields[-1] = "; ".join(newfields) + ';'
        elif (fields[2] == "transcript"):
            subfields = fields[-1].split(';')
            newfields = []
            if (lineNum == giNum + 1):
                newfields.append(prev)
            for subfield in subfields:
                info = subfield.split('=')
                if (info[0] == "ID"):
                    info[0] = "transcript_id"
                elif (info[0] == "Name"):
                    info[0] = "gene_name"
                elif (info[0] == "Parent"):
                    info[1] = gi
                info[1] = f"\"{info[1]}\""
                info = ' '.join(info)
                newfields.append(info)
            fields[-1] = "; ".join(newfields) + ';'
        ### mRNA -> transcriptq
        elif (fields[2] == "mRNA"):
            fields[2] = "transcript"
            subfields = fields[-1].split(';')
            newfields = []
            for subfield in subfields:
                info = subfield.split('=')
                if (info[0] == "ID"):
                    info[0] = "transcript_id"
                    if (tid != ''):
                        info[1] = tid
                elif (info[0] == "Name"):
                    info[0] = "gene_name"
                if (info[1].startswith("mRNA.")):
                    info[1] = info[1].replace("mRNA.", '')
                info[1] = f"\"{info[1]}\""
                info = ' '.join(info)
                newfields.append(info)
            newfields.append(f"gene_id \"{gi}\";")
            fields[-1] = "; ".join(newfields)
        ### CDS -> exon
        elif (fields[2].lower() in ["cds", "exon"]):
            fields[2] = "exon" if fields[2].lower()[0] == 'e' else "CDS"
            subfields = fields[-1].split(';')
            newfields = []
            gene_id   = False
            for subfield in subfields:
                info = subfield.split('=')
                if (info[0] == "ID"):
                    info[0] = "exon_id" if (fields[2].lower()[0] == 'e') else "cds_id"
                    if (info[1].lower() in ["mrna", "cds", "exon"]):
                        info[1] = gi
                elif (info[0] == "Name"):
                    info[0] = "gene_name"
                if (info[1].startswith("cds.")):
                    info[1] = info[1].replace("cds.", '')
                if (info[0] == "gene_id"):
                    gene_id = True
                if (info[0] == "Parent" and tid != ''):
                    info[1] = tid
                for m in [".t1.exon", ".t1.cds", ".exon", ".cds"]:
                    if (info[1].endswith(m)):
                        info[1] = info[1].replace(m, '')
                    elif (m in info[1]):
                        idx = info[1].find(m)
                        info[1] = info[1][:idx]
                if (info[1].startswith("mRNA.")):
                    info[1] = info[1].replace("mRNA.", '')
                info[1] = f"\"{info[1]}\""
                info = ' '.join(info)
                newfields.append(info)
            if (gene_id == False):
                newfields.append(f"gene_id \"{gi}\";")
            fields[-1] = "; ".join(newfields)
        elif (fields[2].endswith("prime_UTR")):
            subfields = fields[-1].split(';')
            newfields = []
            for subfield in subfields:
                info = subfield.split('=')
                if (info[0] == "ID"):
                    info[0] = "transcript_id"
                if (info[1].startswith("mRNA.")):
                    info[1] = info[1].replace("mRNA.", '')
                if (".cds" in info[1]):
                    info[1] = info[1].replace(".cds", '')
                if (".exon" in info[1]):
                    info[1] = info[1].replace(".exon", '')
                info[1] = f"\"{info[1]}\""
                info = ' '.join(info)
                newfields.append(info)
            newfields.append(f"gene_id \"{gi}\";")
            fields[-1] = "; ".join(newfields)
        newline = '\t'.join(fields)
        newline = newline + '\n'
        ofh.write(newline)

    fh.close()
    ofh.close()

def main():
    """start the program"""

    # get input
    gff = get_arguments()

    # write output
    make_gtf(gff)
    
if __name__ == "__main__":
    main()
