# synolog.helper.scripts
This repository contains a bunch of helper scripts used during my journey on developing and working with [synolog](https://catchenlab.life.illinois.edu/synolog/)

`synolog` is not hosted here, but instead on the Catchen's lab website. So I am not sure what scripts can be used by other users for their own analyses, but I can list a few that have the potential to.

For example, maybe one would like to merge a bunch of annotations in *GTF* format into a single file? I wrote `mergeGtfs.py` and used it when initially having non-coding and coding gene annotations in seperate files. That code is now embedded in `synolog`'s `synolog_ncgenes.py` utility.
```
./mergeGTFs.py --gtfs file1.gtf.gz file2.gtf.gz ...
```

Another potential useful script is `ExtractLongestTranscript.py`. Several tools assume a gene will be represented by a single isoform/transcript. A common approach is to use the longest (based on sequence length) and so this script does that by parsing the corresponding annotations to identify the longest isoform/transcript for each gene and grabs them out of the protein fasta file.
```
./ExtractLongestTranscript.py --ann ann.file.gtf.gz --fasta prot.file.fa.gz
```

Several of the scripts/analyses I did relied on a *gpMAP* file, which are just two column *.tsv* files with one column being the gene ID, and the second column the protein ID. These were constructed by parsing *GTF* files, with their purpose simply for speeding up parsing out this information in downstream analyses. These can be constructed using `pyscripts/make_gtf_map.py`
```
./make_gtf_map.py -g file.gtf.gz
```

Another script that I used for the study is `treeclusters2ALGs.py`, which takes the *tree_clusters_membership.tsv* and counts the number of genes these clusters occupy after merging the clusters by the chromosomes they occupy.
```
./treeclusters2ALGs.py -t tree_clusters_membership.tsv
```

Since Synolog does not know which sequences are actually chromosomes, I manually filtered these results to identify which clusters represented ALGs.