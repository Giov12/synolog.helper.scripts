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

More for the paper, but I wrote several scripts to compare the results of `synolog` to that from [orthofinder3](https://github.com/OrthoFinder/OrthoFinder).
The first is to convert (i.e., reformat) the isoform IDs to gene IDs, since `synolog` reports gene IDs in its outputs.

The first was to create a simple two-column tsv file that maps isoform id to gene id.
```
./make_gtf_map.py -g file.gtf.gz
```

In the next step, I fed the `reformatOrthoFinder3.py` script the output map, the results from `orthofinder3`, a single column list with the species ids to use instead of the file ids `orthofinder3` uses to label its columns, and the path to all the gene/protein maps created with `make_gtf_map.py`. This will create a new reformated NO.tsv file that uses the *HOG* ID to identify each group of genes.
```
./reformatOrthoFinder3.py -o Orthogroups.tsv -s spp_list.txt -m ./path/to/gene_maps/
```

I then take the resulting file, along with `synolog`'s orthologs.tsv output file, and process them through `compareToReformated.py` to compare their results (e.g., how many orthogroups are the same). I also give this script the path to the annotations since `synolog`'s ortholog inference method is synteny-based. This helped me tease apart why the two tools would produce different results.
```
./compareToReformated.py -m reformated_orthofinder3.tsv -g ./path/to/annotations/ -s orthologs.tsv
```

Another script that I used for the study is `treeclusters2ALGs.py`, which takes the *tree_clusters_membership.tsv* and counts the number of genes these clusters occupy after merging the clusters by the chromosomes they occupy.
```
./treeclusters2ALGs.py -t tree_clusters_membership.tsv
```

Since Synolog does not know which sequences are actually chromosomes, I manually filtered these results to identify which clusters represented ALGs.

Also, for the orthobench step, since `synolog` reports orthologs using gene IDs instead of isoform/transcript IDs (given that `synolog` uses all available isoforms/transcripts), I wrote a python script to change recode the reference orthogroups used for benchmarking (these, by default use protein IDs). Also, not all protein IDs in these reference orthogroups were found in the gene annotations I used, so I had to search [ENSEMBL](https://www.ensembl.org/) to find gene IDs when possible. I used these as an extra mapping file for recoding.

```
./recodeRefOGs.py -g ./path/to/gpMaps/ -r ./path/to/RefOGs/ -e extra_gpMap.tsv -o ./outpath/RefOGs.Recoded/
```

I also modified the (benchmarking script)[https://github.com/davidemms/Open_Orthobench/blob/master/BENCHMARKS/benchmark.py] so that it can work at the gene-level instead of protein-level (i,e., adjusted it to the recorded RefOGs). Since I can't simply run the *orthologs.tsv* file from `synolog` through this script, I have another helper script to reformat it

```
./synolog_to_benchmark.py -o orthologs.tsv # generates benchmarking.orthologs.txt
```