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

More for the paper, but I wrote several scripts to compare the results of `synolog` to that from [orthofinder2](https://github.com/davidemms/OrthoFinder).
The first is to convert (i.e., reformat) the isoform IDs to gene IDs, since `synolog` reports gene IDs in its outputs.

The first was to create a simple two-column tsv file that maps isoform id to gene id.
```
./make_gtf_map.py -g file.gtf.gz
```

In the next step, I fed the `reformatOrthoFinder.py` script the output map, the results from `orthofinder2`, a single column list with the species ids to use instead of the file ids `orthofinder2` uses to label its columns, and the path to all the gene/protein maps created with `make_gtf_map.py`. This will create a new reformated NO.tsv file that uses the *HOG* ID to identify each group of genes.
```
./reformatOrthoFinder.py -o N0.tsv -s spp_list.txt -m ./path/to/gene_maps/
```

I then take the resulting file, along with `synolog`'s orthologs.tsv output file, and process them through `compareToReformated.py` to compare their results (e.g., how many orthogroups are the same). I also give this script the path to the annotations since `synolog`'s ortholog inference method is synteny-based. This helped me tease apart why the two tools would produce different results.
```
./compareToReformated.py -m reformated_orthofinder2.tsv -g ./path/to/annotations/ -s orthologs.tsv
```