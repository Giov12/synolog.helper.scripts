# Comparison to OrthoFinder3

For the paper, but I wrote several scripts to compare the results of `synolog` to that from [orthofinder3](https://github.com/OrthoFinder/OrthoFinder).
The first is to convert (i.e., reformat) the isoform IDs to gene IDs, since `synolog` reports gene IDs in its outputs. For this step, I used *gpMAP* files, which can be created as shown below
```
./make_gtf_map.py -g file.gtf.gz
```

To reformat, I fed the `reformatOrthoFinder3.py` script the output map, the results from `orthofinder3`, a single column list with the species ids to use instead of the file ids `orthofinder3` uses to label its columns, and the path to all the gene/protein maps created with `make_gtf_map.py`. This will create a new reformated orthogroups file that uses gene IDs instead of protein sequence IDs.
```
./reformatOrthoFinder3.py -o Orthogroups.tsv -s spp_list.txt -m ./path/to/gene_maps/
```

I then took the resulting file, along with `synolog`'s orthologs.tsv output file, and process them through `compareToReformated.py` to compare their results (e.g., how many orthogroups are the same). I also gave this script the path to the annotations since `synolog`'s ortholog inference method is synteny-based. This helped me tease apart why the two tools would produce different results (e.g., were non-syntenic paralogs causing discrepancies).
```
./compareToReformated.py -m reformated_orthofinder3.tsv -g ./path/to/annotations/ -s orthologs.tsv
```
