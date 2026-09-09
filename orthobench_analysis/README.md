# This folder contains the scripts used in benchmarking Synolog's orthogroup using the Orthobench dataset

For this part of the analysis, I used the data available [here](https://github.com/davidemms/Open_Orthobench). Since `synolog` reports orthologs using gene IDs instead of isoform/transcript IDs (given that `synolog` uses all available isoforms/transcripts), I wrote a python script to change/recode the reference orthogroups used for benchmarking (these, by default use protein IDs). Also, not all protein IDs in these reference orthogroups were found in the gene annotations I used (Ensembl Release 99), so I had to search [ENSEMBL](https://www.ensembl.org/) to find gene IDs when possible. I used these as an extra mapping file for recoding.

```
./recodeRefOGs.py -g ./path/to/gpMaps/ -r ./path/to/RefOGs/ -e extra_gpMap.tsv -o ./outpath/RefOGs.Recoded/
```

I also modified the [benchmarking script](https://github.com/davidemms/Open_Orthobench/blob/master/BENCHMARKS/benchmark.py) so that it can work at the gene-level instead of protein-level (i,e., adjusted it to the recorded RefOGs) and can leverage the gpMap files instead of the *GTF* files (speeds up parsing). Since I can't simply run the *orthologs.tsv* file from `synolog` through this script, I have another helper script to reformat it

```
# reformat
./synolog_to_benchmark.py -o orthologs.tsv # generates benchmarking.orthologs.txt

# benchmark
./benchmarking_gene_level.py ../synolog/bit.std.2.0/benchmarking.orthologs.txt
```

**NOTE**: The benchmarking script assumes the directory it is executed in contains the path to the recoded RefOGs and the path to the gpMaps.
