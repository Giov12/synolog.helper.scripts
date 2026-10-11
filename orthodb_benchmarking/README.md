# OrthoDB v11 Benchmarking
This repository contains the workflow used when benchmarking the Synolog's orthogroups to those downloaded in our benchmarking analysis using OrthoDB v11 orthogroups as a reference.

The orthoDB data was obtained from [here](https://data.orthodb.org/v12/download/odb_data_dump/) with much of the code used for filtering the data found in `../cpp.files/`. For our first step, we filtered the  `odb12v2_OGs.tab.gz` file to retain orthogroups at the vertebrata taxonomical level
```
# find desired level ID
 $ zcat odb12v2_levels.tab.gz | awk -F'\t' 'tolower($2)=="vertebrata"'
7742    Vertebrata      21291707        56137   1155

# filter f
zcat odb12v2_OGs.tab.gz | awk -F'\t' '{if ($2==7742) print}' | gzip > vertebrata_OGs.tab.gz
```

For this case study, the following species were selected

| Species | Taxon ID |
|---|---|
| *Lepisosteus oculatus* | 7918_0 |
| *Danio rerio* | 7955_0 |
| *Oryzias latipes* | 8090_0 |
| *Xenopus tropicalis* | 8364_0 |
| *Anolis carolinensis* | 28377_0 |
| *Gallus gallus* | 9031_0 |
| *Monodelphis domestica* | 13616_0 |
| *Bos taurus* | 9913_0 |
| *Canis lupus familiaris* | 9615_1 |
| *Mus musculus* | 10090_1 |
| *Homo sapiens* | 9606_1 |

Orthogroups were further filtered to only retain orthologs from these species. For this procedure, we filtered the `odb12v2_OG2genes.tab.gz` file.
```
$ head -n 3 species_ids.txt
9606_1
7955_0
7918_0

# grab orthogroup IDs
$ zcat vertebrata_OGs.tab.gz | cut -f1 | gzip > vertebrata_OGs_ids.txt.gz

# filter to the specific set of orthogroups
$ ./filter_og2_genes -s species_ids.txt -g vertebrata_OGs_ids.txt -o odb12v2_OG2genes.tab.gz

# repeat with the orthogroup genes file
$ ./filter_odb_genes_by_spp -s species_ids.txt -g odb12v2_genes.tab.gz
```

The two filtered files were merged into a single `.tsv` file, which the taxon IDs replaced with the same species IDs used with synolog.
```
$ head -n 3 spp_map.tsv
9606_1  hsap
7955_0  drer
7918_0  locu

$ ./make_orthodb_refOGs -G filtered.odb12v2_OG2genes.tab.gz -g filtered.odb12v2_genes.tab.gz -o orthodb_orthogroups.tsv.gz -s spp_map.tsv
```
Lastly, we renamed the gene IDs by swapping out the OrthoDB gene IDs with those found in each species respective NCBI genome annotations

```
# directory containing all genome annotations with files renamed following Synolog's species cache's naming convention
tree gtfs/
gtfs/
├── acar.def.gtf.gz
├── btau.def.gtf.gz
├── cfam.def.gtf.gz
├── drer.def.gtf.gz
├── ggal.def.gtf.gz
├── hsap.def.gtf.gz
├── locu.def.gtf.gz
├── mdom.def.gtf.gz
├── mmus.def.gtf.gz
├── olat.def.gtf.gz
└── xtro.def.gtf.gz

# rename orthologs
$ ./find_ensembl_ids -o orthodb_orthogroups.tsv.gz -a gtfs/
```

The renamed orthogroups were then compared to Synolog's `orthologs.tsv` file using `../pyscripts/benchmark_synolog.py`
```
$ ./benchmark_synolog.py -o renamed.orthodb_orthogroups.tsv.gz -s orthologs.tsv
```
