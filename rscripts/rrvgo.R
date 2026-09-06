# fix the working directory
wd = "/path/to/synolog/output"
setwd(wd)

library(biomaRt)
library(tidyverse)
library(rrvgo)

go_analysis = function(input_file, out_file){
  
  orthos = read.csv(input_file, sep = '\t', header = FALSE, col.names = c("id", "genes"))
  # map_chr pulls the first available ID
  orthos = orthos %>% filter(!is.na(genes) & genes != '') %>% mutate(gene = str_split(genes, ',') %>% map_chr(1)) %>% pull(gene)
  
  #
  # create our connection to the db
  # use listDatasets(useMart("ensembl")) to see what's available
  # 43           cpbellii_gene_ensembl          Painted turtle genes (Chrysemys_picta_bellii-3.0.3)      Chrysemys_picta_bellii-3.0.3
  # 73            ggallus_gene_ensembl                  Chicken genes (bGalGal1.mat.broiler.GRCg7b)       bGalGal1.mat.broiler.GRCg7b
  # 71          gevgoodei_gene_ensembl             Goodes thornscrub tortoise genes (rGopEvg1_v1.p)                     rGopEvg1_v1.p
  # 157         psinensis_gene_ensembl                  Chinese softshell turtle genes (PelSin_1.0)                        PelSin_1.0
  # 
  mart      = useMart("ensembl", dataset = "gevgoodei_gene_ensembl")
  attrb     = c("external_gene_name", "go_id", "name_1006", "namespace_1003")
  db        = "org.Gg.eg.db"
  goBM      = getBM(attributes = attrb, filters = "external_gene_name", values = orthos, mart = mart) %>% as_tibble()
  terms     = goBM %>% filter(namespace_1003 == "biological_process" & !is.na(go_id) & go_id != '')
  counts    = terms %>% count(go_id) %>% deframe()
  simMatrix = calculateSimMatrix(terms$go_id, orgdb = db, method="Rel")
  reduced   = reduceSimMatrix(simMatrix, counts, threshold = 0.7, orgdb = db)
  
  #
  # set up the layout for plotting
  #
  scatterPlot(simMatrix, reduced, onlyParents = FALSE)
  heatmapPlot(simMatrix, reduced, annotateParent = TRUE, annotationLabel = "parentTerm", fontsize=8)
  treemapPlot(reduced)
  wordcloudPlot(reduced, min.freq = 1, colors = "black")
  
  #
  # create an output file
  #
  out  = terms %>% group_by(go_id, name_1006) %>% summarise(count = n(), .groups = "drop")
  write.table(out, out_file, sep = '\t', row.names = FALSE, quote = FALSE)
  full = sub(".tsv$", ".full.tsv", out_file)
  write.table(terms, full, sep = '\t', row.names = FALSE, quote = FALSE)
  
}

go_analysis("Group.Helpers-Marine.Expansions.tsv", "Marine.GOs.tsv")
go_analysis("Group.Helpers-Salt.Expansions.tsv", "Salinity.GOs.tsv")
go_analysis("Group.Helpers-Agig.Expansions.tsv", "Gigantism.GOs.tsv")
go_analysis("Group.Helpers-Desert.Expansions.tsv", "Desert.GOs.tsv")
go_analysis("Group.Helpers-Water.Expansions.tsv", "Water.GOs.tsv")
go_analysis("Group.Helpers-Land.Expansions.tsv", "Land.GOs.tsv")

















