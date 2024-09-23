pdata <- list()

# +++++++++++++++++++
# load in the genes annotation for use later

fdata <- list()

fdata$genes <- read.table(paste0(workspace$miscDir, "Generic_human_ensemblIds_bioProcess.an.txt"), 
                          skip = 8, sep = "\t", header = TRUE, stringsAsFactors = FALSE, 
                          quote = "") %>% as_tibble() %>% 
  dplyr::select(ensembl = ProbeName, 
                gene = GeneSymbols, 
                entrez = NCBIids, 
                gene_name = GeneNames)

fdata$genes <- fdata$genes %>% filter(gene != "", entrez != "", gene_name != "")

fdata$genes <- fdata$genes %>% group_by(ensembl) %>% mutate(n = n()) %>% ungroup() %>% filter(n == 1) %>% dplyr::select(-n)

fdata$genes <- fdata$genes %>% mutate(gene_id = ensembl) %>% dplyr::select(gene_id, everything())

# =================================  
# COMMIT ==========================
# fdata$genes %>% saveRDS(paste0(workspace$outputDir, "genes_metadata.rds"))
pdata$genes <- read_rds(paste0(workspace$outputDir, "genes_metadata.rds"))
# =================================


# +++++++++++++++++++
# load in marker genes lists

fdata <- list()

fdata$mkrs <- read_rds(paste0(workspace$miscDir, "mancarci_mouseMarkerGenesNCBI.rds"))$Cortex

fdata$mkrs <- fdata$mkrs %>% 
  session$collectionUtils$lapplyWithName(function(cluster, currMkrs) { tibble(cluster = cluster, gene = currMkrs) }) %>% 
  session$dataWrangler$rbind()

fdata$mkrs <- fdata$mkrs %>% mutate(cell_type = cluster %>% sapply(function(currCluster) {
  if (grepl("^Astrocyte", currCluster)) { "astrocyte" }
  else if (grepl("^Endothelial", currCluster)) { "endothelial" }
  else if (grepl("^Gaba", currCluster)) { "inhibitory" } 
  else if (grepl("^Microglia$", currCluster)) { "microglia" } 
  else if (grepl("^Oligo$", currCluster)) { "oligodendrocyte" } 
  else if (grepl("^OligoPrecursors", currCluster)) { "opc" } 
  else if (grepl("^Pyramidal", currCluster)) { "excitatory" } 
  else { NA }
})) %>% na.omit() %>% unique()

fdata$ortholog <- mouse2human(fdata$mkrs$gene) %>% 
  as_tibble() %>% 
  dplyr::select(gene_mouse = mouseGene, gene_human = humanGene)

# map over to human genes
fdata$mkrs <- fdata$mkrs %>% 
  dplyr::select(cell_type, gene_mouse = gene) %>% 
  left_join(fdata$ortholog, by = "gene_mouse") %>% 
  dplyr::select(cell_type, gene = gene_human) %>% 
  na.omit() %>% unique()

# convert everything to ensembl ids
fdata$mkrs <- fdata$mkrs %>% 
  left_join(pdata$genes %>% dplyr::select(gene, gene_id)) %>% na.omit() %>% unique()

fdata$mkrs %>% group_by(cell_type) %>% summarize(n = n()) # final tally

# produce a flat version 
fdata$celltypes <- fdata$mkrs$cell_type %>% unique() %>% sort()
names(fdata$celltypes) <- fdata$celltypes 

fdata$mkrsFlat <- fdata$celltypes %>% lapply(function(currCelltype) {
  fdata$mkrs %>% filter(cell_type == currCelltype) %>% session$dataWrangler$extractColumn("gene_id")
})

# also a matrix form? maybe later.. if needed

fdata$output$tbl <- fdata$mkrs
fdata$output$flat <- fdata$mkrsFlat

# =================================  
# COMMIT ==========================
fdata$output %>% saveRDS(paste0(workspace$outputDir, "mkrs_mancarci.rds"))
# =================================

























