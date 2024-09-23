# in this script, use co-expression networks as inputs to compute GO AUROCs using EGAD

# LOAD IN ALL THE PIECES OF DATA NEEDED
pdata <- list()
pdata$genes <- readRDS(paste0(workspace$outputDir, "genes_metadata.rds")) # genes

pdata$coexmatFiles <- tibble(outfile = list.files(workspace$outputDir)) %>% 
  filter(grepl("coexmats", outfile)) %>% 
  mutate(outfile_name = gsub(".rds", "", outfile)) %>% 
  mutate(strs = strsplit(outfile_name, "_")) %>% 
  mutate(level = strs %>% sapply(function(currStr) { currStr[1] })) %>% 
  mutate(type = strs %>% sapply(function(currStr) { currStr[2] })) %>% 
  mutate(dataset = strs %>% sapply(function(currStr) { currStr[3] })) %>% 
  mutate(cell_type = strs %>% sapply(function(currStr) { currStr[4] })) %>% 
  dplyr::select(outfile, level, type, dataset, cell_type)

# +++++++++++++++++++
# load in the GO annotation for use later

fdata <- list()

fdata$go <- read.delim(file = paste0(workspace$miscDir, "pro_GO.csv"), sep = ",", stringsAsFactors = FALSE) %>% as_tibble()
fdata$go <- fdata$go %>% dplyr::select(gene = DB_Object_Symbol, go_id = GO.ID) %>% unique()
fdata$go <- fdata$go %>% inner_join(pdata$genes %>% dplyr::select(gene_id, gene) %>% unique(), by = "gene") 
fdata$go <- fdata$go %>% dplyr::select(go_id, gene_id) %>% unique() %>% arrange(go_id, gene_id)

fdata$goCounts <- fdata$go %>% 
  group_by(go_id) %>% summarize(n_gene = n()) %>% 
  arrange(desc(n_gene))

fdata$goValid <- fdata$goCounts %>% filter(n_gene >= 20) # filter for these once there so we never bother with them (because <20 are numerous)

# filter for valid go terms in go 
fdata$go <- fdata$go %>% filter(go_id %in% fdata$goValid$go_id)

# turn it into "list" so each geneset is easily accessible
fdata$goFlat <- fdata$go$go_id %>% unique()
names(fdata$goFlat) <- fdata$goFlat
fdata$goFlat <- fdata$goFlat %>% lapply(function(currGoId) { unique((fdata$go %>% filter(go_id == currGoId))$gene_id) })

# construct spread version for EGAD annotation format
fdata$goSpread <- fdata$go %>% mutate(count = 1) %>% spread(go_id, count)
fdata$goSpread <- fdata$goSpread %>% session$dataWrangler$setColAsRownames("gene_id")
fdata$goSpread <- fdata$goSpread %>% as.matrix()
fdata$goSpread[is.na(fdata$goSpread)] <- 0

fdata$goSpread %>% dim()

fdata$goLabels <- select(GO.db, columns = c("GOID","TERM"), keys = fdata$goValid$go_id, keytype = "GOID")
fdata$goLabels <- fdata$goLabels %>% as_tibble() %>% dplyr::select(go_id = GOID, go_label = TERM)

fdata$output <- list()
fdata$output$tbl <- fdata$go
fdata$output$flat <- fdata$goFlat
fdata$output$mat <- fdata$goSpread

fdata$output$labels <- fdata$goLabels

# =================================  
# COMMIT ==========================
# fdata$output %>% saveRDS(paste0(workspace$outputDir, "gene_ontology.rds"))
pdata$go <- readRDS(paste0(workspace$outputDir, "gene_ontology.rds"))
# =================================


# SANDBOX HERE...


fdata <- list()

fdata$coexmat <- pdata$coexmatFiles$outfile[5] 
fdata$coexmat <- readRDS(paste0(workspace$outputDir, fdata$coexmat))

fdata$genes <- rownames(fdata$coexmat) %>% intersect(rownames(pdata$go))
fdata$go <- pdata$go[fdata$genes, ]

fdata$goCounts <- fdata$go %>% apply(2, sum) %>% session$dataWrangler$vectorToTibble() %>% dplyr::select(go_id = variable, n_gene = value)
fdata$goValid <- fdata$goCounts %>% filter(n_gene >= 20) # filter for these again at the time of computation so we get the final valid set we want, per cell type
fdata$go <- fdata$go[, fdata$goValid$go_id]

fdata$results <- neighbor_voting(genes.labels = fdata$go,
                                 network = fdata$coexmat,
                                 nFold = 3, output = "AUROC")

fdata$results <- fdata$results %>% session$dataWrangler$setRownameAsColumn("go_id")

fdata$results %>% left_join(fdata$goLabs, by = "go_id") -> fdata$results

fdata$results %>% arrange(desc(auc)) %>% View()


