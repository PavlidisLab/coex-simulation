# in this script, load in the raw single cell data
# perform QC, normalization, so that the expression data are ready to be analyzed

# single cell data dir: /cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets

# LOAD IN ALL THE PIECES OF DATA NEEDED
pdata <- list() # ; none in this script
pdata$genes <- read_rds(paste0(workspace$outputDir, "genes_metadata.rds"))

# how many unique ensembl genes are there?
pdata$genes$gene_id %>% length()
pdata$genes$gene_id %>% unique() %>% length()

# unique gene symbols? 
pdata$genes$gene %>% length()
pdata$genes$gene %>% unique() %>% length() # 100 or so are duplicates - I think I filter these out mostly. 

# unique entrez ids? 
pdata$genes$entrez %>% length()
pdata$genes$entrez %>% unique() %>% length() # 1-to-1 mapping between entrez and gene symbols? 

pdata$celltypes <- c("excitatory", "inhibitory", "opc", "oligodendrocyte", "astrocyte", "microglia")
names(pdata$celltypes) <- pdata$celltypes

# +++++++++++++++++++++
# load in Velmeshev single cell dataset
# 1. clean (filter for genes detected in at least 10% of all cells, filter samples 2 SD away), per cell type
# 2. normalize
# 3. get the expression matrices to a usable state

fdata <- list()

fdata$samples <- read.table(paste0(workspace$dataDir, "velmeshev/", "meta.txt"), 
                            header = TRUE,
                            sep = "\t", 
                            stringsAsFactors = FALSE) %>% as_tibble()

fdata$samples <- fdata$samples %>% mutate(cell_type = cluster %>% sapply(function(currCluster) { 
  if (grepl("^AST", currCluster)) { "astrocyte" } 
  else if (grepl("^Endothelial", currCluster)) { "endothelial" } 
  else if (grepl("^IN", currCluster)) { "inhibitory" } 
  else if (grepl("^L", currCluster)) { "excitatory" } 
  else if (grepl("^Microglia", currCluster)) { "microglia" } 
  else if (grepl("^Neu-NRGN", currCluster)) { "nrgn_neu" } 
  else if (grepl("^Neu-mat", currCluster)) { "mat_neu" } 
  else if (grepl("^Olig", currCluster)) { "oligodendrocyte" } 
  else if (grepl("^OPC", currCluster)) { "opc" } 
  else { NA }
})) %>% na.omit()

fdata$samples <- fdata$samples %>% 
  dplyr::select(sample_ = sample, sample = cell, everything()) %>% 
  dplyr::select(sample, everything())

fdata$samples <- fdata$samples %>% 
  mutate(sample_bk = paste0(individual, "_", region)) %>% 
  dplyr::select(sample, sample_bk, everything())

fdata$genes <- read.table(paste0(workspace$dataDir, "velmeshev/", "genes.tsv"), 
                          header = FALSE,
                          sep = "\t", 
                          stringsAsFactors = FALSE) %>% as_tibble() %>% 
  dplyr::select(gene_id = V1, gene = V2)


fdata$cells <- read.table(paste0(workspace$dataDir, "velmeshev/", "barcodes.tsv"), 
                          header = FALSE,
                          sep = "\t", 
                          stringsAsFactors = FALSE) %>% as_tibble() %>% 
  dplyr::select(cell = V1)

fdata$ctmat <- workspace$utils$readMatSparse(paste0(workspace$dataDir, "velmeshev/", "matrix.mtx")) 
rownames(fdata$ctmat) <- fdata$genes$gene_id
colnames(fdata$ctmat) <- fdata$cells$cell

# filter for genes only present in the gemma annotation file
fdata$genes <- pdata$genes %>% filter(gene_id %in% fdata$ctmat@Dimnames[[1]]) %>% arrange(gene_id)
fdata$ctmat <- fdata$ctmat[fdata$genes$gene_id, ]

# per cell type, filter for 1. genes that are detected in at least 5% of all cells, and 2. filter out cells in the top or bottom 5% N genes detected
fdata$celltypes <- pdata$celltypes %>% intersect(fdata$samples$cell_type) 
names(fdata$celltypes) <- fdata$celltypes

fdata$exprmatsV <- fdata$celltypes %>% mclapply(function(currCelltype) {
  
  currSamples <- fdata$samples %>% filter(cell_type == currCelltype)
  currCtmat <- fdata$ctmat[, currSamples$sample]
  currCtmatV <- currCtmat %>% workspace$utils$cleanCtmat()
  currExprmatV <- currCtmatV %>% workspace$utils$ctmatToExprmat()
  
  return(currExprmatV)
  
}, mc.cores = length(fdata$celltypes))

fdata$exprmatsV %>% lapply(function(item) { item@Dim })

fdata$samplesV <- fdata$exprmatsV %>% session$collectionUtils$lapply(function(dat) {
  fdata$samples %>% filter(sample %in% colnames(dat))
})

# dat - exprmat; genes; samples
# one "all", and then one per cell type

fdata$output <- list()

fdata$output$raw$ctmat <- fdata$ctmat
fdata$output$raw$samples <- fdata$samples %>% filter(sample %in% colnames(fdata$ctmat))

fdata$output$exprmats <- fdata$exprmatsV
fdata$output$samples <- fdata$samplesV

# =================================  
# COMMIT ==========================
fdata$output %>% saveRDS(paste0(workspace$outputDir, "sc_exprmats_velmeshev.rds"))
# =================================

# read depth? number of cells? number of cell types? cell type mappings? 
# pdata$exprdats$sc$velmeshev$exprmats$astrocyte[1:10, 1:10]
fdata <- list()

fdata$samples <- do.call(rbind, pdata$exprdats$sc$velmeshev$samples %>% 
                           lapply(function(currSamples) { currSamples %>% dplyr::select(sample, sample_bk, individual, age, sex, diagnosis, cell_type, cluster) }))

fdata$samples %>%
  group_by(sample_bk, cell_type) %>% 
  summarize(n_cell = n()) %>% 
  session$graphingUtils$ggplot(aes(x = sample_bk, y = n_cell, fill = cell_type)) +
  geom_bar(stat = "identity") + 
  session$graphingUtils$tiltX(90)

fdata$samples %>%
  group_by(cell_type, cluster) %>% 
  summarize(n_cell = n()) %>% 
  session$graphingUtils$ggplot(aes(x = cluster, y = n_cell, fill = cell_type)) +
  geom_bar(stat = "identity") + 
  session$graphingUtils$tiltX(90)


# +++++++++++++++++++++
# compute subject level expression data by aggregating cells per subject * cell type

# Velemesehv 

fdata <- list()

fdata$exprdat_sc <- readRDS(paste0(workspace$outputDir, "sc_exprmats_velmeshev.rds"))

fdata$celltypes <- names(fdata$exprdat_sc$samples) %>% unique() %>% sort()
names(fdata$celltypes) <- fdata$celltypes

fdata$samples <- fdata$celltypes %>% lapply(function(currCelltype) {
  fdata$exprdat_sc$samples[[currCelltype]] %>%
    dplyr::select(sample, sample_bk) %>% 
    mutate(cell_type = currCelltype)
})

fdata$exprmats <- fdata$celltypes %>% mclapply(function(currCelltype) {
  currSamples <- fdata$samples[[currCelltype]]
  currExprmat <- fdata$exprdat_sc$exprmats[[currCelltype]]
  currExprmatSbj <- workspace$utils$aggregSbj(currExprmat, currSamples)
  return(currExprmatSbj)
}, mc.cores = length(fdata$celltypes))

fdata$exprmats %>% lapply(dim)

# visualize to see whether quantile normalization is needed; appears so
fdata$exprmats$inhibitory %>%
  session$dataWrangler$setRownameAsColumn("gene_id") %>% 
  gather(gene_id, subject) %>% 
  session$graphingUtils$ggplot(aes(x = subject, y = gene_id)) + 
  geom_boxplot()

fdata$exprmats <- fdata$exprmats %>% session$collectionUtils$lapply(workspace$utils$logTrans)
fdata$exprmats <- fdata$exprmats %>% session$collectionUtils$lapply(workspace$utils$quantileNormalize) # DECIDED TO QUANTILE NORMALIZE

fdata$samples <- do.call(rbind, fdata$samples) %>% 
  group_by(cell_type, sample_bk) %>% summarize(n_cell = n()) %>% ungroup() %>% 
  mutate(sample = sample_bk) %>% 
  dplyr::select(sample, everything())

fdata$output <- list()
fdata$output$samples <- fdata$samples
fdata$output$exprmats <- fdata$exprmats

# =================================  
# COMMIT ==========================
# fdata$output %>% saveRDS(paste0(workspace$outputDir, "sbj_exprmats_velmeshev.rds"))
# =================================

fdata <- list()

fdata$samples <- readRDS(paste0(workspace$outputDir, "sbj_exprmats_velmeshev.rds"))$samples

fdata$samples %>% 
  dplyr::select(sample, cell_type, n_cell) %>% 
  mutate(n_cell = log10(n_cell + 1)) %>% 
  spread(cell_type, n_cell) %>% 
  session$dataWrangler$fillNa(value = 0) %>%
  session$dataWrangler$setColAsRownames("sample") %>% 
  session$graphingUtils$heatmap(decimalNums = 1)


# +++++++++++++++++++++
# load in Rosmap single cell dataset
# 1. clean (filter for genes detected in at least 5% of all cells, filter samples 1.96 SD away), per cell type
# 2. normalize
# 3. get the expression matrices to a usable state

fdata <- list(); gc()

fdata$samples <- read.table(paste0(workspace$dataDir, "rosmap/", "filtered_column_metadata.txt"), 
                            header = TRUE,
                            sep = "\t", 
                            stringsAsFactors = FALSE) %>% as_tibble()

fdata$samples <- fdata$samples %>% mutate(cell_type = broad.cell.type %>% sapply(function(currCluster) { 
  if (grepl("^Ast", currCluster)) { "astrocyte" } 
  else if (grepl("^End", currCluster)) { "endothelial" } 
  else if (grepl("^Ex", currCluster)) { "excitatory" } 
  else if (grepl("^In", currCluster)) { "inhibitory" } 
  else if (grepl("^Mic", currCluster)) { "microglia" } 
  else if (grepl("^Oli", currCluster)) { "oligodendrocyte" } 
  else if (grepl("^Opc", currCluster)) { "opc" } 
  else if (grepl("^Per", currCluster)) { "per" } 
  else { NA }
})) %>% na.omit()

fdata$samples <- fdata$samples %>% 
  mutate(projid = as.character(projid)) %>% 
  dplyr::select(sample = TAG, everything()) 

fdata$genes <- read.table(paste0(workspace$dataDir, "rosmap/", "filtered_gene_row_names.txt"), 
                          header = FALSE,
                          sep = "\t", 
                          stringsAsFactors = FALSE) %>% as_tibble() %>% 
  dplyr::select(gene_id = V1)


fdata$cells <- read.table(paste0(workspace$dataDir, "rosmap/", "filtered_column_metadata.txt"), 
                          header = TRUE,
                          sep = "\t", 
                          stringsAsFactors = FALSE) %>% as_tibble() %>% 
  dplyr::select(cell = TAG)


fdata$ctmat <- workspace$utils$readMatSparse(paste0(workspace$dataDir, "rosmap/", "filtred_count_matrix.mtx")) 
rownames(fdata$ctmat) <- fdata$genes$gene_id
colnames(fdata$ctmat) <- fdata$cells$cell

# filter for genes only present in the gemma annotation file
fdata$genes <- pdata$genes %>% filter(gene %in% fdata$ctmat@Dimnames[[1]]) %>% arrange(gene_id)
fdata$genes <- fdata$genes %>% group_by(gene) %>% mutate(n = n()) %>% ungroup() %>% filter(n == 1) %>% dplyr::select(-n) # get rid of duplicates (one symbol mapping to multiple ensembl ids)
fdata$genes <- fdata$genes %>% arrange(gene_id)
fdata$ctmat <- fdata$ctmat[fdata$genes$gene, ] 
rownames(fdata$ctmat) <- fdata$genes$gene_id # map gene symbols to ensembl ids

# per cell type, filter for 1. genes that are detected in at least 5% of all cells, and 2. cells that top or bottom 5% of N genes detected
fdata$celltypes <- pdata$celltypes %>% intersect(fdata$samples$cell_type) 
names(fdata$celltypes) <- fdata$celltypes

fdata$exprmatsV <- fdata$celltypes %>% mclapply(function(currCelltype) {
  
  currSamples <- fdata$samples %>% filter(cell_type == currCelltype)
  currCtmat <- fdata$ctmat[, currSamples$sample]
  currCtmatV <- currCtmat %>% workspace$utils$cleanCtmat()
  currExprmatV <- currCtmatV %>% workspace$utils$ctmatToExprmat()
  
  return(currExprmatV)
  
}, mc.cores = length(fdata$celltypes))

fdata$exprmatsV %>% lapply(function(item) { item@Dim })

fdata$samplesV <- fdata$exprmatsV %>% session$collectionUtils$lapply(function(dat) {
  fdata$samples %>% filter(sample %in% colnames(dat))
})

# dat - exprmat; genes; samples
# one "all", and then one per cell type

fdata$output <- list()

fdata$output$raw$ctmat <- fdata$ctmat
fdata$output$raw$samples <- fdata$samples %>% filter(sample %in% colnames(fdata$ctmat))

fdata$output$exprmats <- fdata$exprmatsV
fdata$output$samples <- fdata$samplesV %>% lapply(function(currSamples) { currSamples %>% mutate(sample_bk = projid) })

# =================================  
# COMMIT ==========================
fdata$output %>% saveRDS(paste0(workspace$outputDir, "sc_exprmats_rosmap.rds"))
# =================================

fdata <- list()

fdata$samples <- do.call(rbind, readRDS((paste0(workspace$outputDir, "sc_exprmats_rosmap.rds")))$samples %>% 
                           lapply(function(currSamples) { currSamples %>% dplyr::select(sample, sample_bk = projid, cell_type, cluster = Subcluster, broad.cell.type) }))

fdata$samples %>%
  group_by(sample_bk, cell_type) %>% 
  summarize(n_cell = n()) %>% 
  session$graphingUtils$ggplot(aes(x = sample_bk, y = n_cell, fill = cell_type)) +
  geom_bar(stat = "identity") + 
  session$graphingUtils$tiltX(90)

fdata$samples %>%
  group_by(cell_type, cluster) %>% 
  summarize(n_cell = n()) %>% 
  session$graphingUtils$ggplot(aes(x = cluster, y = n_cell, fill = cell_type)) +
  geom_bar(stat = "identity") + 
  session$graphingUtils$tiltX(90)

# +++++++++++++++++++++
# compute subject level expression data by aggregating cells per subject * cell type

# Rosmap 

fdata <- list(); gc()

fdata$exprdat_sc <- readRDS(paste0(workspace$outputDir, "sc_exprmats_rosmap.rds"))

fdata$celltypes <- names(fdata$exprdat_sc$samples) %>% unique() %>% sort()
names(fdata$celltypes) <- fdata$celltypes

fdata$samples <- fdata$celltypes %>% lapply(function(currCelltype) {
  fdata$exprdat_sc$samples[[currCelltype]] %>%
    mutate(sample_bk = projid) %>% 
    dplyr::select(sample, sample_bk) %>% 
    mutate(cell_type = currCelltype)
})

fdata$exprmats <- fdata$celltypes %>% mclapply(function(currCelltype) {
  currSamples <- fdata$samples[[currCelltype]]
  currExprmat <- fdata$exprdat_sc$exprmats[[currCelltype]]
  currExprmatSbj <- workspace$utils$aggregSbj(currExprmat, currSamples)
  return(currExprmatSbj)
}, mc.cores = length(fdata$celltypes))

fdata$exprmats %>% lapply(function(item) { dim(item) })

# visualize to see whether quantile normalization is needed; appears so
fdata$exprmats$inhibitory %>%
  session$dataWrangler$setRownameAsColumn("gene_id") %>% 
  gather(gene_id, subject) %>% 
  session$graphingUtils$ggplot(aes(x = subject, y = gene_id)) + 
  geom_boxplot()

fdata$exprmats <- fdata$exprmats %>% session$collectionUtils$lapply(workspace$utils$logTrans)
fdata$exprmats <- fdata$exprmats %>% session$collectionUtils$lapply(workspace$utils$quantileNormalize) # DECIDED TO QUANTILE NORMALIZE

fdata$samples <- do.call(rbind, fdata$samples) %>% 
  group_by(cell_type, sample_bk) %>% summarize(n_cell = n()) %>% ungroup() %>% 
  mutate(sample = sample_bk) %>% 
  dplyr::select(sample, everything())

fdata$output <- list()
fdata$output$samples <- fdata$samples
fdata$output$exprmats <- fdata$exprmats

# =================================  
# COMMIT ==========================
fdata$output %>% saveRDS(paste0(workspace$outputDir, "sbj_exprmats_rosmap.rds"))
# =================================

fdata <- list()

fdata$samples <- readRDS(paste0(workspace$outputDir, "sbj_exprmats_rosmap.rds"))$samples

fdata$samples %>% 
  dplyr::select(sample, cell_type, n_cell) %>% 
  mutate(n_cell = log10(n_cell + 1)) %>% 
  spread(cell_type, n_cell) %>% 
  session$dataWrangler$fillNa(value = 0) %>%
  session$dataWrangler$setColAsRownames("sample") %>% 
  session$graphingUtils$heatmap(decimalNums = 1)

# +++++++++++++++++++++
# load in Velmeshev bulk tissue dataset
# 1. clean (filter for genes detected in at least 5% of all cells), per cell type
# 2. normalize the full table; quantile normalize this as well

fdata <- list()

fdata$samples <- read.table(paste0(workspace$dataDir, "velmeshev/", "ASD_BULK_meta.txt"),
                            header = TRUE,
                            sep = "\t",
                            stringsAsFactors = FALSE) %>% as_tibble() %>%
  dplyr::select(sample = Sample.name, everything()) %>%
  mutate(sample_bk = sample) %>%
  dplyr::select(sample, sample_bk, everything()) %>% 
  arrange(sample)

# load in bulk tissue ctmat
fdata$ctmat <- workspace$utils$readCounts(paste0(workspace$dataDir, "velmeshev/", "ASD_BULK.txt")) %>% 
  workspace$utils$asNumeric()

# filter for genes only present in the gemma annotation file
fdata$genes <- pdata$genes %>% filter(gene_id %in% rownames(fdata$ctmat)) %>% arrange(gene_id)
fdata$ctmat <- fdata$ctmat[fdata$genes$gene_id, fdata$samples$sample]

# clean it
fdata$ctmat <- fdata$ctmat %>% workspace$utils$cleanCtmat2(sampleThr = 0)

# normalize the whole expression matrix by library size + log transform
fdata$exprmat <- fdata$ctmat %>% workspace$utils$ctmatToExprmat2()

# visualize to see whether quantile normalization is needed; appears so
fdata$exprmat %>%
  session$dataWrangler$setRownameAsColumn("gene_id") %>% 
  gather(gene_id, subject) %>% 
  session$graphingUtils$ggplot(aes(x = subject, y = gene_id)) + 
  geom_boxplot()

# quantile normalize
fdata$exprmat <- fdata$exprmat %>% workspace$utils$logTrans()
fdata$exprmat <- fdata$exprmat %>% workspace$utils$quantileNormalize() # DECIDED TO QUANTILE NORMALIZE

fdata$output <- list()
fdata$output$samples <- fdata$samples
fdata$output$exprmat <- fdata$exprmat

# =================================  
# COMMIT ==========================
fdata$output %>% saveRDS(paste0(workspace$outputDir, "bk_exprmats_velmeshev.rds"))
# =================================

# +++++++++++++++++++++
# load in Rosmap bulk tissue dataset
# 1. clean (filter for genes detected in at least 5% of all cells, filter samples 2 SD away), per cell type
# 2. normalize the full table, missing fill in with 0 

fdata <- list()

fdata$samples <- (1 %>% lapply(function(i) {
  
  samples1 <- read.table(paste0(workspace$dataDir, "rosmap/", "ROSMAP_biospecimen_metadata.csv"), 
                         header = TRUE,
                         sep = ",", 
                         stringsAsFactors = FALSE) %>% as_tibble() %>% 
    dplyr::select(subject = individualID, specimen = specimenID, organ, tissue, tissueWeight, nucleicAcidSource, notes)
  
  samples2 <- read.table(paste0(workspace$dataDir, "rosmap/", "ROSMAP_Clinical_2019-05_v3.csv"), 
                         header = TRUE,
                         sep = ",", 
                         stringsAsFactors = FALSE) %>% as_tibble() %>% 
    mutate(projid = as.character(projid)) %>% 
    dplyr::select(subject = individualID, projid, everything())
  
  samples <- samples1 %>% inner_join(samples2, by = "subject") %>% 
    dplyr::select(subject, specimen, projid, everything())
  
  return(samples)
  
}))[[1]] %>% 
  dplyr::select(sample = specimen, subject, projid, everything()) %>% 
  arrange(sample)

fdata$ctmat <- (1 %>% lapply(function(i) {
  
  exprmat <- paste0(workspace$dataDir, "rosmap/", "ROSMAP_RNAseq_FPKM_gene.tsv") %>% fread(data.table = FALSE) %>% as.data.frame()
  rownames(exprmat) <- exprmat[, 2]
  exprmat <- exprmat[, 3:ncol(exprmat)]
  
  genes <- rownames(exprmat)
  genes <- genes %>% str_extract(".*(?=\\.)")
  rownames(exprmat) <- genes
  
  samples <- colnames(exprmat)
  samples <- samples %>% str_extract(".+_.+(?=_)")
  colnames(exprmat) <- samples
  
  return(exprmat)
  
}))[[1]] %>% as.matrix()

fdata$samples <- fdata$samples %>% filter(sample %in% colnames(fdata$ctmat)) %>% arrange(sample)

# filter for genes only present in the gemma annotation file
fdata$genes <- pdata$genes %>% filter(gene_id %in% rownames(fdata$ctmat)) %>% arrange(gene_id)
fdata$ctmat <- fdata$ctmat[fdata$genes$gene_id, fdata$samples$sample]

# clean it
fdata$ctmat <- fdata$ctmat %>% workspace$utils$cleanCtmat2(sampleThr = 0)

# normalize the whole expression matrix by library size
fdata$exprmat <- fdata$ctmat %>% workspace$utils$ctmatToExprmat2()

# visualize to see whether quantile normalization is needed; appears so
fdata$exprmat[, 1:30] %>%
  session$dataWrangler$setRownameAsColumn("gene_id") %>% 
  gather(gene_id, subject) %>% 
  session$graphingUtils$ggplot(aes(x = subject, y = gene_id)) + 
  geom_boxplot()

fdata$exprmat <- fdata$exprmat %>% workspace$utils$logTrans()
fdata$exprmat <- fdata$exprmat %>% workspace$utils$quantileNormalize() # DECIDED TO QUANTILE NORMALIZE

fdata$output <- list()
fdata$output$samples <- fdata$samples
fdata$output$exprmat <- fdata$exprmat

# =================================  
# COMMIT ==========================
fdata$output %>% saveRDS(paste0(workspace$outputDir, "bk_exprmats_rosmap.rds"))
# =================================

# +++++++++++++++++++++
# load in Rosmap bulk tissue dataset
# Subset bulk samples for those used to compute the xSubject networks

fdata <- list(); gc()

fdata$scSamples <- readRDS(paste0(workspace$outputDir, "sbj_exprmats_rosmap.rds"))$samples

fdata$bulkDat <- readRDS(paste0(workspace$outputDir, "bk_exprmats_rosmap.rds"))

fdata$samples <- fdata$bulkDat$samples %>% filter(projid %in% fdata$scSamples$sample) # turns out there are only 27 matched bulk tissue samples

fdata$exprmat <- fdata$bulkDat$exprmat[, fdata$samples$sample]

# visualize to see whether quantile normalization is needed; appears so
fdata$exprmat %>%
  session$dataWrangler$setRownameAsColumn("gene_id") %>% 
  gather(gene_id, subject) %>% 
  session$graphingUtils$ggplot(aes(x = subject, y = gene_id)) + 
  geom_boxplot()

fdata$output <- list()
fdata$output$samples <- fdata$samples
fdata$output$exprmat <- fdata$exprmat

# =================================  
# COMMIT ==========================
fdata$output %>% saveRDS(paste0(workspace$outputDir, "bk_exprmats_rosmap-scmatch.rds"))
# =================================


# +++++++++++++++++++++
# load in Nagy single cell dataset
# 1. clean (filter for genes detected in at least 10% of all cells, filter samples 5 percent top or bottom), per cell type
# 2. normalize
# 3. get the expression matrices to a usable state

fdata <- list()

fdata$samples <- read.table(paste0(workspace$dataDir, "nagy/", "GSE144136_CellNames.csv.gz"), 
                            header = TRUE,
                            sep = ",", 
                            stringsAsFactors = FALSE) %>% as_tibble() %>% na.omit() %>% 
  dplyr::select(i = X, sample = x) %>% 
  arrange(i)


fdata$genes <- read.table(paste0(workspace$dataDir, "nagy/", "GSE144136_GeneNames.csv.gz"), 
                          header = TRUE,
                          sep = ",", 
                          stringsAsFactors = FALSE) %>% as_tibble() %>% na.omit() %>% 
  dplyr::select(i = X, gene = x) %>% 
  arrange(i)

fdata$ctmat <- workspace$utils$readMatSparse(paste0(workspace$dataDir, "nagy/", "GSE144136_GeneBarcodeMatrix_Annotated.mtx.gz")) 
rownames(fdata$ctmat) <- fdata$genes$gene
colnames(fdata$ctmat) <- fdata$samples$sample

fdata$samples <- fdata$samples %>% 
  mutate(sample_strs = strsplit(sample, "\\.")) %>%
  mutate(cluster = sample_strs %>% sapply(function(str) { str[1] })) %>% 
  mutate(cell = sample_strs %>% sapply(function(str) { str[2] })) %>% 
  mutate(sample_bk = cell %>% sapply(function(currCell) { str_extract(currCell, "[0-9]+_.+_.+(?=_)")  })) %>% 
  dplyr::select(sample, sample_bk, cluster)

# now map over the cell types
fdata$samples <- fdata$samples %>% mutate(cell_type = cluster %>% sapply(function(currCluster) { 
  if (grepl("^Astros", currCluster)) { "astrocyte" } 
  else if (grepl("^Endo", currCluster)) { "endothelial" } 
  else if (grepl("^Inhib", currCluster)) { "inhibitory" } 
  else if (grepl("^Ex", currCluster)) { "excitatory" } 
  else if (grepl("^Micro", currCluster)) { "microglia" } 
  else if (grepl("^Mix", currCluster)) { "mix" } 
  else if (grepl("^Oligos", currCluster)) { "oligodendrocyte" } 
  else if (grepl("^OPC", currCluster)) { "opc" } 
  else { NA }
})) %>% na.omit()

fdata$samples %>% dplyr::select(cell_type, cluster) %>% unique() %>% arrange(cell_type, cluster) 

# filter for the samples with valid cell type annotations & ordering
fdata$samples <- fdata$samples %>% arrange(cell_type, sample_bk)
fdata$ctmat <- fdata$ctmat[, fdata$samples$sample] 

# filter for genes only present in the gemma annotation file
fdata$genes <- pdata$genes %>% filter(gene %in% fdata$ctmat@Dimnames[[1]]) %>% arrange(gene_id)
fdata$genes <- fdata$genes %>% group_by(gene) %>% mutate(n = n()) %>% ungroup() %>% filter(n == 1) %>% dplyr::select(-n) # get rid of duplicates (one symbol mapping to multiple ensembl ids)
fdata$genes <- fdata$genes %>% arrange(gene_id)
fdata$ctmat <- fdata$ctmat[fdata$genes$gene, ] 
rownames(fdata$ctmat) <- fdata$genes$gene_id # map gene symbols to ensembl ids

# per cell type, filter for 1. genes that are detected in at least 5% of all cells, and 2. filter out cells in the top or bottom 5% N genes detected
fdata$celltypes <- pdata$celltypes %>% intersect(fdata$samples$cell_type) 
names(fdata$celltypes) <- fdata$celltypes

fdata$exprmatsV <- fdata$celltypes %>% mclapply(function(currCelltype) {
  
  currSamples <- fdata$samples %>% filter(cell_type == currCelltype)
  currCtmat <- fdata$ctmat[, currSamples$sample]
  currCtmatV <- currCtmat %>% workspace$utils$cleanCtmat()
  currExprmatV <- currCtmatV %>% workspace$utils$ctmatToExprmat()
  
  return(currExprmatV)
  
}, mc.cores = length(fdata$celltypes))

fdata$exprmatsV %>% lapply(function(item) { item@Dim })

fdata$samplesV <- fdata$exprmatsV %>% session$collectionUtils$lapply(function(dat) {
  fdata$samples %>% filter(sample %in% colnames(dat))
})

# dat - exprmat; genes; samples
# one "all", and then one per cell type

fdata$output <- list()

fdata$output$raw$ctmat <- fdata$ctmat
fdata$output$raw$samples <- fdata$samples %>% filter(sample %in% colnames(fdata$ctmat))

fdata$output$exprmats <- fdata$exprmatsV
fdata$output$samples <- fdata$samplesV

# =================================  
# COMMIT ==========================
fdata$output %>% saveRDS(paste0(workspace$outputDir, "sc_exprmats_nagy.rds"))
# pdata$exprdats$sc$nagy <- readRDS(paste0(workspace$outputDir, "sc_exprmats_nagy.rds"))
# =================================

# read depth? number of cells? number of cell types? cell type mappings? 
fdata <- list()

fdata$samples <- do.call(rbind, readRDS(paste0(workspace$outputDir, "sc_exprmats_nagy.rds"))$samples %>% 
                           lapply(function(currSamples) { currSamples %>% dplyr::select(sample, sample_bk, cell_type, cluster) }))

fdata$samples %>%
  group_by(sample_bk, cell_type) %>% 
  summarize(n_cell = n()) %>% 
  session$graphingUtils$ggplot(aes(x = sample_bk, y = n_cell, fill = cell_type)) +
  geom_bar(stat = "identity") + 
  session$graphingUtils$tiltX(90)

fdata$samples %>%
  group_by(cell_type, cluster) %>% 
  summarize(n_cell = n()) %>% 
  session$graphingUtils$ggplot(aes(x = cluster, y = n_cell, fill = cell_type)) +
  geom_bar(stat = "identity") + 
  session$graphingUtils$tiltX(90)


# +++++++++++++++++++++
# compute subject level expression data by aggregating cells per subject * cell type

# Nagy 

fdata <- list()

fdata$exprdat_sc <- readRDS(paste0(workspace$outputDir, "sc_exprmats_nagy.rds"))

fdata$celltypes <- names(fdata$exprdat_sc$samples) %>% unique() %>% sort()
names(fdata$celltypes) <- fdata$celltypes

fdata$samples <- fdata$celltypes %>% lapply(function(currCelltype) {
  fdata$exprdat_sc$samples[[currCelltype]] %>%
    dplyr::select(sample, sample_bk) %>% 
    mutate(cell_type = currCelltype)
})

fdata$exprmats <- fdata$celltypes %>% session$collectionUtils$lapply(function(currCelltype) {
  currSamples <- fdata$samples[[currCelltype]]
  currExprmat <- fdata$exprdat_sc$exprmats[[currCelltype]]
  currExprmatSbj <- workspace$utils$aggregSbj(currExprmat, currSamples)
  return(currExprmatSbj)
})

fdata$exprmats %>% lapply(dim)

# visualize to see whether quantile normalization is needed; appears so
fdata$exprmats$inhibitory %>%
  session$dataWrangler$setRownameAsColumn("gene_id") %>% 
  gather(gene_id, subject) %>% 
  session$graphingUtils$ggplot(aes(x = subject, y = gene_id)) + 
  geom_boxplot()

fdata$exprmats <- fdata$exprmats %>% session$collectionUtils$lapply(workspace$utils$logTrans)
fdata$exprmats <- fdata$exprmats %>% session$collectionUtils$lapply(workspace$utils$quantileNormalize) # DECIDED TO QUANTILE NORMALIZE

fdata$samples <- do.call(rbind, fdata$samples) %>% 
  group_by(cell_type, sample_bk) %>% summarize(n_cell = n()) %>% ungroup() %>% 
  mutate(sample = sample_bk) %>% 
  dplyr::select(sample, everything())

fdata$output <- list()
fdata$output$samples <- fdata$samples
fdata$output$exprmats <- fdata$exprmats

# =================================  
# COMMIT ==========================
fdata$output %>% saveRDS(paste0(workspace$outputDir, "sbj_exprmats_nagy.rds"))
# =================================

# +++++++++++++++++++++
# load in Pineda single cell dataset
# 0. load in the data by merging multiple datasets!

fdata <- list()

fdata$dirs <- list.dirs(paste0(workspace$dataDir, "pineda/GSE174332_RAW"))[-1]

fdata$samples <- fdata$dirs %>% session$collectionUtils$lapply(function(currDir) {
  read.table(paste0(currDir, "/col_metadata.tsv"), 
             header = TRUE,
             sep = "\t", 
             stringsAsFactors = FALSE) %>% as_tibble() %>% na.omit()
}) %>% session$dataWrangler$rbind() %>% 
  mutate(sample = Barcode) %>% 
  dplyr::select(sample, everything())

fdata$ctmat <- fdata$dirs %>% session$collectionUtils$lapply(function(currDir) {
  
  genes <- read.table(paste0(currDir, "/row_metadata.tsv"),
                      header = TRUE,
                      sep = "\t",
                      stringsAsFactors = FALSE) %>% as_tibble() %>% na.omit()
  
  samples <- read.table(paste0(currDir, "/col_metadata.tsv"), 
             header = TRUE,
             sep = "\t", 
             stringsAsFactors = FALSE) %>% as_tibble() %>% na.omit()
  
  ctmat <- workspace$utils$readMatSparse(paste0(currDir, "/counts_fil.mtx")) 
  rownames(ctmat) <- genes$ENSEMBL
  colnames(ctmat) <- samples$Barcode
  
  return(ctmat)
}) %>% session$dataWrangler$cbind()


# dat - exprmat; genes; samples
# one "all", and then one per cell type

fdata$output <- list()

fdata$output$raw$ctmat <- fdata$ctmat
fdata$output$raw$samples <- fdata$samples %>% filter(sample %in% colnames(fdata$ctmat))

# =================================  
# COMMIT ==========================
fdata$output %>% saveRDS(paste0(workspace$outputDir, "sc_exprmats_pineda_RAW.rds"))
# =================================

# still working on the pineda dataset here
# 1. clean (filter for genes detected in at least 10% of all cells, filter samples 5 percent top or bottom), per cell type
# 2. normalize
# 3. get the expression matrices to a usable state

fdata <- list(); gc()

fdata$exprdat <- readRDS(paste0(workspace$outputDir, "sc_exprmats_pineda_RAW.rds"))

fdata$samples <- fdata$exprdat$raw$samples

fdata$samples <- fdata$samples %>% mutate(cell_type = CellType %>% sapply(function(currCluster) { 
  if (grepl("^Astrocyte", currCluster)) { "astrocyte" } 
  else if (grepl("^Endothelial", currCluster)) { "endothelial" } 
  else if (grepl("^Ex", currCluster)) { "excitatory" } 
  else if (grepl("^In", currCluster)) { "inhibitory" } 
  else if (grepl("^Fibroblast", currCluster)) { "fibroblast" } 
  else if (grepl("^Microglia", currCluster)) { "microglia" } 
  else if (grepl("^Mural", currCluster)) { "mural" } 
  else if (grepl("^Oligodendrocyte", currCluster)) { "oligodendrocyte" } 
  else if (grepl("^OPC", currCluster)) { "opc" } 
  else if (grepl("^T_Cell", currCluster)) { "tcell" } 
  else { NA }
})) %>% na.omit()

fdata$samples <- fdata$samples %>% 
  dplyr::select(sample, sample_bk = Sample_ID, cell_type, everything())

fdata$ctmat <- fdata$exprdat$raw$ctmat %>% as("dgTMatrix")

# filter for genes only present in the gemma annotation file
fdata$genes <- pdata$genes %>% filter(gene_id %in% fdata$ctmat@Dimnames[[1]]) %>% arrange(gene_id)
fdata$ctmat <- fdata$ctmat[fdata$genes$gene_id, ]

# per cell type, filter for 1. genes that are detected in at least 5% of all cells, and 2. filter out cells in the top or bottom 5% N genes detected
fdata$celltypes <- pdata$celltypes %>% intersect(fdata$samples$cell_type) 
names(fdata$celltypes) <- fdata$celltypes

fdata$exprmatsV <- fdata$celltypes %>% mclapply(function(currCelltype) {
  
  currSamples <- fdata$samples %>% filter(cell_type == currCelltype)
  currCtmat <- fdata$ctmat[, currSamples$sample]
  currCtmatV <- currCtmat %>% workspace$utils$cleanCtmat()
  currExprmatV <- currCtmatV %>% workspace$utils$ctmatToExprmat()
  
  return(currExprmatV)
  
}, mc.cores = length(fdata$celltypes))

fdata$exprmatsV %>% lapply(function(item) { item@Dim })

fdata$samplesV <- fdata$exprmatsV %>% session$collectionUtils$lapply(function(dat) {
  fdata$samples %>% filter(sample %in% colnames(dat))
})

# dat - exprmat; genes; samples
# one "all", and then one per cell type

fdata$output <- list()

fdata$output$raw$ctmat <- fdata$ctmat
fdata$output$raw$samples <- fdata$samples %>% filter(sample %in% colnames(fdata$ctmat))

fdata$output$exprmats <- fdata$exprmatsV
fdata$output$samples <- fdata$samplesV

# =================================  
# COMMIT ==========================
fdata$output %>% saveRDS(paste0(workspace$outputDir, "sc_exprmats_pineda.rds"))
# =================================

# +++++++++++++++++++++
# compute subject level expression data by aggregating cells per subject * cell type

# Pineda 

fdata <- list()

fdata$exprdat_sc <- readRDS(paste0(workspace$outputDir, "sc_exprmats_pineda.rds"))

fdata$celltypes <- names(fdata$exprdat_sc$samples) %>% unique() %>% sort()
names(fdata$celltypes) <- fdata$celltypes

fdata$samples <- fdata$celltypes %>% lapply(function(currCelltype) {
  fdata$exprdat_sc$samples[[currCelltype]] %>%
    dplyr::select(sample, sample_bk) %>% 
    mutate(cell_type = currCelltype)
})

fdata$exprmats <- fdata$celltypes %>% session$collectionUtils$lapply(function(currCelltype) {
  currSamples <- fdata$samples[[currCelltype]]
  currExprmat <- fdata$exprdat_sc$exprmats[[currCelltype]]
  currExprmatSbj <- workspace$utils$aggregSbj(currExprmat, currSamples)
  return(currExprmatSbj)
})

fdata$exprmats %>% lapply(dim)

# visualize to see whether quantile normalization is needed; appears so
fdata$exprmats$inhibitory %>%
  session$dataWrangler$setRownameAsColumn("gene_id") %>% 
  gather(gene_id, subject) %>% 
  session$graphingUtils$ggplot(aes(x = subject, y = gene_id)) + 
  geom_boxplot()


fdata$exprmats <- fdata$exprmats %>% session$collectionUtils$lapply(workspace$utils$logTrans)
fdata$exprmats <- fdata$exprmats %>% session$collectionUtils$lapply(workspace$utils$quantileNormalize) # DECIDED TO QUANTILE NORMALIZE

fdata$samples <- do.call(rbind, fdata$samples) %>% 
  group_by(cell_type, sample_bk) %>% summarize(n_cell = n()) %>% ungroup() %>% 
  mutate(sample = sample_bk) %>% 
  dplyr::select(sample, everything())

fdata$output <- list()
fdata$output$samples <- fdata$samples
fdata$output$exprmats <- fdata$exprmats

# =================================  
# COMMIT ==========================
fdata$output %>% saveRDS(paste0(workspace$outputDir, "sbj_exprmats_pineda.rds"))
# =================================

# +++++++++++++++++++++
# load in Lau single cell dataset
# 0. load in the data by merging multiple datasets!

fdata <- list(); gc()

fdata$dataDir <- paste0(workspace$dataDir, "lau/unzip/")

fdata$sampleIds <- list.files(fdata$dataDir) %>% 
  str_extract("GSM[0-9]+_[A-Z]+[0-9]+") %>% unique() %>% sort()
names(fdata$sampleIds) <- fdata$sampleIds 

fdata$samples <- fdata$sampleIds %>% session$collectionUtils$lapply(function(currId) {
  read.table(paste0(fdata$dataDir, currId, "_barcodes.tsv.gz"), 
             header = FALSE,
             sep = "\t", 
             stringsAsFactors = FALSE) %>% as_tibble() %>% na.omit() %>% 
    dplyr::select(sample = V1) %>% 
    mutate(sample = gsub("-", "_", sample)) %>% 
    mutate(sample_bk = currId) %>% 
    mutate(sample = paste0(sample_bk, "_", sample))
})

fdata$ctmat <- fdata$sampleIds %>% session$collectionUtils$lapply(function(currId) {
  
  genes <- read.table(paste0(fdata$dataDir, currId, "_features.tsv.gz"),
                      header = FALSE,
                      sep = "\t",
                      stringsAsFactors = FALSE) %>% as_tibble() %>% na.omit() %>% 
    dplyr::select(gene_id = V1, gene = V2)
  
  samples <- fdata$samples[[currId]]
  
  ctmat <- workspace$utils$readMatSparse(paste0(fdata$dataDir, currId, "_matrix.mtx.gz")) 
  rownames(ctmat) <- genes$gene_id
  colnames(ctmat) <- samples$sample
  
  return(ctmat)
}) %>% session$dataWrangler$cbind()

fdata$samples <- fdata$samples %>% session$dataWrangler$rbind()


# dat - exprmat; genes; samples
# one "all", and then one per cell type

fdata$output <- list()

fdata$output$raw$ctmat <- fdata$ctmat
fdata$output$raw$samples <- fdata$samples %>% filter(sample %in% colnames(fdata$ctmat))

# =================================  
# COMMIT ==========================
fdata$output %>% saveRDS(paste0(workspace$outputDir, "sc_exprmats_lau_RAW.rds"))
# =================================

# still working on the lau dataset here
# 1. clean (filter for genes detected in at least 10% of all cells, filter samples 5 percent top or bottom), per cell type
# 2. normalize
# 3. get the expression matrices to a usable state

fdata <- list(); gc()

fdata$exprdat <- readRDS(paste0(workspace$outputDir, "sc_exprmats_lau_RAW.rds"))

fdata$samples <- fdata$exprdat$raw$samples

fdata$samples <- fdata$samples %>% 
  mutate(barcode = sample %>% sapply(function(currSample) { str_extract(currSample, "(?<=_)[A-Z]{5,}") }))

fdata$samplesCT <- read.table(paste0(workspace$miscDir, "lau_celltypes.tsv"), header = TRUE) %>% 
  as_tibble() %>% 
  dplyr::select(barcode = ID, cell_type = Cell_type) %>% 
  mutate(barcode = barcode %>% sapply(function(currSample) { str_extract(currSample, "^[A-Z]{5,}") }))

fdata$samplesCT <- fdata$samplesCT %>% 
  unique() %>% group_by(barcode) %>% mutate(n = n()) %>% filter(n == 1) %>% ungroup() %>% dplyr::select(-n)

fdata$samples <- fdata$samples %>% right_join(fdata$samplesCT, by = "barcode")

fdata$samples <- fdata$samples %>% mutate(cell_type = cell_type %>% sapply(function(currCluster) { 
  if (grepl("^Astro", currCluster)) { "astrocyte" } 
  else if (grepl("^Endo", currCluster)) { "endothelial" } 
  else if (grepl("^Excit", currCluster)) { "excitatory" } 
  else if (grepl("^Inhit", currCluster)) { "inhibitory" } 
  else if (grepl("^Mic", currCluster)) { "microglia" } 
  else if (grepl("^Oligo", currCluster)) { "oligodendrocyte" } 
  else { NA }
})) %>% na.omit()

fdata$ctmat <- fdata$exprdat$raw$ctmat %>% as("dgTMatrix")

# filter for genes only present in the gemma annotation file
fdata$genes <- pdata$genes %>% filter(gene_id %in% fdata$ctmat@Dimnames[[1]]) %>% arrange(gene_id)
fdata$ctmat <- fdata$ctmat[fdata$genes$gene_id, ]

# per cell type, filter for 1. genes that are detected in at least 5% of all cells, and 2. filter out cells in the top or bottom 5% N genes detected
fdata$celltypes <- pdata$celltypes %>% intersect(fdata$samples$cell_type) 
names(fdata$celltypes) <- fdata$celltypes

fdata$exprmatsV <- fdata$celltypes %>% mclapply(function(currCelltype) {
  
  currSamples <- fdata$samples %>% filter(cell_type == currCelltype)
  currCtmat <- fdata$ctmat[, currSamples$sample]
  currCtmatV <- currCtmat %>% workspace$utils$cleanCtmat()
  currExprmatV <- currCtmatV %>% workspace$utils$ctmatToExprmat()
  
  return(currExprmatV)
  
}, mc.cores = length(fdata$celltypes))

fdata$exprmatsV %>% lapply(function(item) { item@Dim })

fdata$samplesV <- fdata$exprmatsV %>% session$collectionUtils$lapply(function(dat) {
  fdata$samples %>% filter(sample %in% colnames(dat))
})

# dat - exprmat; genes; samples
# one "all", and then one per cell type

fdata$output <- list()

fdata$output$raw$ctmat <- fdata$ctmat
fdata$output$raw$samples <- fdata$samples %>% filter(sample %in% colnames(fdata$ctmat))

fdata$output$exprmats <- fdata$exprmatsV
fdata$output$samples <- fdata$samplesV

# =================================
# COMMIT ==========================
fdata$output %>% saveRDS(paste0(workspace$outputDir, "sc_exprmats_lau.rds"))
# =================================

# +++++++++++++++++++++
# compute subject level expression data by aggregating cells per subject * cell type

# Lau 

fdata <- list()

fdata$exprdat_sc <- readRDS(paste0(workspace$outputDir, "sc_exprmats_lau.rds"))

fdata$celltypes <- names(fdata$exprdat_sc$samples) %>% unique() %>% sort()
names(fdata$celltypes) <- fdata$celltypes

fdata$samples <- fdata$celltypes %>% lapply(function(currCelltype) {
  fdata$exprdat_sc$samples[[currCelltype]] %>%
    dplyr::select(sample, sample_bk) %>% 
    mutate(cell_type = currCelltype)
})

fdata$exprmats <- fdata$celltypes %>% session$collectionUtils$lapply(function(currCelltype) {
  currSamples <- fdata$samples[[currCelltype]]
  currExprmat <- fdata$exprdat_sc$exprmats[[currCelltype]]
  currExprmatSbj <- workspace$utils$aggregSbj(currExprmat, currSamples)
  return(currExprmatSbj)
})

fdata$exprmats %>% lapply(dim)

# visualize to see whether quantile normalization is needed; appears so
fdata$exprmats$excitatory %>%
  session$dataWrangler$setRownameAsColumn("gene_id") %>% 
  gather(gene_id, subject) %>% 
  session$graphingUtils$ggplot(aes(x = subject, y = gene_id)) + 
  geom_boxplot()

fdata$exprmats <- fdata$exprmats %>% session$collectionUtils$lapply(workspace$utils$logTrans)
fdata$exprmats <- fdata$exprmats %>% session$collectionUtils$lapply(workspace$utils$quantileNormalize) # DECIDED TO QUANTILE NORMALIZE

fdata$samples <- do.call(rbind, fdata$samples) %>% 
  group_by(cell_type, sample_bk) %>% summarize(n_cell = n()) %>% ungroup() %>% 
  mutate(sample = sample_bk) %>% 
  dplyr::select(sample, everything())

fdata$output <- list()
fdata$output$samples <- fdata$samples
fdata$output$exprmats <- fdata$exprmats

# =================================  
# COMMIT ==========================
fdata$output %>% saveRDS(paste0(workspace$outputDir, "sbj_exprmats_lau.rds"))
# =================================

# +++++++++++++++++++++
# load in Lim single cell dataset
# 1. clean (filter for genes detected in at least 10% of all cells, filter samples 2 SD away), per cell type
# 2. normalize
# 3. get the expression matrices to a usable state

fdata <- list()

fdata$samples <- read.table(paste0(workspace$dataDir, "lim/", "GSE180928_metadata.csv.gz"), 
                            header = TRUE,
                            sep = ",", 
                            stringsAsFactors = FALSE) %>% as_tibble() %>% 
  dplyr::select(sample = X, sample_bk = case_num, everything()) %>% 
  mutate(sample = gsub("-", ".", sample)) 

fdata$ctmat <- paste0(workspace$dataDir, "lim/", "GSE180928_filtered_cell_counts.csv.gz") %>% fread(data.table = FALSE) %>% as.data.frame()
rownames(fdata$ctmat) <- fdata$ctmat[, "V1"]
fdata$ctmat <- fdata$ctmat[, -1]

fdata$samples <- fdata$samples %>% mutate(cell_type.x = Cluster %>% sapply(function(currCluster) { 
  if (grepl("^Astro", currCluster)) { "astrocyte" } 
  else if (grepl("^Endothelial", currCluster)) { "endothelial" } 
  else if (grepl("^Ependymal", currCluster)) { "Ependymal" } 
  else if (grepl("^GABA_", currCluster)) { "inhibitory" } 
  else if (grepl("^N_inh", currCluster)) { "inhibitory" } 
  else if (grepl("^Glut", currCluster)) { "excitatory" } 
  else if (grepl("^N_ex", currCluster)) { "excitatory" } 
  else if (grepl("^Microglia", currCluster)) { "microglia" } 
  else { NA }
})) %>% mutate(cell_type.y = Lineage %>% sapply(function(currCluster) {
  if (grepl("^Oligodendrocyte", currCluster)) { "oligodendrocyte" } 
  else if (grepl("^OPC", currCluster)) { "opc" }
  else { NA }
})) %>%
  session$dataWrangler$mergeColumnsXy("cell_type") %>%
  na.omit()

fdata$samples <- fdata$samples %>% arrange(cell_type)
fdata$ctmat <- fdata$ctmat[, fdata$samples$sample]
fdata$ctmat <- fdata$ctmat %>% as.matrix() %>% as("dgTMatrix")

# filter for genes only present in the gemma annotation file 
# filter out duplicates as well --- check ROSMAP
fdata$genes <- pdata$genes %>% filter(gene %in% fdata$ctmat@Dimnames[[1]]) %>% arrange(gene_id)
fdata$genes <- fdata$genes %>% group_by(gene) %>% mutate(n = n()) %>% ungroup() %>% filter(n == 1) %>% dplyr::select(-n) # get rid of duplicates (one symbol mapping to multiple ensembl ids)
fdata$genes <- fdata$genes %>% arrange(gene_id)
fdata$ctmat <- fdata$ctmat[fdata$genes$gene, ] 
rownames(fdata$ctmat) <- fdata$genes$gene_id # map gene symbols to ensembl ids

# per cell type, filter for 1. genes that are detected in at least 5% of all cells, and 2. filter out cells in the top or bottom 5% N genes detected
fdata$celltypes <- pdata$celltypes %>% intersect(fdata$samples$cell_type) 
names(fdata$celltypes) <- fdata$celltypes

fdata$exprmatsV <- fdata$celltypes %>% mclapply(function(currCelltype) {
  
  currSamples <- fdata$samples %>% filter(cell_type == currCelltype)
  currCtmat <- fdata$ctmat[, currSamples$sample]
  currCtmatV <- currCtmat %>% workspace$utils$cleanCtmat()
  currExprmatV <- currCtmatV %>% workspace$utils$ctmatToExprmat()
  
  return(currExprmatV)
  
}, mc.cores = length(fdata$celltypes))

fdata$exprmatsV %>% lapply(function(item) { item@Dim })

fdata$samplesV <- fdata$exprmatsV %>% session$collectionUtils$lapply(function(dat) {
  fdata$samples %>% filter(sample %in% colnames(dat))
})

# dat - exprmat; genes; samples
# one "all", and then one per cell type

fdata$output <- list()

fdata$output$raw$ctmat <- fdata$ctmat
fdata$output$raw$samples <- fdata$samples %>% filter(sample %in% colnames(fdata$ctmat))

fdata$output$exprmats <- fdata$exprmatsV
fdata$output$samples <- fdata$samplesV

# =================================  
# COMMIT ==========================
fdata$output %>% saveRDS(paste0(workspace$outputDir, "sc_exprmats_lim.rds"))
# pdata$exprdats$sc$lim <- readRDS(paste0(workspace$outputDir, "sc_exprmats_lim.rds"))
# =================================

# +++++++++++++++++++++
# compute subject level expression data by aggregating cells per subject * cell type

# Lim 

fdata <- list()

fdata$exprdat_sc <- readRDS(paste0(workspace$outputDir, "sc_exprmats_lim.rds"))

fdata$celltypes <- names(fdata$exprdat_sc$samples) %>% unique() %>% sort()
names(fdata$celltypes) <- fdata$celltypes

fdata$samples <- fdata$celltypes %>% lapply(function(currCelltype) {
  fdata$exprdat_sc$samples[[currCelltype]] %>%
    dplyr::select(sample, sample_bk) %>% 
    mutate(cell_type = currCelltype)
})

fdata$exprmats <- fdata$celltypes %>% session$collectionUtils$lapply(function(currCelltype) {
  currSamples <- fdata$samples[[currCelltype]]
  currExprmat <- fdata$exprdat_sc$exprmats[[currCelltype]]
  currExprmatSbj <- workspace$utils$aggregSbj(currExprmat, currSamples)
  return(currExprmatSbj)
})

fdata$exprmats %>% lapply(dim)

# visualize to see whether quantile normalization is needed; appears so
fdata$exprmats$excitatory %>%
  session$dataWrangler$setRownameAsColumn("gene_id") %>% 
  gather(gene_id, subject) %>% 
  session$graphingUtils$ggplot(aes(x = subject, y = gene_id)) + 
  geom_boxplot()

fdata$exprmats <- fdata$exprmats %>% session$collectionUtils$lapply(workspace$utils$logTrans)
fdata$exprmats <- fdata$exprmats %>% session$collectionUtils$lapply(workspace$utils$quantileNormalize) # DECIDED TO QUANTILE NORMALIZE

fdata$samples <- do.call(rbind, fdata$samples) %>% 
  group_by(cell_type, sample_bk) %>% summarize(n_cell = n()) %>% ungroup() %>% 
  mutate(sample = sample_bk) %>% 
  dplyr::select(sample, everything())

fdata$output <- list()
fdata$output$samples <- fdata$samples
fdata$output$exprmats <- fdata$exprmats

# =================================  
# COMMIT ==========================
fdata$output %>% saveRDS(paste0(workspace$outputDir, "sbj_exprmats_lim.rds"))
# pdata$exprdats$sbj$lim <- readRDS(paste0(workspace$outputDir, "sbj_exprmats_lim.rds"))
# =================================

# +++++++++++++++++++++
# load in Ramos single cell dataset
# 1. clean (filter for genes detected in at least 10% of all cells, filter samples 2 SD away), per cell type
# 2. normalize
# 3. get the expression matrices to a usable state

fdata <- list()

fdata$samples <- read.table(paste0(workspace$dataDir, "ramos/", "clean_metadata.csv"), 
                            header = TRUE,
                            sep = ",", 
                            stringsAsFactors = FALSE) %>% as_tibble() %>% 
  dplyr::select(sample = cellids, sample_bk = sample, celltypes)

fdata$samples <- fdata$samples %>%
  unique() %>% group_by(sample) %>% mutate(n = n()) %>% ungroup() %>% filter(n == 1) %>% dplyr::select(-n)

load(file = paste0(paste0(workspace$dataDir, "ramos/", "ctmat.rvar"))) # load as "ctmat"
fdata$ctmat <- ctmat %>% as("dgTMatrix")
rm(ctmat)

# filter for cells where we have both metadata and count data for
fdata$samples <- fdata$samples %>% filter(sample %in% colnames(fdata$ctmat))

# filter ctmat for the same samples
fdata$ctmat <- fdata$ctmat[, fdata$samples$sample]

fdata$samples <- fdata$samples %>% mutate(cell_type = celltypes %>% sapply(function(currCluster) { 
  if (grepl("^AC", currCluster)) { "astrocyte" } 
  else if (grepl("^EN", currCluster)) { "excitatory" } 
  else if (grepl("^IN", currCluster)) { "inhibitory" } 
  else if (grepl("^L[0-9]+", currCluster)) { "excitatory" } 
  else if (grepl("^MG", currCluster)) { "microglia" } 
  else if (grepl("^MSN", currCluster)) { "inhibitory" } 
  else if (grepl("^OL", currCluster)) { "oligodendrocyte" } 
  else if (grepl("^OPC", currCluster)) { "opc" } 
  else if (grepl("^preOL", currCluster)) { "oligodendrocyte" } 
  else if (grepl("^SPN", currCluster)) { "inhibitory" } 
  else { NA }
})) %>% na.omit()

fdata$samples <- fdata$samples %>% arrange(cell_type)
fdata$ctmat <- fdata$ctmat[, fdata$samples$sample]

# filter for genes only present in the gemma annotation file 
# filter out duplicates as well --- check ROSMAP
fdata$genes <- pdata$genes %>% filter(gene %in% fdata$ctmat@Dimnames[[1]]) %>% arrange(gene_id)
fdata$genes <- fdata$genes %>% group_by(gene) %>% mutate(n = n()) %>% ungroup() %>% filter(n == 1) %>% dplyr::select(-n) # get rid of duplicates (one symbol mapping to multiple ensembl ids)
fdata$genes <- fdata$genes %>% arrange(gene_id)
fdata$ctmat <- fdata$ctmat[fdata$genes$gene, ] 
rownames(fdata$ctmat) <- fdata$genes$gene_id # map gene symbols to ensembl ids

# per cell type, filter for 1. genes that are detected in at least 5% of all cells, and 2. filter out cells in the top or bottom 5% N genes detected
fdata$celltypes <- pdata$celltypes %>% intersect(fdata$samples$cell_type) 
names(fdata$celltypes) <- fdata$celltypes

fdata$exprmatsV <- fdata$celltypes %>% mclapply(function(currCelltype) {
  
  currSamples <- fdata$samples %>% filter(cell_type == currCelltype)
  currCtmat <- fdata$ctmat[, currSamples$sample]
  currCtmatV <- currCtmat %>% workspace$utils$cleanCtmat()
  currExprmatV <- currCtmatV %>% workspace$utils$ctmatToExprmat()
  
  return(currExprmatV)
  
}, mc.cores = length(fdata$celltypes))

fdata$exprmatsV %>% lapply(function(item) { item@Dim })

fdata$samplesV <- fdata$exprmatsV %>% session$collectionUtils$lapply(function(dat) {
  fdata$samples %>% filter(sample %in% colnames(dat))
})

# dat - exprmat; genes; samples
# one "all", and then one per cell type

fdata$output <- list()

fdata$output$raw$ctmat <- fdata$ctmat
fdata$output$raw$samples <- fdata$samples %>% filter(sample %in% colnames(fdata$ctmat))

fdata$output$exprmats <- fdata$exprmatsV
fdata$output$samples <- fdata$samplesV

# =================================  
# COMMIT ==========================
fdata$output %>% saveRDS(paste0(workspace$outputDir, "sc_exprmats_ramos.rds"))
# pdata$exprdats$sc$ramos <- readRDS(paste0(workspace$outputDir, "sc_exprmats_ramos.rds"))
# =================================

# +++++++++++++++++++++
# compute subject level expression data by aggregating cells per subject * cell type
# Ramos

fdata <- list()

fdata$exprdat_sc <- readRDS(paste0(workspace$outputDir, "sc_exprmats_ramos.rds"))

fdata$celltypes <- names(fdata$exprdat_sc$samples) %>% unique() %>% sort()
names(fdata$celltypes) <- fdata$celltypes

fdata$samples <- fdata$celltypes %>% lapply(function(currCelltype) {
  fdata$exprdat_sc$samples[[currCelltype]] %>%
    dplyr::select(sample, sample_bk) %>% 
    mutate(cell_type = currCelltype)
})

fdata$exprmats <- fdata$celltypes %>% session$collectionUtils$lapply(function(currCelltype) {
  currSamples <- fdata$samples[[currCelltype]]
  currExprmat <- fdata$exprdat_sc$exprmats[[currCelltype]]
  currExprmatSbj <- workspace$utils$aggregSbj(currExprmat, currSamples)
  return(currExprmatSbj)
})

fdata$exprmats %>% lapply(dim)

# visualize to see whether quantile normalization is needed; appears so
fdata$exprmats$inhibitory %>%
  session$dataWrangler$setRownameAsColumn("gene_id") %>% 
  gather(gene_id, subject) %>% 
  session$graphingUtils$ggplot(aes(x = subject, y = gene_id)) + 
  geom_boxplot()

fdata$exprmats <- fdata$exprmats %>% session$collectionUtils$lapply(workspace$utils$logTrans)
fdata$exprmats <- fdata$exprmats %>% session$collectionUtils$lapply(workspace$utils$quantileNormalize) # DECIDED TO QUANTILE NORMALIZE

fdata$samples <- do.call(rbind, fdata$samples) %>% 
  group_by(cell_type, sample_bk) %>% summarize(n_cell = n()) %>% ungroup() %>% 
  mutate(sample = sample_bk) %>% 
  dplyr::select(sample, everything())

fdata$output <- list()
fdata$output$samples <- fdata$samples
fdata$output$exprmats <- fdata$exprmats

# =================================  
# COMMIT ==========================
fdata$output %>% saveRDS(paste0(workspace$outputDir, "sbj_exprmats_ramos.rds"))
# pdata$exprdats$sbj$ramos <- readRDS(paste0(workspace$outputDir, "sbj_exprmats_ramos.rds"))
# =================================

pdata$exprmatFiles <- tibble(outfile = list.files(workspace$outputDir)) %>% 
  filter(grepl("exprmats", outfile)) %>% 
  mutate(outfile_name = gsub(".rds", "", outfile)) %>% 
  filter(!grepl("RAW", outfile_name)) %>% 
  mutate(strs = strsplit(outfile_name, "_")) %>% 
  mutate(level = strs %>% sapply(function(currStr) { currStr[1] })) %>% 
  mutate(dataset = strs %>% sapply(function(currStr) { currStr[3] })) %>% 
  dplyr::select(outfile, level, dataset)

# +++++++++++++++++++++++++++++++++++++++++++=
# work out the cel type expression profiles
fdata <- list()

fdata$datasets <- grep("^sc_exprmats_[a-z]+.rds$", list.files(workspace$outputDir), value = TRUE) %>% sort()
names(fdata$datasets) <- fdata$datasets %>% str_extract("[a-z]+(?=.rds)")

# loop through each dataset, compute one expression vector per cell type

fdata$ctps <- fdata$datasets %>% mclapply(function(currDataset) {
  
  dat <- readRDS(paste0(workspace$outputDir, currDataset))$exprmats
  
  celltypes <- pdata$celltypes %>% intersect(names(dat))
  names(celltypes) <- celltypes
  
  dat <- dat[celltypes]
  ctp <- dat %>% lapply(function(currExprmat) { rowMeans(currExprmat) })
  
  ctp %>% session$collectionUtils$lapplyWithName(function(currCelltype, currCtp) {
    currCtp %>% 
      session$dataWrangler$vectorToTibble() %>% 
      dplyr::select(gene_id = variable, expr = value) %>% 
      mutate(cell_type = currCelltype)
  }) %>% session$dataWrangler$rbind()
  
}, mc.cores = length(fdata$datasets)) %>% 
  session$collectionUtils$lapplyWithName(function(currName, currDat) {
    currDat %>% mutate(dataset = currName)
  }) %>% session$dataWrangler$rbind()

# =================================  
# COMMIT ==========================
fdata$ctps %>% saveRDS(paste0(workspace$outputDir, "sc_stats_ctprofiles.rds"))
pdata$ctprofiles <- read_rds(paste0(workspace$outputDir, "sc_stats_ctprofiles.rds"))
# =================================

# +++++++++++++++++++++
# load in Rosmap bulk tissue dataset
# Subset bulk samples for those with ihc
# fit IHC model and regress it out

fdata <- list(); gc()

fdata$ihc <- paste0(workspace$miscDir, "patrick_ihc") %>% list.files(full.names = TRUE) 
names(fdata$ihc) <- c("astrocyte", "endothelial", "microglia", "neuron", "oligodendrocyte")

fdata$ihc <- fdata$ihc %>% session$collectionUtils$lapplyWithName(function(currCelltype, currFile) {
  dat <- read.table(currFile, header = TRUE) %>% t() %>% as_tibble(rownames = "projid") 
  names(dat) <- c("projid", currCelltype)
  return(dat)
}) %>% reduce(full_join, by = "projid")
fdata$ihc <- fdata$ihc %>% mutate(projid = projid %>% str_extract("[0-9]+$"))
fdata$ihc <- fdata$ihc %>% na.omit() # just 49 samples where the proportions for all 5 cell types are known

fdata$ihc <- fdata$ihc %>% 
  session$dataWrangler$setColAsRownames("projid") %>% as.matrix() %>% 
  apply(1, function(vals) { vals / sum(vals) }) %>% t() %>%  # normalize 
  session$dataWrangler$setRownameAsColumn("projid")

# plot the ihc data across samples to get an idea of what the data looks like
fdata$ihc %>% 
  gather(cell_type, fraction, -projid) %>% 
  session$graphingUtils$ggplot(aes(x = projid, y = fraction)) + 
  geom_bar(stat = "identity", aes(fill = cell_type)) + 
  session$graphingUtils$tiltX(angle = 90)

fdata$bulkDat <- readRDS(paste0(workspace$outputDir, "bk_exprmats_rosmap.rds"))
fdata$samples <- fdata$bulkDat$samples %>% filter(projid %in% fdata$ihc$projid) %>% dplyr::select(sample, projid)
fdata$ihcMat <- fdata$samples %>% inner_join(fdata$ihc, by = "projid") %>% arrange(sample)
fdata$ihcMat <- fdata$ihcMat %>% dplyr::select(-projid)
fdata$ihcMat %>% dplyr::select(-sample) %>% GGally::ggpairs() # take a look at inter cell type proportions correlations
fdata$ihcMat <- fdata$ihcMat %>% session$dataWrangler$setColAsRownames("sample") %>% as.matrix()

fdata$exprmat <- fdata$bulkDat$exprmat[, rownames(fdata$ihcMat)]
fdata$ccvModel <- workspace$utils$fitCCVModels(fdata$ihcMat, fdata$exprmat)

# visualize to see whether quantile normalization is needed; appears so
fdata$ccvModel$exprmats$residual %>%
  session$dataWrangler$setRownameAsColumn("gene_id") %>% 
  gather(gene_id, subject) %>% 
  session$graphingUtils$ggplot(aes(x = subject, y = gene_id)) + 
  geom_boxplot()

fdata$output1 <- list()
fdata$output1$samples <- fdata$bulkDat$samples %>% filter(sample %in% colnames(fdata$ccvModel$exprmats$orig))
fdata$output1$exprmat <- fdata$ccvModel$exprmats$orig

fdata$output2 <- list()
fdata$output2$samples <- fdata$bulkDat$samples %>% filter(sample %in% colnames(fdata$ccvModel$exprmats$residual))
fdata$output2$exprmat <- fdata$ccvModel$exprmats$residual
fdata$output2$ccvModel <- fdata$ccvModel

# =================================  
# COMMIT ==========================
fdata$output1 %>% saveRDS(paste0(workspace$outputDir, "bk_exprmats_rosmap-ihc.rds"))
fdata$output2 %>% saveRDS(paste0(workspace$outputDir, "bk_exprmats_rosmap-ihcres.rds"))
# =================================


# +++++++++++++++++++++++++++++++
# for rosmap and velmeshev, construct new single cell and subject level exprmats using the union of genes detected in all cell types (call them uexprmats)
# use uexprmats for 1. idnetifying marker genes and 2. fitting the dilution model
# identify marker genes for neurons, not excitatory & inhibitory neurons by collapsing all neuronal cells for computing cell type profiles

# identify the union of genes & samples for each dataset

fdata <- list(); gc()

fdata$dataset <- "velmeshev"

fdata$scdat <- read_rds(paste0(workspace$outputDir, "sc_exprmats_", fdata$dataset,".rds"))

fdata$celltypes <- pdata$celltypes %>% intersect(names(fdata$scdat$exprmats)) %>% session$dataWrangler$attachNames()

fdata$genesV <- fdata$celltypes %>% # identify the union of valid genes across all cell type
  lapply(function(celltype) { fdata$scdat$exprmats[[celltype]] %>% rownames() }) %>% 
  unname() %>% unlist() %>% unique() %>% sort()

fdata$samplesV <- fdata$celltypes %>% 
  lapply(function(celltype) { fdata$scdat$exprmats[[celltype]] %>% colnames() }) %>% 
  unname() %>% unlist() %>% unique() %>% sort()

fdata$ctmatV <- fdata$scdat$raw$ctmat[fdata$genesV, fdata$samplesV]
fdata$exprmatV <- fdata$ctmatV %>% workspace$utils$ctmatToExprmat()

# now break them up into different cell types
fdata$exprmatV <- fdata$celltypes %>% mclapply(function(celltype) {
  samples <- fdata$scdat$samples[[celltype]]$sample
  fdata$exprmatV[, samples]
}, mc.cores = length(fdata$celltypes))

fdata$samplesV <- fdata$scdat$samples

# compute cell type profiles

fdata$ctp <- fdata$exprmatV %>% sapply(function(currExprmat) { rowMeans(currExprmat) })

fdata$ctp[rownames(fdata$ctp) != "ENSG00000251562", ] %>% as.data.frame() %>% GGally::ggpairs() # strangely there is a single gene that is extremely highly expressed
# correlation of expressions between inhibitory and excitatoryo neurons is about 0.8; still worth trying to get a set of markers for every cell type

fdata$ctpTbl <- fdata$ctp %>% 
  session$dataWrangler$setRownameAsColumn("gene_id") %>% 
  gather(cell_type, expr, -gene_id)

fdata$ctpTbl %>% 
  filter(gene_id != "ENSG00000251562") %>% 
  session$graphingUtils$ggplot(aes(x = log2(expr + 1))) + 
  geom_density(aes(color = cell_type))

# per cell type also compute subject level expression profiles

fdata$exprmatVSbj <- fdata$celltypes %>% mclapply(function(celltype) {
  workspace$utils$aggregSbj(fdata$exprmatV[[celltype]], samples = fdata$samplesV[[celltype]])
}, mc.cores = length(fdata$celltypes))


# =================================  
# COMMIT ==========================
fdata$output$sc$ctmat <- fdata$ctmatV
fdata$output$sc$exprmats <- fdata$exprmatV
fdata$output$sc$samples <- fdata$samplesV
fdata$output$sbj$exprmats <- fdata$exprmatVSbj
fdata$output$ctp$mat <- fdata$ctp
fdata$output$ctp$tbl <- fdata$ctpTbl
fdata$output %>% saveRDS(paste0(workspace$outputDir, "uexprdats_", fdata$dataset, ".rds"))
# =================================  

# construct the union genes dataset for rosmap as well

fdata <- list(); gc()

fdata$dataset <- "rosmap"

fdata$scdat <- read_rds(paste0(workspace$outputDir, "sc_exprmats_", fdata$dataset,".rds"))

fdata$celltypes <- pdata$celltypes %>% intersect(names(fdata$scdat$exprmats)) %>% session$dataWrangler$attachNames()

fdata$genesV <- fdata$celltypes %>% # identify the union of valid genes across all cell type
  lapply(function(celltype) { fdata$scdat$exprmats[[celltype]] %>% rownames() }) %>% 
  unname() %>% unlist() %>% unique() %>% sort()

fdata$samplesV <- fdata$celltypes %>% 
  lapply(function(celltype) { fdata$scdat$exprmats[[celltype]] %>% colnames() }) %>% 
  unname() %>% unlist() %>% unique() %>% sort()

fdata$ctmatV <- fdata$scdat$raw$ctmat[fdata$genesV, fdata$samplesV]
fdata$exprmatV <- fdata$ctmatV %>% workspace$utils$ctmatToExprmat()

# now break them up into different cell types
fdata$exprmatV <- fdata$celltypes %>% mclapply(function(celltype) {
  samples <- fdata$scdat$samples[[celltype]]$sample
  fdata$exprmatV[, samples]
}, mc.cores = length(fdata$celltypes))

fdata$samplesV <- fdata$scdat$samples

# compute cell type profiles

fdata$ctp <- fdata$exprmatV %>% sapply(function(currExprmat) { rowMeans(currExprmat) })

fdata$ctp %>% as.data.frame() %>% GGally::ggpairs() 
# correlation of expressions between inhibitory and excitatory neurons is about 0.8; still worth trying to get a set of markers for every cell type

fdata$ctpTbl <- fdata$ctp %>% 
  session$dataWrangler$setRownameAsColumn("gene_id") %>% 
  gather(cell_type, expr, -gene_id)

fdata$ctpTbl %>% 
  session$graphingUtils$ggplot(aes(x = log2(expr + 1))) + 
  geom_density(aes(color = cell_type))

# per cell type also compute subject level expression profiles

fdata$exprmatVSbj <- fdata$celltypes %>% mclapply(function(celltype) {
  workspace$utils$aggregSbj(fdata$exprmatV[[celltype]], samples = fdata$samplesV[[celltype]])
}, mc.cores = length(fdata$celltypes))


# =================================  
# COMMIT ==========================
fdata$output$sc$ctmat <- fdata$ctmatV
fdata$output$sc$exprmats <- fdata$exprmatV
fdata$output$sc$samples <- fdata$samplesV
fdata$output$sbj$exprmats <- fdata$exprmatVSbj
fdata$output$ctp$mat <- fdata$ctp
fdata$output$ctp$tbl <- fdata$ctpTbl
fdata$output %>% saveRDS(paste0(workspace$outputDir, "uexprdats_", fdata$dataset, ".rds"))
# =================================  

# +++++++++++++++++++++++++++++++++
# identify the marker genes (possibly intersect per cell type between the two datasets: rosmap & velmeshev)
# decided on top 100 minfc intersecting the two datasets

fdata <- list(); gc()

fdata$datasets <- c("velmeshev", "rosmap") %>% session$dataWrangler$attachNames()

fdata$ctp <- fdata$datasets %>% lapply(function(dataset) { read_rds(paste0(workspace$outputDir, "uexprdats_", dataset, ".rds"))$ctp })

# for every gene in every cell type and dataset, compute its expr, othr_expr_max, min_fc

fdata$minFc <- fdata$datasets %>% lapply(function(dataset) {
  pdata$celltypes %>% lapply(function(celltype) {
    
    ctp <- fdata$ctp[[dataset]]$mat
    targetExpr <- ctp[, celltype]
    othrExprMax <- ctp[, colnames(ctp) != celltype] %>% apply(1, max) 
    
    tibble(gene_id = names(targetExpr)) %>% 
      mutate(expr = targetExpr[gene_id], 
             othr_expr_max = othrExprMax[gene_id]) %>% 
      mutate(min_fc = expr / othr_expr_max) %>% 
      mutate(dataset = dataset, cell_type = celltype)
    
  }) %>% session$dataWrangler$rbind()
}) %>% session$dataWrangler$rbind()


# per dataset and cell type, take the top 100 genes by min_fc
fdata$minFc <- fdata$minFc %>% group_by(dataset, cell_type) %>% mutate(min_fc_rank = rank(desc(min_fc))) %>% ungroup() # compute the ranks

fdata$minFc %>% mutate(top100 = (min_fc_rank <= 100)) %>% 
  session$graphingUtils$ggplot(aes(x = min_fc + 1)) + 
  geom_density(aes(color = top100)) + 
  facet_wrap(~cell_type * dataset, ncol = 2) + 
  scale_x_continuous(trans = "log2")

fdata$mkrs <- pdata$celltypes %>% lapply(function(celltype) {
  fdata$datasets %>% lapply(function(currDataset) {
    fdata$minFc %>% filter(cell_type == celltype, dataset == currDataset, min_fc_rank <= 100) %>% session$dataWrangler$extractColumn("gene_id")
  })
})

fdata$mkrs <- fdata$mkrs %>% lapply(function(mkrs) {
  mkrs$rosmap %>% intersect(mkrs$velmeshev)
})

fdata$mkrs %>% sapply(length)

fdata$minFc %>% filter(gene_id %in% fdata$mkrs$excitatory, cell_type == "excitatory") %>% 
  dplyr::select(gene_id, dataset, min_fc) %>% 
  spread(dataset, min_fc) %>% 
  left_join(pdata$genes %>% dplyr::select(gene_id, gene)) %>% 
  arrange(desc(velmeshev))

fdata$mkrs %>% lapply(function(mkrs) { pdata$genes %>% filter(gene_id %in% mkrs) %>% session$dataWrangler$extractColumn("gene") %>% sort() })

# as a sanity check, let's plot the marker genes in each cell type using just the cell type expression profiles

fdata$ctp$velmeshev$mat[unlist(unname(fdata$mkrs)),] %>% 
  apply(1, function(vals) { (vals - mean(vals)) / sd(vals) }) %>% t() %>% 
  workspace$utils$logTrans() %>% 
  session$graphingUtils$heatmap(size = session$graphingUtils$LARGE, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE)

fdata$ctp$velmeshev$mat[unlist(unname(fdata$mkrs)),] %>% 
  workspace$utils$logTrans() %>% 
  session$graphingUtils$heatmap(size = session$graphingUtils$LARGE, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE)


fdata$ctp$rosmap$mat[unlist(unname(fdata$mkrs)),] %>% 
  apply(1, function(vals) { (vals - mean(vals)) / sd(vals) }) %>% t() %>% 
  workspace$utils$logTrans() %>% 
  session$graphingUtils$heatmap(size = session$graphingUtils$LARGE, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE)

fdata$ctp$rosmap$mat[unlist(unname(fdata$mkrs)),] %>% 
  workspace$utils$logTrans() %>% 
  session$graphingUtils$heatmap(size = session$graphingUtils$LARGE, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE)

# also make mkrsMat available for ease of use

fdata$mkrsMat <- fdata$mkrs %>% session$collectionUtils$lapplyWithName(function(celltype, mkrs) {
  tibble(cell_type = celltype, gene_id = mkrs)
}) %>% session$dataWrangler$rbind() %>% 
  mutate(in_set = 1) %>% 
  spread(cell_type, in_set) %>% 
  session$dataWrangler$fillNa(pdata$celltypes, value = 0) %>% 
  session$dataWrangler$setColAsRownames("gene_id")


# =================================  
# COMMIT ==========================
fdata$output$minfc <- fdata$minFc
fdata$output$mkrsFlat <- fdata$mkrs
fdata$output$mkrsMat <- fdata$mkrsMat
fdata$output %>% saveRDS(paste0(workspace$outputDir, "umkrs.rds"))
# =================================  

# +++++++++++++++++++++++++++++++++++++++++++++++++++
# next, use the identified marker genes to construct marker gene profiles in 
# 1. velmeshev

fdata <- list(); gc()

fdata$exprdat <- read_rds(paste0(workspace$outputDir, "bk_exprmats_velmeshev.rds"))

fdata$mkrs <- readRDS(paste0(workspace$outputDir, "umkrs.rds"))

fdata$mkrsV <- rownames(fdata$mkrs$mkrsMat) %>% intersect(rownames(fdata$exprdat$exprmat))

fdata$mkrsTbl <- fdata$mkrs$mkrsFlat %>% session$collectionUtils$lapplyWithName(function(celltype, mkrs) {
  tibble(cell_type = celltype, gene_id = mkrs) 
}) %>% session$dataWrangler$rbind() %>% 
  filter(gene_id %in% fdata$mkrsV) 

fdata$exprmatMkrs <- fdata$exprdat$exprmat[fdata$mkrsV, ]

xdata <- list()

xdata$anno <- fdata$mkrsTbl %>% session$dataWrangler$setColAsRownames("gene_id")

fdata$exprmatMkrs %>% 
  apply(1, function(vals) { (vals - mean(vals)) / sd(vals) }) %>% t() %>% 
  session$graphingUtils$heatmap(size = session$graphingUtils$SMALL, annotation_row = xdata$anno,
                                show_rownames = FALSE)

# plot out the pairwise co-expression distribution within and across cell types

xdata$coexmat <- fdata$exprmatMkrs %>% t() %>% cor()

xdata$coextbl <- xdata$coexmat %>% workspace$utils$coexMatToTbl() %>% 
  left_join(fdata$mkrsTbl %>% dplyr::select(gene_a = gene_id, cell_type_a = cell_type)) %>% 
  left_join(fdata$mkrsTbl %>% dplyr::select(gene_b = gene_id, cell_type_b = cell_type)) %>% 
  mutate(cell_type_in = cell_type_a == cell_type_b)
  
xdata$coextbl %>% 
  session$graphingUtils$ggplot(aes(x = cor_coef)) +
  geom_density(aes(color = cell_type_in)) + 
  facet_wrap(~cell_type_a)

# now, let's perform PCA to get the marker genes profiles
fdata$mgps <- workspace$utils$computeMgps(fdata$mkrs$mkrsFlat, fdata$exprmatMkrs)

# now, apply to the whole transcriptome, correct it out!
fdata$ccvModels <- workspace$utils$fitCCVModels(fdata$mgps$main, exprmat = fdata$exprdat$exprmat)

# save 
# mgp; ccv models; expression matrix + samples

# =================================  
# COMMIT ==========================
fdata$output$mgps <- fdata$mgps
fdata$output$ccvModels <- fdata$ccvModels
fdata$output$exprdat$exprmat <- fdata$ccvModels$exprmats$residual
fdata$output$exprdat$samples <- fdata$exprdat$samples
fdata$output %>% saveRDS(paste0(workspace$outputDir, "bk_exprmats_velmeshev-mgpres.rds"))
# =================================  

# TODO - TO BE DELETED
fdata$ccvModel <- readRDS(paste0(workspace$outputDir, "bk_exprmats_velmeshev-mgpres.rds"))$ccvModels

fdata$ccvModel$stats$lmStats %>% session$graphingUtils$ggplot(aes(x = rsqr)) + geom_histogram()
fdata$ccvModel$stats$lmStats$rsqr %>% mean(na.rm = TRUE)
fdata$ccvModel$stats$lmStats %>% session$graphingUtils$ggplot(aes(x = pvalue)) + geom_histogram()
fdata$ccvModel$stats$lmStats %>% filter(qvalue < 0.1)

# TODO - TO BE DELETED AS WELL
# SANDBOX - exploring the relationship between xSubject varaibility & CCV

fdata$mgpres <- readRDS(paste0(workspace$outputDir, "bk_exprmats_rosmap-mgpres.rds"))

fdata$ccvModel <- fdata$mgpres$ccvModels

fdata$main <- fdata$mgpres$mgps$rotations$excitatory[, "PC1"] %>% 
  session$dataWrangler$vectorToTibble() %>% 
  dplyr::select(gene_id = variable, loading = value)

fdata$exprdat$sc <- read_rds(paste0(workspace$outputDir, "uexprdats_rosmap.rds"))

fdata$exprmatSbj <- fdata$exprdat$sc$sbj$exprmats 

fdata$variability <- fdata$exprmatSbj %>% apply(1, sd) %>% 
  session$dataWrangler$vectorToTibble() %>% 
  dplyr::select(gene_id = variable, sd = value)

fdata$sd <- fdata$exprdat$sc$sbj$exprmats %>% sapply(function(exprmat) { exprmat %>% apply(1, sd) })

fdata$minFcSd <- pdata$celltypes %>% lapply(function(celltype) {
    
    targetExpr <- fdata$sd[, celltype]
    othrExprMax <- fdata$sd[, colnames(fdata$sd) != celltype] %>% apply(1, max) 
    
    tibble(gene_id = names(targetExpr)) %>% 
      mutate(expr = targetExpr[gene_id], 
             othr_expr_max = othrExprMax[gene_id]) %>% 
      mutate(min_fc = expr / othr_expr_max) %>% 
      mutate(cell_type = celltype)
    
}) %>% session$dataWrangler$rbind()

ct <- "oligodendrocyte"

fdata$mkrs$minfc %>% 
  filter(cell_type == ct, dataset == "rosmap") %>% 
  dplyr::select(gene_id, min_fc) %>% 
  left_join(fdata$minFcSd %>% filter(cell_type == ct) %>% dplyr::select(gene_id, min_fc_sd = min_fc), by = "gene_id") %>% 
  left_join(fdata$ccvModel$stats$coefStats %>% filter(cell_type == ct) %>% dplyr::select(gene_id, stat = rsqr_indep), by = "gene_id") %>% 
  session$graphingUtils$ggplot(aes(x = log2(min_fc + 1), y = log2(min_fc_sd + 1))) + 
  geom_point(aes(color = stat)) + 
  geom_abline(slope = 1, linetype = "dashed") 



fdata$main %>% left_join(fdata$variability) %>% session$graphingUtils$ggplot(aes(x = sd, y = loading)) + geom_point()
fdata$main %>% left_join(fdata$variability) -> x; cor(x$loading, x$sd) # no fucking way

fdata$ccvModel$stats$coefStats %>% filter(cell_type == "inhibitory") %>% dplyr::select(gene_id, rsqr_indep) %>% 
  inner_join(fdata$variability) %>% session$graphingUtils$ggplot(aes(x = sd, y = rsqr_indep)) + geom_point(shape = 1)

fdata$ccvModel$stats$coefStats %>% filter(cell_type == "inhibitory") %>% dplyr::select(gene_id, rsqr_indep) %>% 
  inner_join(fdata$variability) -> x; cor(x$rsqr_indep, x$sd, use = 'pairwise.complete.obs')

# next, use the identified marker genes to construct marker gene profiles in 
# 2. rosmap

fdata <- list(); gc()

fdata$exprdat <- read_rds(paste0(workspace$outputDir, "bk_exprmats_rosmap.rds"))

fdata$mkrs <- readRDS(paste0(workspace$outputDir, "umkrs.rds"))

fdata$mkrsV <- rownames(fdata$mkrs$mkrsMat) %>% intersect(rownames(fdata$exprdat$exprmat))

fdata$mkrsTbl <- fdata$mkrs$mkrsFlat %>% session$collectionUtils$lapplyWithName(function(celltype, mkrs) {
  tibble(cell_type = celltype, gene_id = mkrs) 
}) %>% session$dataWrangler$rbind() %>% 
  filter(gene_id %in% fdata$mkrsV) 

fdata$exprmatMkrs <- fdata$exprdat$exprmat[fdata$mkrsV, ]

xdata <- list()

xdata$anno <- fdata$mkrsTbl %>% session$dataWrangler$setColAsRownames("gene_id")

fdata$exprmatMkrs %>% 
  apply(1, function(vals) { (vals - mean(vals)) / sd(vals) }) %>% t() %>% 
  session$graphingUtils$heatmap(size = session$graphingUtils$SMALL, annotation_row = xdata$anno, 
                                show_rownames = FALSE, show_colnames = FALSE)


# plot out the pairwise co-expression distribution within and across cell types

xdata$coexmat <- fdata$exprmatMkrs %>% t() %>% cor()

xdata$coextbl <- xdata$coexmat %>% workspace$utils$coexMatToTbl() %>% 
  left_join(fdata$mkrsTbl %>% dplyr::select(gene_a = gene_id, cell_type_a = cell_type)) %>% 
  left_join(fdata$mkrsTbl %>% dplyr::select(gene_b = gene_id, cell_type_b = cell_type)) %>% 
  mutate(cell_type_in = cell_type_a == cell_type_b)

xdata$coextbl %>% 
  session$graphingUtils$ggplot(aes(x = cor_coef)) +
  geom_density(aes(color = cell_type_in)) + 
  facet_wrap(~cell_type_a)

# now, let's perform PCA to get the marker genes profiles
fdata$mgps <- workspace$utils$computeMgps(fdata$mkrs$mkrsFlat, fdata$exprmatMkrs)

# now, apply to the whole transcriptome, correct it out!
fdata$ccvModels <- workspace$utils$fitCCVModels(fdata$mgps$main, exprmat = fdata$exprdat$exprmat)

# save 
# mgp; ccv models; expression matrix + samples

# =================================  
# COMMIT ==========================
fdata$output$mgps <- fdata$mgps
fdata$output$ccvModels <- fdata$ccvModels
fdata$output$exprdat$exprmat <- fdata$ccvModels$exprmats$residual
fdata$output$exprdat$samples <- fdata$exprdat$samples
fdata$output %>% saveRDS(paste0(workspace$outputDir, "bk_exprmats_rosmap-mgpres.rds"))
# =================================  

# TODO - TO BE DELETED
fdata$ccvModel <- readRDS(paste0(workspace$outputDir, "bk_exprmats_rosmap-mgpres.rds"))$ccvModels

fdata$ccvModel$stats$lmStats %>% session$graphingUtils$ggplot(aes(x = rsqr)) + geom_histogram()
fdata$ccvModel$stats$lmStats$rsqr %>% mean(na.rm = TRUE)
fdata$ccvModel$stats$lmStats %>% session$graphingUtils$ggplot(aes(x = pvalue)) + geom_histogram()
fdata$ccvModel$stats$lmStats %>% filter(qvalue < 0.1)

pdata$files$uexprdats <- tibble(outfile = list.files(workspace$outputDir)) %>% 
  filter(grepl("uexprdats", outfile)) %>% 
  mutate(outfile_name = gsub(".rds", "", outfile)) %>% 
  mutate(strs = strsplit(outfile_name, "_")) %>% 
  mutate(type = strs %>% sapply(function(currStr) { currStr[1] })) %>% 
  mutate(dataset = strs %>% sapply(function(currStr) { currStr[2] })) %>% 
  dplyr::select(outfile, type, dataset)


# +++++++++++++++++++++++++++++++++++++++++++++++++++
# next, develop the code to fit the dilution model per gene in Velmeshev and ROSMAP
# 1. Velmeshev 

fdata <- list(); gc()

fdata$exprdat$sc <- read_rds(paste0(workspace$outputDir, "uexprdats_velmeshev.rds"))
fdata$exprdat$bk <- read_rds(paste0(workspace$outputDir, "bk_exprmats_velmeshev.rds"))

fdata$exprmats$sbj <- fdata$exprdat$sc$sbj$exprmats %>% lapply(workspace$utils$logTrans) # log transform

fdata$models <- workspace$utils$fitModels(fdata$exprdat$bk$exprmat, fdata$exprmats$sbj) 

# add another one to use the MPG corrected version

fdata$models$lmStats %>% session$graphingUtils$ggplot(aes(x = rsqr)) + geom_histogram()
fdata$models$lmStats$rsqr %>% mean()

fdata$models$lmStats %>% session$graphingUtils$ggplot(aes(pvalue)) + geom_histogram()
fdata$models$lmStats$pvalue %>% mean() 
fdata$models$lmStats %>% filter(qvalue < 0.1) # how many pass qvalue < 0.1

# let's take a look at the relationship between variance and predictability  
fdata$variability <- fdata$exprdat$sc$sbj$exprmats %>% sapply(function(exprmat) {
  exprmat %>% apply(1, function(vals) { sd(vals) })
})

fdata$variabilityFc <- pdata$celltypes %>% lapply(function(celltype) {
  
  targetExpr <- fdata$variability[, celltype]
  othrExprMax <- fdata$variability[, colnames(fdata$variability) != celltype] %>% apply(1, max) 
  
  tibble(gene_id = names(targetExpr)) %>% 
    mutate(expr = targetExpr[gene_id], 
           othr_expr_max = othrExprMax[gene_id]) %>% 
    mutate(min_fc = expr / othr_expr_max) %>% 
    mutate(cell_type = celltype)
  
}) %>% session$dataWrangler$rbind()

fdata$models$coefStats %>% dplyr::select(gene_id, cell_type, rsqr_indep) %>% 
  left_join(fdata$variabilityFc, by = c("gene_id", "cell_type")) %>% 
  session$graphingUtils$ggplot(aes(x = log2(min_fc + 1), y = rsqr_indep)) + 
  geom_point(shape = 1) + 
  facet_wrap(~cell_type)
 
# =================================  
# COMMIT ==========================
fdata$output$lmModel <- fdata$models 
fdata$output$exprmat$sbj <- fdata$exprmats$sbj
fdata$output$exprmat$bkOrig <- fdata$models$exprmats$orig
fdata$output$exprmat$bkRes <- fdata$models$exprmats$residual
fdata$output$xSbjSd <- list(sd = fdata$variability, minfc = fdata$variabilityFc)
fdata$output %>% saveRDS(paste0(workspace$outputDir, "bk_exprmats_velmeshev-psdbkres.rds"))
# =================================  


# +++++++++++++++++++++++++++++++++++++++++++++++++++
# next, develop the code to fit the dilution model per gene in Velmeshev and ROSMAP
# 2. ROSMAP 

fdata <- list(); gc()

fdata$exprdat$sc <- read_rds(paste0(workspace$outputDir, "uexprdats_rosmap.rds"))
fdata$exprdat$bk <- read_rds(paste0(workspace$outputDir, "bk_exprmats_rosmap.rds"))

fdata$exprmats$sbj <- fdata$exprdat$sc$sbj$exprmats %>% lapply(workspace$utils$logTrans) # log transform

# map the sbj names to bulk version
fdata$projids <- fdata$exprdat$bk$samples %>% dplyr::select(projid, sample) %>% filter(projid %in% colnames(fdata$exprmats$sbj$excitatory)) %>% 
  arrange(sample)
fdata$exprmats$sbj <- fdata$exprmats$sbj %>% lapply(function(exprmat) {
  copy <- exprmat[, fdata$projids$projid]
  colnames(copy) <- fdata$projids$sample
  return(copy)
})

fdata$models <- workspace$utils$fitModels(fdata$exprdat$bk$exprmat, fdata$exprmats$sbj) 

# add another one to use the MPG corrected version

fdata$models$lmStats %>% session$graphingUtils$ggplot(aes(x = rsqr)) + geom_histogram()
fdata$models$lmStats$rsqr %>% mean()

fdata$models$lmStats %>% session$graphingUtils$ggplot(aes(pvalue)) + geom_histogram()
fdata$models$lmStats$pvalue %>% mean() 
fdata$models$lmStats %>% filter(qvalue < 0.1) # how many pass qvalue < 0.1

# let's take a look at the relationship between variance and predictability  
fdata$variability <- fdata$exprdat$sc$sbj$exprmats %>% sapply(function(exprmat) {
  exprmat %>% apply(1, function(vals) { sd(vals) })
})

fdata$variabilityFc <- pdata$celltypes %>% lapply(function(celltype) {
  
  targetExpr <- fdata$variability[, celltype]
  othrExprMax <- fdata$variability[, colnames(fdata$variability) != celltype] %>% apply(1, max) 
  
  tibble(gene_id = names(targetExpr)) %>% 
    mutate(expr = targetExpr[gene_id], 
           othr_expr_max = othrExprMax[gene_id]) %>% 
    mutate(min_fc = expr / othr_expr_max) %>% 
    mutate(cell_type = celltype)
  
}) %>% session$dataWrangler$rbind()

fdata$models$coefStats %>% dplyr::select(gene_id, cell_type, rsqr_indep) %>% 
  left_join(fdata$variabilityFc, by = c("gene_id", "cell_type")) %>% 
  session$graphingUtils$ggplot(aes(x = log2(min_fc + 1), y = rsqr_indep)) + 
  geom_point(shape = 1) + 
  facet_wrap(~cell_type)

# =================================  
# COMMIT ==========================
fdata$output$lmModel <- fdata$models 
fdata$output$exprmat$sbj <- fdata$exprmats$sbj
fdata$output$exprmat$bkOrig <- fdata$models$exprmats$orig
fdata$output$exprmat$bkRes <- fdata$models$exprmats$residual
fdata$output$xSbjSd <- list(sd = fdata$variability, minfc = fdata$variabilityFc)
fdata$output %>% saveRDS(paste0(workspace$outputDir, "bk_exprmats_rosmap-psdbkres.rds"))
# =================================  





























