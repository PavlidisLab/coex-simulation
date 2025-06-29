pdata <- list()

# +++++++++++++++++++++
# load in allen institute data
# 1. clean (filter for genes detected in at least 10% of all cells, filter samples 2 SD away), per cell type
# 2. normalize the full table, missins fill in with 0 
# 3. find marker genes, which are * 3-fold more highly expressed than all others
# 4. sample 1k non-marker genes to include + all the markers; correlations only imposed for the non-markers

fdata <- list()

fdata$samples <- (1 %>% lapply(function(i) {
  samples <- read_csv(paste0(workspace$dataDir, "allen_sn_all/", "metadata.csv")) %>% 
    dplyr::select(sample = sample_name, 
                  subject = external_donor_name_label, 
                  region = region_label, 
                  cluster = cluster_label)
  mapping <- read_csv(paste0(workspace$dataDir, "allen_sn_all/", "sample_mapping.csv"))
  samples <- samples %>% left_join(mapping %>% dplyr::select(sample = sample_name, sample_id = exp_component_name), by = "sample")
  samples <- samples %>% dplyr::select(sample = sample_id, subject, region, cluster)
  return(samples)
}))[[1]]

fdata$samples <- fdata$samples %>% mutate(cell_type = cluster %>% sapply(function(currCluster) { 
  if (grepl("^Inh", currCluster)) { "inhibitory" } 
  else if (grepl("^Exc", currCluster)) { "excitatory" } 
  else if (grepl("^Astro", currCluster)) { "astrocyte" } 
  else if (grepl("^Micro", currCluster)) { "microglia" } 
  else if (grepl("^Oligo", currCluster)) { "oligodendrocyte" } 
  else if (grepl("^OPC", currCluster)) { "opc" } 
  else { NA }
})) %>% na.omit()

fdata$ctmat <- workspace$utils$readCounts(paste0(workspace$dataDir, "allen_sn_all/", "exon.csv")) 
fdata$ctmat <- fdata$ctmat[, fdata$samples$sample]

fdata$celltypes <- fdata$samples$cell_type %>% unique()
names(fdata$celltypes) <- fdata$celltypes

# per cell type, filter for 1. genes that are detected in at least 10% of all cells, and 2. cells that are beyond +-2SD * N genes detected
fdata$vDims <- fdata$ctmat %>% workspace$utils$cleanCtmat()
fdata$ctmat <- fdata$ctmat[fdata$vDims$vGenes, fdata$vDims$vSamples]
fdata$samples <- fdata$samples %>% filter(sample %in% colnames(fdata$ctmat))

# normalize
fdata$exprmat <- fdata$ctmat %>% workspace$utils$normalizeCtmat()
fdata$logexprmat <- log2(fdata$exprmat + 1)

fdata$allenSn <- list()
fdata$allenSn$samples <- fdata$samples
fdata$allenSn$ctmat <- fdata$ctmat
fdata$allenSn$exprmat <- fdata$exprmat
fdata$allenSn$logexprmat <- fdata$logexprmat

# =================================  
# COMMIT ==========================
fdata$allenSn %>% saveRDS(paste0(workspace$outputDir, "ref_data_allen_sn.rds"))
# =================================

# +++++++++++++++++++++
# Now use the loaded data to look at cell type expression profiles
# 1. identify the marker genes
# 2. sample 1k random genes to include
# 3. prefix mkr (marker) and nam (not a marker) for each gene, concatenate them >> ctprf

fdata <- list(); gc()

fdata$allenSn <- read_rds(paste0(workspace$outputDir, "ref_data_allen_sn.rds"))

fdata$samples <- fdata$allenSn$samples

fdata$logexprmat <- fdata$allenSn$logexprmat

fdata$celltypes <- fdata$samples$cell_type %>% unique() %>% sort()
names(fdata$celltypes) <- fdata$celltypes

fdata$cteprfLog <- fdata$celltypes %>% sapply(function(currCelltype) {
  currSamples <- fdata$samples %>% filter(cell_type == currCelltype)
  currExprmat <- fdata$logexprmat[, currSamples$sample]
  result <- currExprmat %>% apply(1, mean)
  return(result)
})

# let's filter for genes that are detected in all cell types - limit the simulation to only those expressed at some baseline level
fdata$minmax <- fdata$cteprfLog %>% apply(1, function(vals) { c(min = min(vals), max = max(vals)) }) %>% t() %>% 
  session$dataWrangler$setRownameAsColumn("gene")
fdata$minmax <- fdata$minmax %>% filter(min > 0, max <= 10)
fdata$cteprfLog <- fdata$cteprfLog[intersect(rownames(fdata$cteprfLog), fdata$minmax$gene), ] # filter out extreme values

# insepct overall distribution of these genes; I think they need to be qunatile normalized
fdata$cteprfLog %>% 
  session$dataWrangler$setRownameAsColumn("gene_id") %>% 
  gather(cell_type, expr, -gene_id) %>% 
  session$graphingUtils$ggplot(aes(x = expr)) + 
  geom_density(aes(color = cell_type))

# quantile normalize
fdata$cteprfLog <- fdata$cteprfLog %>% workspace$utils$quantileNormalize()

fdata$mkrs <- fdata$celltypes %>% sapply(function(currCelltype) {
  currCt <- fdata$cteprfLog[, currCelltype]
  othrCtMax <- fdata$cteprfLog[, colnames(fdata$cteprfLog) != currCelltype] %>% apply(1, max)
  (currCt / othrCtMax)
}) %>% 
  as_tibble(rownames = "gene_id") %>% 
  gather(cell_type, minfc, -gene_id)

fdata$mkrs <- fdata$mkrs %>% group_by(cell_type) %>% mutate(minfc_rank = rank(-minfc)) %>% ungroup()

fdata$mkrs <- fdata$mkrs %>% 
  dplyr::select(gene_id, mkr_ct = cell_type, minfc, minfc_rank) %>% 
  filter(minfc_rank <= 20) # pick 20 mkrs per cell_type

fdata$mkrs %>% group_by(mkr_ct) %>% summarize(n = n())

# now sample 1k regular genes and then all the marker genes
# we will simulate expression for these genes, preserving their cell type expression profiles

set.seed(0); fdata$geneIds <- rownames(fdata$cteprfLog) %>% 
  setdiff(fdata$mkrs$gene_id) %>% 
  sample((3000 - nrow(fdata$mkrs))) %>% 
  c(fdata$mkrs$gene_id)

# get the exprmat (cpm) cteprf because this is the scale on which we are fitting the models

fdata$cteprf <- fdata$celltypes %>% sapply(function(currCelltype) {
  currSamples <- fdata$samples %>% filter(cell_type == currCelltype)
  currExprmat <- fdata$allenSn$exprmat[fdata$geneIds, currSamples$sample] # filter for only the genes we want
  result <- currExprmat %>% apply(1, mean)
  return(result)
})

# insepct overall distribution of these genes; I think they need to be qunatile normalized
fdata$cteprf %>% 
  session$dataWrangler$setRownameAsColumn("gene_id") %>% 
  gather(cell_type, expr, -gene_id) %>% 
  session$graphingUtils$ggplot(aes(x = log2(expr + 1))) + 
  geom_density(aes(color = cell_type))

# quantile normalize
fdata$cteprf <- fdata$cteprf %>% workspace$utils$quantileNormalize()

# ================================= 
# COMMIT ==========================
fdata$output <- list()
fdata$output$cteprf <- fdata$cteprf
fdata$output$mkrs <- fdata$mkrs # only include 60 markers, 10 top markers per cell type
fdata$output %>% saveRDS(paste0(workspace$outputDir, "seed_cteprf.rds"))
# =================================

# +++++++++++++++++++++
# Now load in the LCL data for subject level variability
# 1. filter out genes not detected
# 2. filter out extreme samples
# 3. normalize / log transform

fdata <- list()

fdata$samples <- read.table(paste0(workspace$dataDir, "lappalainen/", "E-GEUV-3.sdrf.txt"), 
                            sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% as_tibble()

fdata$samples <- fdata$samples %>% dplyr::select(sample =  Assay.Name, 
                                                 subject = Source.Name, 
                                                 performer = Performer) %>% unique()

fdata$subjects <- read.table(paste0(workspace$dataDir, "lappalainen/", "E-GEUV-1.sdrf.txt"), 
                             sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% as_tibble()

fdata$subjects <- fdata$subjects %>% dplyr::select(subject = Characteristics.individual., sex = Characteristics.sex.) %>% unique()

fdata$samples <- fdata$samples %>% left_join(fdata$subjects, by = "subject") %>% na.omit()

fdata$ctmat <- read.table(paste0(workspace$dataDir, "lappalainen/", "GD660.GeneQuantCount.txt"), 
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE) 

fdata$samples <- fdata$samples %>% filter(sample %in% colnames(fdata$ctmat)) %>% arrange(sample)

fdata$geneIds <- fdata$ctmat$TargetID
fdata$ctmat <- fdata$ctmat[, fdata$samples$sample] %>% as.matrix()
rownames(fdata$ctmat) <- fdata$geneIds

# per cell type, filter for 1. genes that are detected in at least 10% of all cells, and 2. cells that are beyond +-2SD * N genes detected
fdata$vDims <- fdata$ctmat %>% workspace$utils$cleanCtmat()
fdata$ctmat <- fdata$ctmat[fdata$vDims$vGenes, fdata$vDims$vSamples]
fdata$samples <- fdata$samples %>% filter(sample %in% colnames(fdata$ctmat))

fdata$samples <- fdata$samples %>% 
  group_by(subject) %>% 
  mutate(n_sample = n()) %>% # there are 5 subjects with 7-8 samples! - might be good to look at technical variability
  arrange(desc(n_sample)) %>% 
  ungroup() %>% 
  filter(n_sample == 1) %>% 
  dplyr::select(-n_sample)

fdata$ctmat <- fdata$ctmat[, fdata$samples$sample]

# normalize
fdata$exprmat <- fdata$ctmat %>% workspace$utils$normalizeCtmat() %>% workspace$utils$quantileNormalize()
fdata$logexprmat <- log2(fdata$exprmat + 1)

# yeah, let's quantile normalize
fdata$exprmat[, 1:30] %>% 
  session$dataWrangler$setRownameAsColumn("gene_id") %>% 
  gather(sample, expr, -gene_id) %>% 
  session$graphingUtils$ggplot(aes(x = sample, y = log2(expr + 1))) + 
  geom_boxplot()

# ================================= let's try this; now every time you commit to pdata, write it to snapshot so this point can be retrieved easily. 
# COMMIT ==========================
fdata$lclBk <- list()
fdata$lclBk$samples <- fdata$samples
fdata$lclBk$ctmat <- fdata$ctmat
fdata$lclBk$exprmat <- fdata$exprmat
fdata$lclBk$logexprmat <- fdata$logexprmat

fdata$lclBk %>% saveRDS(paste0(workspace$outputDir, "ref_data_lcl_bk.rds"))
# =================================

# +++++++++++++++++++++
# Set up all the input data for the simulation here
# 1. cell level reference data
# 2. subject level reference data
# 3. cell type expression profile (annotated as marker vs. non-marker)
# 4. cell type specific correlations

fdata <- list()

fdata$allenSn <- read_rds(paste0(workspace$outputDir, "ref_data_allen_sn.rds"))
fdata$lclBk <- read_rds(paste0(workspace$outputDir, "ref_data_lcl_bk.rds"))
fdata$cteprf <- read_rds(paste0(workspace$outputDir, "seed_cteprf.rds"))

fdata$simulator <- workspace$initSimulator()
fdata$this <- fdata$simulator$this

# cell level reference data

# ** prepare the cell level md reference data

xdata <- list()

# H200.1030 MTG excitatory 4674 has the highest number of cells - also from surgical specimen so probably very good quality
# H200.1023 MTG excitatory  4480 has the second most number of cells
fdata$allenSn$samples %>% group_by(subject, region, cell_type) %>% summarize(n_cel = n()) %>% ungroup() %>% arrange(desc(n_cel))

xdata$samples <- fdata$allenSn$samples %>% filter(subject == "H200.1030", region == "MTG", cell_type == "excitatory") 
xdata$exprmat <- fdata$allenSn$exprmat[, xdata$samples$sample]

fdata$this <- fdata$this %>% fdata$simulator$api$setRefDatMd(xdata$exprmat, level = "cel")

# ** prepare the subject level md reference data

xdata <- list()

xdata$samples <- fdata$lclBk$samples

xdata$samples <- xdata$samples %>% dplyr::select(sample, subject) %>% unique()

# set.seed(0); xdata$samples <- xdata$samples %>% group_by(subject) %>% summarize(sample = sample(sample, 1)) # retain only a single sample per subject

xdata$exprmat <- xdata$exprmat <- fdata$lclBk$exprmat[, xdata$samples$sample]

fdata$this <- fdata$this %>% fdata$simulator$api$setRefDatMd(xdata$exprmat, level = "sbj")

# now prepare the cell type expression profiles

fdata$this <- fdata$this %>% fdata$simulator$api$setRefDatCteprf(fdata$cteprf$cteprf)

# now prepare the cormats

xdata <- list()

xdata$cteprf <- fdata$cteprf$cteprf

xdata$genes <- xdata$cteprf %>% rownames()

# generate length(xdata$genes) of size 3, using rnorm(), in order to generate a random cormat consisting of all correlation values
# do this per cell type

xdata$celltypes <- colnames(xdata$cteprf)
names(xdata$celltypes) <- xdata$celltypes

set.seed(0); xdata$genesCoex_1 <- xdata$genes %>% sample(300, replace = FALSE)
set.seed(2); xdata$genesCoex_2 <- xdata$genes %>% sample(150, replace = FALSE)

xdata$coexPrgms <- list()

xdata$coexPrgms$prgm_1 <- list(
  celltypes = c("astrocyte", "excitatory"), 
  genes = xdata$genesCoex_1, 
  value = 0.9
)

xdata$coexPrgms$prgm_2 <- list(
  celltypes = c("microglia"), 
  genes = xdata$genesCoex_2, 
  value = 0.8
)

fdata$this <- fdata$this %>% fdata$simulator$api$setRefDatCoexPrgms(xdata$coexPrgms, level = "sbj") # use exactly the same cormats for both sbj and cel levels
fdata$this <- fdata$this %>% fdata$simulator$api$setRefDatCoexPrgms(xdata$coexPrgms, level = "cel") # use exactly the same cormats for both sbj and cel levels

# make sure you have all the ingredients in place

fdata$this$refdat$exprmat$cel %>% dim()
fdata$this$refdat$exprmat$sbj %>% dim()

fdata$this$refdat$cteprf %>% dim()

fdata$this$refdat$coexPrograms$cel %>% lapply(function(program) { program$celltypes })
fdata$this$refdat$coexPrograms$cel %>% lapply(function(program) { length(program$genes) })
fdata$this$refdat$coexPrograms$cel %>% lapply(function(program) { program$value })

fdata$this$refdat$coexPrograms$sbj %>% lapply(function(program) { program$celltypes })
fdata$this$refdat$coexPrograms$sbj %>% lapply(function(program) { length(program$genes) })
fdata$this$refdat$coexPrograms$sbj %>% lapply(function(program) { program$value })

# fit the models to get the marginal distributions for the specified genes
fdata$this <- fdata$this %>% fdata$simulator$api$fitMds()


# =================================  
# COMMIT ==========================
fdata$this %>% saveRDS(paste0(workspace$outputDir, "sim_test.rds"))
# =================================  

# ++++++++++++++++++
# now simulate cells

fdata <- list()

fdata$this <- read_rds(paste0(workspace$outputDir, "sim_test.rds"))

fdata$cellsBsln <- fdata$this %>% workspace$api$simBseLnCels(nSubject = 100, nCell = 100) 

fdata$sbjMeans <- fdata$this %>% workspace$api$simSbjLvMeans(100) 

fdata$celParamsSbjv <- fdata$this %>% workspace$api$computeCelParams(fdata$sbjMeans$simSbjs) 

fdata$cellsSbjv <- fdata$cellsBsln$simCells %>% workspace$api$convertCelLvDist(fdata$cellsBsln$params$md, fdata$celParamsSbjv)

# =================================  
# COMMIT ==========================
fdata$simObjs <- list()
fdata$simObjs$cellsBsln <- list(simdat = fdata$cellsBsln$simCells, params = fdata$cellsBsln$params)
fdata$simObjs$cellsSbjv <- list(simdat = fdata$cellsSbjv, params = fdata$celParamsSbjv)
fdata$simObjs$subjects <- list(simdat = fdata$sbjMeans$simSbjs, params = fdata$sbjMeans$params)

fdata$simObjs %>% saveRDS(paste0(workspace$outputDir, "sim_test_simObjs.rds"))
# =================================  

# SOFT COMMIT ====================
fdata <- list(); gc()
fdata$simObjs <- read_rds(paste0(workspace$outputDir, "sim_test_simObjs.rds"))
pdata$simObjs <- fdata$simObjs
# ================================

# ++++++++++++++++++++++=============================
# sniff test of all the mean vs variance levels here

# start with mean & variance of the subject levels
fdata <- list()
fdata$celltype <- "excitatory"
fdata$samples <- pdata$simObjs$subjects$simdat$samples %>% filter(cell_type == fdata$celltype)
fdata$exprmat <- pdata$simObjs$subjects$simdat$exprmat[, fdata$samples$sample]

fdata$mm <- fdata$exprmat %>% 
  apply(1, function(vals) { c(m = mean(vals), v = var(vals)) }) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble()

fdata$mm %>% session$graphingUtils$ggplot(aes(x = log2(m + 1), y = log2(v + 1))) + geom_point()
fdata$mm %>% arrange(desc(m))

# baseline cell level mean & variance

fdata <- list()
fdata$celltype <- "excitatory"
fdata$subject <- "sbj_001"
fdata$samples <- pdata$simObjs$cellsBsln$simdat$samples %>% filter(cell_type == fdata$celltype, subject == fdata$subject)
fdata$exprmat <- pdata$simObjs$cellsBsln$simdat$exprmat[, fdata$samples$sample]

fdata$mm <- fdata$exprmat %>% 
  apply(1, function(vals) { c(m = mean(vals), v = var(vals)) }) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble()

fdata$mm %>% session$graphingUtils$ggplot(aes(x = log2(m + 1), y = log2(v + 1))) + geom_point()
fdata$mm %>% arrange(desc(m))

# sbjv cell level mean & variance

fdata <- list()
fdata$celltype <- "excitatory"
fdata$subject <- "sbj_002"
fdata$samples <- pdata$simObjs$cellsBsln$simdat$samples %>% filter(cell_type == fdata$celltype, subject == fdata$subject)
fdata$exprmat <- pdata$simObjs$cellsSbjv$simdat$exprmat[, fdata$samples$sample]

fdata$mm <- fdata$exprmat %>% 
  apply(1, function(vals) { c(m = mean(vals), s = sd(vals)) }) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble()

fdata$mm %>% session$graphingUtils$ggplot(aes(x = log2(m + 1), y = log2(s + 1))) + geom_point()


# check the co-expressions of the intended genes

fdata <- list()

fdata$coexProgram <- pdata$simObjs$cellsBsln$params$coexPrograms$prgm_2
fdata$coexProgram$celltypes

fdata$celltype <- "microglia"
fdata$subject <- "sbj_013"
fdata$samples <- pdata$simObjs$cellsBsln$simdat$samples %>% filter(cell_type == fdata$celltype, subject == fdata$subject)
fdata$exprmat <- pdata$simObjs$cellsBsln$simdat$exprmat[, fdata$samples$sample]
fdata$coexmat <- fdata$exprmat %>% workspace$utils$computeCoexmat()
fdata$coextbl <- fdata$coexmat %>% 
  workspace$utils$getCoords() %>% 
  mutate(coex = fdata$coexmat %>% workspace$utils$vectorize())

fdata$coextbl <- fdata$coextbl %>% 
  mutate(gene_a_coex = (gene_a %in% fdata$coexProgram$genes), 
         gene_b_coex = (gene_b %in% fdata$coexProgram$genes)) %>% 
  mutate(n_genes_coex = as.numeric(gene_a_coex) + as.numeric(gene_b_coex))

fdata$coextbl %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_density(aes(color = as.character(n_genes_coex)))

fdata$coextbl %>% arrange(desc(coex)) %>% filter(n_genes_coex == 1) %>% arrange(desc(coex))

# fdata$exprmat %>% workspace$utils$plotScatter("PHF20L1", "NADSYN1")
  

# check the synchrony of the intended genes

fdata <- list()

fdata$coexProgram <- pdata$simObjs$cellsBsln$params$coexPrograms$prgm_1
fdata$coexProgram$celltypes

fdata$celltype <- "astrocyte"
fdata$samples <- pdata$simObjs$subjects$simdat$samples %>% filter(cell_type == fdata$celltype)
fdata$exprmat <- pdata$simObjs$subjects$simdat$exprmat[, fdata$samples$sample]
fdata$coexmat <- fdata$exprmat %>% workspace$utils$computeCoexmat()
fdata$coextbl <- fdata$coexmat %>% 
  workspace$utils$getCoords() %>% 
  mutate(coex = fdata$coexmat %>% workspace$utils$vectorize())

fdata$coextbl <- fdata$coextbl %>% 
  mutate(gene_a_coex = (gene_a %in% fdata$coexProgram$genes), 
         gene_b_coex = (gene_b %in% fdata$coexProgram$genes)) %>% 
  mutate(n_genes_coex = as.numeric(gene_a_coex) + as.numeric(gene_b_coex))

fdata$coextbl %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_density(aes(color = as.character(n_genes_coex)))

fdata$coextbl %>% arrange(desc(coex)) %>% filter(n_genes_coex == 2) %>% arrange(desc(coex))

fdata$exprmat %>% workspace$utils$plotScatter("HERPUD2", "ARFRP1")

# let's check the synchrony among cell types

fdata <- list()

fdata$coexProgram <- pdata$simObjs$cellsBsln$params$coexPrograms$prgm_1
fdata$coexProgram$celltypes
fdata$coexProgram$genes

fdata$gene <- "TMEM43"

fdata$celltypes <- pdata$simObjs$subjects$simdat$samples$cell_type %>% unique() %>% sort() %>% 
  session$dataWrangler$attachNames()

fdata$exprGene <- fdata$celltypes %>% sapply(function(celltype) {
  samples <- pdata$simObjs$subjects$simdat$samples %>% filter(cell_type == celltype)
  exprmat <- pdata$simObjs$subjects$simdat$exprmat[, samples$sample]
  exprmat[fdata$gene, ]
})

fdata$exprGene %>% cor() %>% session$graphingUtils$heatmap()

fdata$exprGene %>% 
  session$dataWrangler$setRownameAsColumn("subject") %>% 
  session$graphingUtils$ggplot(aes(x = astrocyte, y = excitatory)) + 
  geom_point()

# check the marginal distribution per gene, compared against the real data?

fdata <- list(); gc()

fdata$this <- read_rds(paste0(workspace$outputDir, "sim_test.rds"))

fdata$celltype <- "excitatory"
fdata$subject <- "sbj_013"
fdata$samples <- pdata$simObjs$cellsBsln$simdat$samples %>% filter(cell_type == fdata$celltype, subject == fdata$subject)
fdata$exprmat <- pdata$simObjs$cellsBsln$simdat$exprmat[, fdata$samples$sample]

fdata$mm <- fdata$exprmat %>% 
  apply(1, function(vals) { c(m = mean(vals), s = sd(vals)) }) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble()

fdata$gene <- "RTN1"

fdata$exprmat[fdata$gene, ] %>% session$dataWrangler$vectorToTibble() %>% 
  session$graphingUtils$ggplot(aes(x = value)) + geom_density()

fdata$this$refdat$exprmat$cel[fdata$gene, ] %>% session$dataWrangler$vectorToTibble() %>% 
  session$graphingUtils$ggplot(aes(x = value)) + geom_density()



