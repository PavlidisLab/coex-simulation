pdata <- list()

pdata$figuresDir <- paste0(workspace$workspaceDir, "042_ANALYSIS_05_ccv_FIGURES/")

pdata$sfilesDir <- paste0(workspace$workspaceDir, "042_S_FILES/")

# just reuse the simulations from the previous part

pdata$cteprf <- read_rds(paste0(workspace$outputDir, "seed_cteprf.rds"))

pdata$this <- read_rds(paste0(workspace$outputDir, "this_cell_sampling.rds"))

pdata$exprmats <- read_rds(paste0(workspace$outputDir, "exprmats_agg_cell_sampling.rds"))$exprmats

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++===========
# simulate co-expressed mkrs, assign new gene ids and attach to the rest of the expression matrix

fdata <- list(); gc()

fdata$this <- read_rds(paste0(workspace$outputDir, "this_cell_sampling.rds"))
fdata$mkrs <- pdata$cteprf$mkrs

fdata$celltypes <- fdata$mkrs$mkr_ct %>% unique() %>% sort()
names(fdata$celltypes) <- fdata$celltypes

set.seed(0); fdata$genes <- pdata$cteprf$cteprf %>% rownames() %>% setdiff(fdata$mkrs$gene_id) %>% sample(600)
fdata$genes <- tibble(gene_id = fdata$genes) %>% mutate(gene_id_coex = paste0("coex.", gene_id)) %>% mutate(is_mkr = FALSE)
fdata$genes <- fdata$genes %>% 
  rbind(fdata$mkrs %>% dplyr::select(gene_id) %>% mutate(gene_id_coex = paste0("coex.", gene_id)) %>% mutate(is_mkr = TRUE))

# construct the co-expression programs, in each cell type; concatenate 600 genes to be co-expressed, in every cell type with all the markers as well
fdata$coexPrgms <- list(thePrgm = list(celltypes = fdata$celltypes, genes = fdata$genes$gene_id, value = 0.3))

fdata$this <- fdata$this %>% workspace$api$setRefDatCoexPrgms(fdata$coexPrgms, level = "cel")
fdata$this <- fdata$this %>% workspace$api$setRefDatCoexPrgms(fdata$coexPrgms, level = "sbj")

# just strip out the genes that I don't want to simulate from the mdparams
fdata$this$mdParams$cel <- fdata$this$mdParams$cel %>% lapply(function(mdParam) { mdParam[fdata$genes$gene_id, ] })
fdata$this$mdParams$sbj <- fdata$this$mdParams$sbj %>% lapply(function(mdParam) { mdParam[fdata$genes$gene_id, ] })

# =================================  
# COMMIT ==========================
fdata$genes %>% saveRDS(paste0(workspace$outputDir, "coex_genes_ccv.rds"))
fdata$this %>% saveRDS(paste0(workspace$outputDir, "this_coexgenes_ccv.rds"))
# =================================  

# ++++++++++++++++++
# now simulate cells
# in this case you actually need to generate enough cells per cell type 
# go with a large superset to work with, maybe 5k cells per cell type then sample with replacement later
# just 15 samples? 

fdata <- list(); gc()

fdata$this <- read_rds(paste0(workspace$outputDir, "this_coexgenes_ccv.rds"))

fdata$cellsBsln <- fdata$this %>% workspace$api$simBseLnCels(nSubject = 30, nCell = 5000) 

fdata$sbjMeans <- fdata$this %>% workspace$api$simSbjLvMeans(30)

fdata$celParamsSbjv <- fdata$this %>% workspace$api$computeCelParams(fdata$sbjMeans$simSbjs) 

fdata$cellsSbjv <- fdata$cellsBsln$simCells %>% workspace$api$convertCelLvDist(fdata$cellsBsln$params$md, fdata$celParamsSbjv)

# =================================  
# COMMIT ==========================
fdata$simObjs <- list()
fdata$simObjs$cellsBsln <- list(simdat = fdata$cellsBsln$simCells, params = fdata$cellsBsln$params)
fdata$simObjs$cellsSbjv <- list(simdat = fdata$cellsSbjv, params = fdata$celParamsSbjv)
fdata$simObjs$subjects <- list(simdat = fdata$sbjMeans$simSbjs, params = fdata$sbjMeans$params)

fdata$simObjs %>% saveRDS(paste0(workspace$outputDir, "simObjs_coexgenes_ccv.rds"))
# ================================= 

# ++++++===============
# now sample from the 5k to form cell type specific xSbj expression profiles

fdata <- list(); gc()

fdata$simObj <- read_rds(paste0(workspace$outputDir, "simObjs_coexgenes_ccv.rds"))
fdata$samples <- fdata$simObj$cellsBsln$simdat$samples

fdata$celltypes <- c("excitatory", "inhibitory", "oligodendrocyte", "astrocyte", "microglia", "opc") %>% session$dataWrangler$attachNames()
fdata$nCels <- c(9750, 3250, 8775, 5850, 1950, 2925) %>% session$dataWrangler$attachNames(useNames = fdata$celltypes)

fdata$sbjs <- fdata$samples$subject %>% unique() %>% sort() %>% session$dataWrangler$attachNames()

set.seed(0); fdata$cells <- fdata$celltypes %>% session$collectionUtils$lapply(function(celltype) {
  fdata$sbjs %>% session$collectionUtils$lapply(function(currSbj) { 
    currSamples <- fdata$samples %>% filter(subject == currSbj, cell_type == celltype)
    bootstrap <- currSamples$sample %>% sample(fdata$nCels[celltype], replace = TRUE)
    currSamples <- tibble(sample = bootstrap) %>% left_join(currSamples, by = "sample")
    return(currSamples)
  }, verbose = FALSE) %>% session$dataWrangler$rbind()
})

# now compute the sample level expressions

fdata$exprmats$bsln <- fdata$cells %>% session$collectionUtils$lapply(function(currCells) {
  workspace$utils$aggregSbj(fdata$simObj$cellsBsln$simdat$exprmat, samples = currCells)
})

fdata$exprmats$sbjv <- fdata$cells %>% session$collectionUtils$lapply(function(currCells) {
  workspace$utils$aggregSbj(fdata$simObj$cellsSbjv$simdat$exprmat, samples = currCells)
})

fdata$exprmats$islv <- fdata$celltypes %>% session$collectionUtils$lapply(function(celltype) {
  samples <- fdata$simObj$subjects$simdat$samples %>% filter(cell_type == celltype)
  exprmat <- fdata$simObj$subjects$simdat$exprmat[, samples$sample]
  colnames(exprmat) <- samples$subject
  return(exprmat)
})

# COMMIT =====
fdata$output <- list()
fdata$output$cells <- fdata$cells 
fdata$output$exprmats <- fdata$exprmats
fdata$output %>% saveRDS(paste0(workspace$outputDir, "exprmats_agg_coexgenes_ccv.rds"))
# ============

# +=================++++++++++++++===
# now load it back in
# and attach them to the original, with modified gene names
# change name of the markers in the original as well. 

fdata <- list(); gc()

fdata$exprmats <- read_rds(paste0(workspace$outputDir, "exprmats_agg_coexgenes_ccv.rds"))$exprmats

fdata$exprmats <- fdata$exprmats %>% lapply(function(currExprmats) {
  currExprmats %>% lapply(function(exprmat) {
    rownames(exprmat) <- paste0("coex.", rownames(exprmat))
    return(exprmat)
  })
})

fdata$exprmatsConcat <- pdata$exprmats %>% session$collectionUtils$lapplyWithName(function(level, currExprmats) {
  currExprmats %>% session$collectionUtils$lapplyWithName(function(celltype, exprmat) {
    exprmatCoex <- fdata$exprmats[[level]][[celltype]]
    exprmat %>% rbind(exprmatCoex[, colnames(exprmat)])
  })
})

# SOFT COMMITS ============
pdata$exprmats <- fdata$exprmatsConcat
# ==============================================


# +++++++++++++====================
# come up with estimates on how to vary cell type proportions

fdata <- list()

fdata$ihc <- read.table(paste0(workspace$dataDir, "patrick_ihc/ihc_ct_astro.txt"), header = TRUE) %>% t() %>% as_tibble(rownames = "sample") %>% mutate(cell_type = "astrocyte") %>% 
  rbind(read.table(paste0(workspace$dataDir, "patrick_ihc/ihc_ct_endo.txt"), header = TRUE) %>% t() %>% as_tibble(rownames = "sample") %>% mutate(cell_type = "endothelial")) %>% 
  rbind(read.table(paste0(workspace$dataDir, "patrick_ihc/ihc_ct_microglia.txt"), header = TRUE) %>% t() %>% as_tibble(rownames = "sample") %>% mutate(cell_type = "microglia")) %>% 
  rbind(read.table(paste0(workspace$dataDir, "patrick_ihc/ihc_ct_neuro.txt"), header = TRUE) %>% t() %>% as_tibble(rownames = "sample") %>% mutate(cell_type = "neuron")) %>% 
  rbind(read.table(paste0(workspace$dataDir, "patrick_ihc/ihc_ct_oligo.txt"), header = TRUE) %>% t() %>% as_tibble(rownames = "sample") %>% mutate(cell_type = "oligodendrocyte")) %>% 
  dplyr::select(sample, ct_prop = V1, cell_type)

fdata$ihc %>% 
  session$graphingUtils$ggplot(aes(x = sample, y = ct_prop, fill = cell_type)) + geom_bar(stat = "identity") + 
  session$graphingUtils$tiltX(angle = 90)

fdata$ihc %>% group_by(sample) %>% mutate(x = ct_prop / sum(ct_prop)) %>% ungroup() %>% group_by(cell_type) %>% summarize(mean(x)) 

fdata$ihcSmry <- fdata$ihc %>% 
  group_by(cell_type) %>% summarize(ct_prop = mean(ct_prop)) %>% 
  mutate(ct_prop_norm = round(ct_prop / sum(ct_prop), 2)) %>% 
  arrange(desc(ct_prop_norm))

fdata$ihcSmry %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = ct_prop_norm)) + 
  geom_bar(stat = "identity") + 
  geom_text(aes(label = ct_prop_norm), hjust = -0.3, size = 6) + 
  xlim(rev(fdata$ihcSmry$cell_type)) + 
  ylim(0, 0.5) +
  coord_flip() + 
  xlab("Cell type proportion") + 
  ylab("")

fdata$rs <- fdata$ihc %>% group_by(cell_type) %>% summarize(m = mean(ct_prop), s = sd(ct_prop)) %>% mutate(rs = s / m)
# coefficient of variation 

fdata$rs$rs %>% max()
fdata$rs$rs %>% mean() 
fdata$rs$rs %>% min() 

# =============================== FINAL NUMBERS OVER HERE
# 0.1, 0.25,  0.5 ? something like that?
# ========================================================

# +++++++++++++++++++++++++++++++++++++++============
# generate CCV & compute corresponding expression data (with AND without)

fdata <- list()

fdata$celltypes <- c("excitatory", "inhibitory", "oligodendrocyte", "astrocyte", "microglia", "opc") %>% session$dataWrangler$attachNames()
fdata$nCels <- c(9750, 3250, 8775, 5850, 1950, 2925) %>% session$dataWrangler$attachNames(useNames = fdata$celltypes)
fdata$fracCells <- fdata$nCels / sum(fdata$nCels)

fdata$sdCells <- list(cv01 = 0.1, cv03 = 0.3, cv05 = 0.5)
fdata$sdCells <- fdata$sdCells %>% lapply(function(cv) {
  fdata$nCels * cv
})

set.seed(0); fdata$ccvs <- fdata$sdCells %>% lapply(function(currSds) {
  ccv <- fdata$celltypes %>% sapply(function(celltype) {
    rand <- rnorm(30, mean = fdata$nCels[celltype], sd = currSds[celltype]) * (fdata$fracCells[celltype] / fdata$nCels[celltype])
    rand <- abs(rand) # take the absolute fraction to prevent negative numbers; shouldn't have a big impact on SD
  })
  rownames(ccv) <- pdata$exprmats$sbjv$excitatory %>% colnames()
  return(ccv)
})


# now compute the diluted version of the ccvs
fdata$ccvsNull <- fdata$ccvs %>% lapply(function(ccv) {
  libsize <- rep(1, nrow(ccv))
  names(libsize) <- rownames(ccv)
  ccvNull <- fdata$fracCells %>% sapply(function(currFrac) {
    libsize * currFrac
  })
  return(ccvNull)
})

# finally can put together the exprmats

# first the versions with CCV
fdata$exprmats$ccv <- fdata$ccvs %>% lapply(function(currCCV) {
  Reduce("+", fdata$celltypes %>% lapply(function(celltype) {
    ccv <- currCCV[, celltype]
    exprmat <- pdata$exprmats$sbjv[[celltype]] 
    exprmat %>% apply(1, function(vals) { vals * ccv }) %>% t()
  }))
})

# then the nulls
fdata$exprmats$null <- fdata$ccvsNull %>% lapply(function(currCCV) {
  Reduce("+", fdata$celltypes %>% lapply(function(celltype) {
    ccv <- currCCV[, celltype]
    exprmat <- pdata$exprmats$sbjv[[celltype]] 
    exprmat %>% apply(1, function(vals) { vals * ccv }) %>% t()
  }))
})

fdata$exprmats$sbjv <- pdata$exprmats$sbjv # add in the cell type specific expression profiles for reference

fdata$output <- list()
fdata$output$ccvs <- fdata$ccvs
fdata$output$ccvsNull <- fdata$ccvsNull
fdata$output$exprmats <- fdata$exprmats

# COMMIT =====
fdata$output %>% saveRDS(paste0(workspace$outputDir, "exprmats_ccv.rds"))
# ============

# SOFT COMMIT ===
pdata$exprdat <- read_rds(paste0(workspace$outputDir, "exprmats_ccv.rds"))
# ======================

# ++++++++++++++++===============
# first, analysis
# put together the RSQRS (both diluted + CCVs)
# keep the markers in there - get rid of the co-expressed ones

fdata <- list()

fdata$ccvs <- pdata$exprdat$ccvs
fdata$exprmats <- pdata$exprdat$exprmats

fdata$ccvModels <- fdata$ccvs %>% session$collectionUtils$lapplyWithName(function(ccvLvl, ccv) {
  workspace$utils$fitCCVModels(ccv, fdata$exprmats$ccv[[ccvLvl]]) 
})

fdata$nullModels <- fdata$exprmats$ccv %>% session$collectionUtils$lapplyWithName(function(ccvLvl, exprmat) {
  workspace$utils$computeRSqr(exprmat, pdata$exprdat$exprmats$null[[ccvLvl]]) 
})


# let's put together the rsqrs

fdata$rsqrs <- fdata$ccvModels %>% session$collectionUtils$lapplyWithName(function(ccvLvl, ccvModels) {
  ccvModels$stats$lmStats %>% dplyr::select(gene_id, rsqr) %>% mutate(ccv_lvl = ccvLvl, type = "ccv")
}) %>% session$dataWrangler$rbind()

fdata$rsqrs <- fdata$rsqrs %>% rbind(fdata$nullModels %>% session$collectionUtils$lapplyWithName(function(ccvLvl, nullModels) {
  nullModels %>% dplyr::select(gene_id, rsqr) %>% mutate(ccv_lvl = ccvLvl, type = "null")
}) %>% session$dataWrangler$rbind())

fdata$output <- list()
fdata$output$ccvModels <- fdata$ccvModels
fdata$output$nullModels <- fdata$nullModels
fdata$output$rsqrs <- fdata$rsqrs

# COMMIT the models & the rsqrs +==================
fdata$output %>% saveRDS(paste0(workspace$outputDir, "ccvModels_ccv.rds"))
# ===============================

# +++++++++++++++============
# now plot the rsqr against rsqrs

fdata <- list(); gc()

fdata$models <- read_rds(paste0(workspace$outputDir, "ccvModels_ccv.rds"))

fdata$rsqrs <- fdata$models$rsqrs

# remove all the marker genes 
# remove all the coex genes

# coex actualy means coreg

fdata$rsqrs %>% 
  filter(!grepl("^coex", gene_id)) %>% 
  filter(type == "ccv") %>% 
  mutate(is_mkr = gene_id %in% pdata$cteprf$mkrs$gene_id) %>% 
  session$graphingUtils$ggplot(aes(x = is_mkr, y = rsqr)) + 
  geom_boxplot() + 
  facet_wrap(~ccv_lvl)

fdata$rsqrs %>% 
  filter(!grepl("^coex", gene_id)) %>% # get rid of co-regulated genes
  group_by(ccv_lvl, type) %>% summarize(rsqr = mean(rsqr), n_gene = n())


# =================================================
# OUTPUT S_FILE
xdata <- list()

xdata$fileName <- "sfile_07_ccv_effects.tsv" 

xdata$dat <- read_rds(paste0(workspace$outputDir, "ccvModels_ccv.rds"))$rsqrs
xdata$dat <- xdata$dat %>% spread(type, rsqr) %>% filter(!grepl("^coex", gene_id))

colnames(xdata$dat) <- c("gene_id", "ccv_lvl", "rsqr_ccv", "rsqr_baseline")

write.table(xdata$dat, file = paste0(pdata$sfilesDir, xdata$fileName), 
            sep = "\t", row.names = FALSE)
# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_01.eps" # vector file too big just draw bitmap

# let's plot cv01 to cv05 separately
xdata$main <- fdata$rsqrs %>% filter(!grepl("^coex", gene_id), ccv_lvl == "cv01")

xdata$plot <- xdata$main %>% 
  spread(type, rsqr) %>% 
  session$graphingUtils$ggplot(aes(x = null, y = ccv)) + 
  geom_point(shape = 1, color = RColorBrewer::brewer.pal(9, "Blues")[rev(c(8))]) +
  ylim(0, 1) + 
  xlim(0, 1) + 
  xlab(bquote("PVE by baseline bulk expression"~(R^2))) + 
  ylab(bquote("PVE by bulk CCV"~(R^2))) +
  ggtitle("CCV Coefficient of variation = 0.1") + 
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 14))

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 15, height = 15) 
# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_02.eps" # vector file too big just draw bitmap

# let's plot cv01 to cv05 separately
xdata$main <- fdata$rsqrs %>% filter(!grepl("^coex", gene_id), ccv_lvl == "cv03")

xdata$plot <- xdata$main %>% 
  spread(type, rsqr) %>% 
  session$graphingUtils$ggplot(aes(x = null, y = ccv)) + 
  geom_point(shape = 1, color = RColorBrewer::brewer.pal(9, "Blues")[rev(c(6))]) +
  ylim(0, 1) + 
  xlim(0, 1) + 
  xlab(bquote("PVE by baseline bulk expression"~(R^2))) + 
  ylab(bquote("PVE by bulk CCV"~(R^2))) +
  ggtitle("CCV Coefficient of variation = 0.3") + 
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 14))

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 15, height = 15) 
# =================================================

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_03.eps" # vector file too big just draw bitmap

# let's plot cv01 to cv05 separately
xdata$main <- fdata$rsqrs %>% filter(!grepl("^coex", gene_id), ccv_lvl == "cv05")

xdata$plot <- xdata$main %>% 
  spread(type, rsqr) %>% 
  session$graphingUtils$ggplot(aes(x = null, y = ccv)) + 
  geom_point(shape = 1, color = RColorBrewer::brewer.pal(9, "Blues")[rev(c(4))]) +
  ylim(0, 1) + 
  xlim(0, 1) + 
  xlab(bquote("PVE by baseline bulk expression"~(R^2))) + 
  ylab(bquote("PVE by bulk CCV"~(R^2))) +
  ggtitle("CCV Coefficient of variation = 0.5") + 
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 14))

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 15, height = 15) 
# =================================================


# =================================================
# OUTPUT S_FILE
xdata <- list()

xdata$fileName <- "sfile_08_ccv_simulation_genes.tsv"

xdata$dat <- read_rds(paste0(workspace$outputDir, "coex_genes_ccv.rds"))

write.table(xdata$dat, file = paste0(pdata$sfilesDir, xdata$fileName), 
            sep = "\t", row.names = FALSE) 
# =================================================

# +++++++++++++++++++++++================
# show that this induced co-expressions 
# 2 gene sets: 
# 1. 600 co-expressed genes
# 2. 600 non-coexpressed genes (matching gene names)

fdata <- list(); gc()

fdata$genes <- read_rds(paste0(workspace$outputDir, "coex_genes_ccv.rds"))
fdata$genes <- fdata$genes %>% filter(!is_mkr)

fdata$exprmats$uncoex <- pdata$exprdat$exprmats$ccv %>% 
  lapply(function(exprmat) { exprmat[fdata$genes$gene_id, ] }) %>% 
  c(list(null = pdata$exprdat$exprmats$null$cv01[fdata$genes$gene_id, ]))

fdata$exprmats$coex <- pdata$exprdat$exprmats$ccv %>% 
  lapply(function(exprmat) { exprmat[fdata$genes$gene_id_coex, ] }) %>% 
  c(list(null = pdata$exprdat$exprmats$null$cv01[fdata$genes$gene_id_coex, ]))

fdata$coexmats$uncoex <- fdata$exprmats$uncoex %>% session$collectionUtils$lapply(workspace$utils$computeCoexmat)
fdata$coexmats$coex <- fdata$exprmats$coex %>% session$collectionUtils$lapply(workspace$utils$computeCoexmat)

fdata$coextbls$uncoex <- fdata$coexmats$uncoex %>% session$collectionUtils$lapplyWithName(function(ccvlvl, coexmat) {
  tibble(pair_id = workspace$utils$getPairIds(coexmat), coex = workspace$utils$vectorize(coexmat)) %>% 
    mutate(ccv_lvl = ccvlvl)
}) %>% session$dataWrangler$rbind()

fdata$coextbls$coex <- fdata$coexmats$coex %>% session$collectionUtils$lapplyWithName(function(ccvlvl, coexmat) {
  tibble(pair_id = workspace$utils$getPairIds(coexmat), coex = workspace$utils$vectorize(coexmat)) %>% 
    mutate(ccv_lvl = ccvlvl)
}) %>% session$dataWrangler$rbind()

# =========== SOFT COMMIT 
# so I can put the co-expressions together to write into a single supplementary file
pdata$coextbls$uncorrected <- fdata$coextbls
# ========================

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_04.eps" # vector file too big just draw bitmap

# let's plot cv01 to cv05 separately
xdata$main <- fdata$coextbls$uncoex %>% filter(ccv_lvl != "null") 

xdata$main %>% group_by(ccv_lvl) %>% summarize(coex = mean(coex))

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_density(data = fdata$coextbls$uncoex %>% filter(ccv_lvl == "null"), size = 1) +
  geom_density(aes(color = ccv_lvl), size = 1) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  theme(legend.position = "none") +
  xlim(-1, 1) + 
  scale_color_manual(values = RColorBrewer::brewer.pal(9, "Blues")[rev(c(4, 6, 8))]) + 
  xlab("Co-expression (Pearson's r)") + 
  ylab("Density") +
  ggtitle("Distribution of pairwise gene co-expression for\nsimulated genes not undergoing co-regulation") + 
  theme(plot.title = element_text(size = 17), plot.subtitle = element_text(size = 14))

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 15, height = 15) 
# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_05.eps" # vector file too big just draw bitmap

# let's plot cv01 to cv05 separately
xdata$main <- fdata$coextbls$coex %>% filter(ccv_lvl != "null") 

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_density(data = fdata$coextbls$coex %>% filter(ccv_lvl == "null"), size = 1) +
  geom_density(aes(color = ccv_lvl), size = 1) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  theme(legend.position = "none") +
  xlim(-1, 1) + 
  scale_color_manual(values = RColorBrewer::brewer.pal(9, "Blues")[rev(c(4, 6, 8))]) + 
  xlab("Co-expression (Pearson's r)") + 
  ylab("Density") +
  ggtitle("Distribution of pairwise gene co-expression for\nsimulated genes undergoing active co-regulation") + 
  theme(plot.title = element_text(size = 17), plot.subtitle = element_text(size = 14))

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 15, height = 15) 
# =================================================


# +++++++++++++++++++++++================
# show that I can correct out the co-expressions using the ground truth CCV vectors

fdata <- list(); gc()

fdata$genes <- read_rds(paste0(workspace$outputDir, "coex_genes_ccv.rds"))
fdata$genes <- fdata$genes %>% filter(!is_mkr)

fdata$models <- read_rds(paste0(workspace$outputDir, "ccvModels_ccv.rds"))

# get the residuals from the CCV models
fdata$exprmatsCrrted <- fdata$models$ccvModels %>% lapply(function(ccvModel) {
  ccvModel$exprmats$residual
})

fdata$exprmats$uncoex <- fdata$exprmatsCrrted %>% 
  lapply(function(exprmat) { exprmat[fdata$genes$gene_id, ] }) %>% 
  c(list(null = pdata$exprdat$exprmats$null$cv01[fdata$genes$gene_id, ]))

fdata$exprmats$coex <- fdata$exprmatsCrrted %>% 
  lapply(function(exprmat) { exprmat[fdata$genes$gene_id_coex, ] }) %>% 
  c(list(null = pdata$exprdat$exprmats$null$cv01[fdata$genes$gene_id_coex, ]))

fdata$coexmats$uncoex <- fdata$exprmats$uncoex %>% session$collectionUtils$lapply(workspace$utils$computeCoexmat)

fdata$coexmats$coex <- fdata$exprmats$coex %>% session$collectionUtils$lapply(workspace$utils$computeCoexmat)

fdata$coextbls$uncoex <- fdata$coexmats$uncoex %>% session$collectionUtils$lapplyWithName(function(ccvlvl, coexmat) {
  tibble(pair_id = workspace$utils$getPairIds(coexmat), coex = workspace$utils$vectorize(coexmat)) %>% 
    mutate(ccv_lvl = ccvlvl)
}) %>% session$dataWrangler$rbind()

fdata$coextbls$coex <- fdata$coexmats$coex %>% session$collectionUtils$lapplyWithName(function(ccvlvl, coexmat) {
  tibble(pair_id = workspace$utils$getPairIds(coexmat), coex = workspace$utils$vectorize(coexmat)) %>% 
    mutate(ccv_lvl = ccvlvl)
}) %>% session$dataWrangler$rbind()

fdata$coextbls$uncoex %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_density(aes(color = ccv_lvl), size = 1) + 
  theme(legend.position = "none") +
  xlim(-1, 1)

fdata$coextbls$coex %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_density(aes(color = ccv_lvl), size = 1) + 
  theme(legend.position = "none") +
  xlim(-1, 1)

# =========== SOFT COMMIT 
# so I can put the co-expressions together to write into a single supplementary file
pdata$coextbls$grdtruth <- fdata$coextbls
# ========================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_06.eps" # vector file too big just draw bitmap

# let's plot cv01 to cv05 separately
xdata$main <- fdata$coextbls$uncoex %>% filter(ccv_lvl != "null") 

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_density(data = fdata$coextbls$uncoex %>% filter(ccv_lvl == "null"), linewidth = 1) +
  geom_density(aes(color = ccv_lvl), size = 1) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  theme(legend.position = "none") +
  xlim(-1, 1) + 
  scale_color_manual(values = RColorBrewer::brewer.pal(9, "Blues")[rev(c(4, 6, 8))]) + 
  xlab("Co-expression (Pearson's r)") + 
  ylab("Density") +
  ggtitle("Correction by ground truth cell type proportions\nfor genes not undergoing co-regulation") + 
  theme(plot.title = element_text(size = 17), plot.subtitle = element_text(size = 14))

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 15, height = 15) 
# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_07.eps" # vector file too big just draw bitmap

# let's plot cv01 to cv05 separately
xdata$main <- fdata$coextbls$coex %>% filter(ccv_lvl != "null") 

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_density(data = fdata$coextbls$coex %>% filter(ccv_lvl == "null"), linewidth = 1) +
  geom_density(aes(color = ccv_lvl), size = 1) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  theme(legend.position = "none") +
  xlim(-1, 1) + 
  scale_color_manual(values = RColorBrewer::brewer.pal(9, "Blues")[rev(c(4, 6, 8))]) + 
  xlab("Co-expression (Pearson's r)") + 
  ylab("Density") +
  ggtitle("Correction by ground truth cell type proportions\nfor genes undergoing active co-regulation") + 
  theme(plot.title = element_text(size = 17), plot.subtitle = element_text(size = 14))

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 15, height = 15) 
# =================================================

# +++++++++++++++++++++++================
# show that MGPs can capture CCV under certain circumstances

fdata <- list(); gc(); 

fdata$celltypes <- pdata$cteprf$mkrs$mkr_ct %>% unique() %>% sort() %>% session$dataWrangler$attachNames()

fdata$mkrs$indep <- fdata$celltypes %>% lapply(function(celltype) {
  pdata$cteprf$mkrs %>% filter(mkr_ct == celltype) %>% session$dataWrangler$extractColumn("gene_id")
})

fdata$mkrs$coex <- fdata$celltypes %>% lapply(function(celltype) {
  mkrs <- pdata$cteprf$mkrs %>% filter(mkr_ct == celltype) %>% session$dataWrangler$extractColumn("gene_id")
  paste0("coex.", mkrs)
})

fdata$exprmats <- pdata$exprdat$exprmats$ccv

fdata$mgps <- fdata$mkrs %>% session$collectionUtils$lapply(function(mkrs) {
  fdata$exprmats %>% session$collectionUtils$lapply(function(exprmat) {
    currExprmat <- exprmat %>% session$dataWrangler$setRownameAsColumn("gene_id")
    mgp <- markerGeneProfile::mgpEstimate(currExprmat, mkrs, geneColName = "gene_id", geneTransform = NULL)
    mgp$main <- mgp$estimates %>% sapply(function(vals) { vals }) # marker gene profiles (one vector per cell type)
    return(mgp)
  })
})

fdata$mgpsMain <- fdata$mgps %>% session$collectionUtils$lapply(function(currMgps) {
  currMgps %>% session$collectionUtils$lapply(function(mgps) {
    mgps$main
  })
})

fdata$output <- list()
fdata$output$mgps <- fdata$mgps
fdata$output$main <- fdata$mgpsMain

# +====++++++++++===========
# COMMIT
fdata$output %>% saveRDS(paste0(workspace$outputDir, "mgps_ccv.rds"))
# =================================

# +===========================
# now let's check how close the MGPs are to the ground truth CCVs

fdata <- list()

fdata$ccvs <- pdata$exprdat$ccvs
fdata$mgps <- read_rds(paste0(workspace$outputDir, "mgps_ccv.rds"))$main

fdata$celltypes <- pdata$cteprf$mkrs$mkr_ct %>% unique() %>% sort() %>% session$dataWrangler$attachNames()

# first how well do the independent markers correlate with the ccvs?
fdata$mgpccvCorcoef <- fdata$mgps$indep %>% session$collectionUtils$lapplyWithName(function(ccvLvl, mgps) {
  ccvs <- fdata$ccvs[[ccvLvl]]
  fdata$celltypes %>% lapply(function(celltype) {
    mgp <- mgps[, celltype]
    ccv <- ccvs[names(mgp), celltype]
    tibble(cell_type = celltype, ccv_lvl = ccvLvl, cor_coef = cor(mgp, ccv), mkrs_coex = FALSE)
  }) %>% session$dataWrangler$rbind()
}) %>% session$dataWrangler$rbind()

# now what if the mkrs are coexpressed? 
fdata$mgpccvCorcoef <- fdata$mgpccvCorcoef %>% rbind(fdata$mgps$coex %>% session$collectionUtils$lapplyWithName(function(ccvLvl, mgps) {
  ccvs <- fdata$ccvs[[ccvLvl]]
  fdata$celltypes %>% lapply(function(celltype) {
    mgp <- mgps[, celltype]
    ccv <- ccvs[names(mgp), celltype]
    tibble(cell_type = celltype, ccv_lvl = ccvLvl, cor_coef = cor(mgp, ccv), mkrs_coex = TRUE)
  }) %>% session$dataWrangler$rbind()
}) %>% session$dataWrangler$rbind())

fdata$mgpccvCorcoef %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = cor_coef)) + 
  geom_point(aes(color = ccv_lvl)) + 
  geom_line(aes(color = ccv_lvl, group = ccv_lvl)) + 
  facet_wrap(~mkrs_coex) + 
  session$graphingUtils$tiltX(angle = 90) + 
  ylim(0, 1)




# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_12.eps" 

xdata$plot <- fdata$mgpccvCorcoef %>% 
  session$graphingUtils$ggplot(aes(x = ccv_lvl, y = cor_coef)) + 
  geom_boxplot(aes(color = mkrs_coex), position = position_dodge(width = 0)) + 
  geom_point(aes(color = mkrs_coex)) + 
  ylim(0, 1) + 
  scale_color_brewer(palette = "Dark2") + 
  ggtitle("MGPs capture ground truth cell type\nproportions") + 
  theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 14)) + 
  xlab("CCV coefficient of variation") + 
  ylab("Correlation with ground truth\ncell type proportions (Pearson's r)") + 
  theme(legend.position = "none")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 15, height = 15) 
# =================================================


# +===========================
# now let's use the MGPs to correct co-expression and see how well it does

fdata <- list()

fdata$mgps <- read_rds(paste0(workspace$outputDir, "mgps_ccv.rds"))$main

fdata$exprmats <- pdata$exprdat$exprmats

# first the independent mkrs
fdata$ccvModels$indep <- fdata$mgps$indep %>% session$collectionUtils$lapplyWithName(function(ccvLvl, ccv) {
  workspace$utils$fitCCVModels(ccv, fdata$exprmats$ccv[[ccvLvl]]) 
})

# next the coex mkrs
fdata$ccvModels$coex <- fdata$mgps$coex %>% session$collectionUtils$lapplyWithName(function(ccvLvl, ccv) {
  workspace$utils$fitCCVModels(ccv, fdata$exprmats$ccv[[ccvLvl]]) 
})

# let's put together the rsqrs

fdata$rsqrs <- fdata$ccvModels$indep %>% session$collectionUtils$lapplyWithName(function(ccvLvl, ccvModels) {
  ccvModels$stats$lmStats %>% dplyr::select(gene_id, rsqr) %>% mutate(ccv_lvl = ccvLvl, type = "indep")
}) %>% session$dataWrangler$rbind()

fdata$rsqrs <- fdata$rsqrs %>% rbind(fdata$ccvModels$coex %>% session$collectionUtils$lapplyWithName(function(ccvLvl, ccvModels) {
  ccvModels$stats$lmStats %>% dplyr::select(gene_id, rsqr) %>% mutate(ccv_lvl = ccvLvl, type = "coex")
}) %>% session$dataWrangler$rbind())

# fdata$rsqrs %>% session$graphingUtils$ggplot(aes(x = rsqr)) + geom_density(aes(color = ccv_lvl)) + facet_wrap(~type)

fdata$output <- list()
fdata$output$ccvModels <- fdata$ccvModels
fdata$output$rsqrs <- fdata$rsqrs

# COMMIT ========================
fdata$output %>% saveRDS(paste0(workspace$outputDir, "mgpModels_ccv.rds"))
# ==================================

# ++++++++++++++++++++++====================
# now what happens to the co-expressions after corrections?
# let's see what happens to the genes that ARE NOT CO-EXPRESSED

fdata <- list(); gc()

fdata$genes <- read_rds(paste0(workspace$outputDir, "coex_genes_ccv.rds"))
fdata$genes <- fdata$genes %>% filter(!is_mkr)

fdata$models <- read_rds(paste0(workspace$outputDir, "mgpModels_ccv.rds"))

# get the residuals from the CCV models
fdata$exprmatsCrrted <- fdata$models$ccvModels %>% lapply(function(ccvModels) {
  ccvModels %>% lapply(function(model) { model$exprmats$residual })
})

fdata$exprmats$indep <- fdata$exprmatsCrrted$indep %>% 
  lapply(function(exprmat) { exprmat[fdata$genes$gene_id, ] }) %>% 
  c(list(null = pdata$exprdat$exprmats$null$cv01[fdata$genes$gene_id, ]))

fdata$exprmats$depen <- fdata$exprmatsCrrted$coex %>% 
  lapply(function(exprmat) { exprmat[fdata$genes$gene_id, ] }) %>% 
  c(list(null = pdata$exprdat$exprmats$null$cv01[fdata$genes$gene_id, ]))

fdata$coexmats$indep <- fdata$exprmats$indep %>% session$collectionUtils$lapply(workspace$utils$computeCoexmat)

fdata$coexmats$depen <- fdata$exprmats$depen %>% session$collectionUtils$lapply(workspace$utils$computeCoexmat)

fdata$coextbls$indep <- fdata$coexmats$indep %>% session$collectionUtils$lapplyWithName(function(ccvlvl, coexmat) {
  tibble(pair_id = workspace$utils$getPairIds(coexmat), coex = workspace$utils$vectorize(coexmat)) %>% 
    mutate(ccv_lvl = ccvlvl)
}) %>% session$dataWrangler$rbind()

fdata$coextbls$depen <- fdata$coexmats$depen %>% session$collectionUtils$lapplyWithName(function(ccvlvl, coexmat) {
  tibble(pair_id = workspace$utils$getPairIds(coexmat), coex = workspace$utils$vectorize(coexmat)) %>% 
    mutate(ccv_lvl = ccvlvl)
}) %>% session$dataWrangler$rbind()

fdata$coextbls$indep %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_density(aes(color = ccv_lvl), size = 1) + 
  theme(legend.position = "none") +
  xlim(-1, 1)

fdata$coextbls$depen %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_density(aes(color = ccv_lvl), size = 1) + 
  xlim(-1, 1)

# =========== SOFT COMMIT 
# so I can put the co-expressions together to write into a single supplementary file
pdata$coextbls$mgpUncoex <- fdata$coextbls
# ========================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_08.eps" # vector file too big just draw bitmap

# let's plot cv01 to cv05 separately
xdata$main <- fdata$coextbls$indep %>% filter(ccv_lvl != "null") 

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_density(data = fdata$coextbls$indep %>% filter(ccv_lvl == "null"), linewidth = 1) +
  geom_density(aes(color = ccv_lvl), size = 1) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  theme(legend.position = "none") +
  xlim(-1, 1) + 
  scale_color_manual(values = RColorBrewer::brewer.pal(9, "Blues")[rev(c(4, 6, 8))]) + 
  xlab("Co-expression (Pearson's r)") + 
  ylab("Density") +
  ggtitle("Correction by independent markers derived MGPs\nfor genes not undergoing co-regulation") + 
  theme(plot.title = element_text(size = 17), plot.subtitle = element_text(size = 14))

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 15, height = 15) 
# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_09.eps" # vector file too big just draw bitmap

# let's plot cv01 to cv05 separately
xdata$main <- fdata$coextbls$depen %>% filter(ccv_lvl != "null") 

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_density(data = fdata$coextbls$depen %>% filter(ccv_lvl == "null"), linewidth = 1) +
  geom_density(aes(color = ccv_lvl), size = 1) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  theme(legend.position = "none") +
  xlim(-1, 1) + 
  scale_color_manual(values = RColorBrewer::brewer.pal(9, "Blues")[rev(c(4, 6, 8))]) + 
  xlab("Co-expression (Pearson's r)") + 
  ylab("Density") +
  ggtitle("Correction by co-regulated markers derived MGPs\nfor genes undergoing active co-regulation") + 
  theme(plot.title = element_text(size = 17), plot.subtitle = element_text(size = 14))

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 15, height = 15) 
# =================================================


# ++++++++++++++++++++++====================
# now what happens to the co-expressions after corrections?
# let's see what happens to the genes that ARE IN FACT CO-EXPRESSED

fdata <- list(); gc()

fdata$genes <- read_rds(paste0(workspace$outputDir, "coex_genes_ccv.rds"))
fdata$genes <- fdata$genes %>% filter(!is_mkr)

fdata$models <- read_rds(paste0(workspace$outputDir, "mgpModels_ccv.rds"))

# get the residuals from the CCV models
fdata$exprmatsCrrted <- fdata$models$ccvModels %>% lapply(function(ccvModels) {
  ccvModels %>% lapply(function(model) { model$exprmats$residual })
})

fdata$exprmats$indep <- fdata$exprmatsCrrted$indep %>% 
  lapply(function(exprmat) { exprmat[fdata$genes$gene_id_coex, ] }) %>% 
  c(list(null = pdata$exprdat$exprmats$null$cv01[fdata$genes$gene_id_coex, ]))

fdata$exprmats$depen <- fdata$exprmatsCrrted$coex %>% 
  lapply(function(exprmat) { exprmat[fdata$genes$gene_id_coex, ] }) %>% 
  c(list(null = pdata$exprdat$exprmats$null$cv01[fdata$genes$gene_id_coex, ]))

fdata$coexmats$indep <- fdata$exprmats$indep %>% session$collectionUtils$lapply(workspace$utils$computeCoexmat)

fdata$coexmats$depen <- fdata$exprmats$depen %>% session$collectionUtils$lapply(workspace$utils$computeCoexmat)

fdata$coextbls$indep <- fdata$coexmats$indep %>% session$collectionUtils$lapplyWithName(function(ccvlvl, coexmat) {
  tibble(pair_id = workspace$utils$getPairIds(coexmat), coex = workspace$utils$vectorize(coexmat)) %>% 
    mutate(ccv_lvl = ccvlvl)
}) %>% session$dataWrangler$rbind()

fdata$coextbls$depen <- fdata$coexmats$depen %>% session$collectionUtils$lapplyWithName(function(ccvlvl, coexmat) {
  tibble(pair_id = workspace$utils$getPairIds(coexmat), coex = workspace$utils$vectorize(coexmat)) %>% 
    mutate(ccv_lvl = ccvlvl)
}) %>% session$dataWrangler$rbind()

fdata$coextbls$indep %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_density(aes(color = ccv_lvl), size = 1) + 
  theme(legend.position = "none") +
  xlim(-1, 1)

fdata$coextbls$depen %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_density(aes(color = ccv_lvl), size = 1) + 
  xlim(-1, 1)

# =========== SOFT COMMIT 
# so I can put the co-expressions together to write into a single supplementary file
pdata$coextbls$mgpCoex <- fdata$coextbls
# ========================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_10.eps" # vector file too big just draw bitmap

# let's plot cv01 to cv05 separately
xdata$main <- fdata$coextbls$indep %>% filter(ccv_lvl != "null") 

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_density(data = fdata$coextbls$indep %>% filter(ccv_lvl == "null"), linewidth = 1) +
  geom_density(aes(color = ccv_lvl), size = 1) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  theme(legend.position = "none") +
  xlim(-1, 1) + 
  scale_color_manual(values = RColorBrewer::brewer.pal(9, "Blues")[rev(c(4, 6, 8))]) + 
  xlab("Co-expression (Pearson's r)") + 
  ylab("Density") +
  ggtitle("Correction by independent markers derived MGPs\nfor genes undergoing active co-regulation") + 
  theme(plot.title = element_text(size = 17), plot.subtitle = element_text(size = 14))

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 15, height = 15) 
# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_11.eps" # vector file too big just draw bitmap

# let's plot cv01 to cv05 separately
xdata$main <- fdata$coextbls$depen %>% filter(ccv_lvl != "null") 

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_density(data = fdata$coextbls$depen %>% filter(ccv_lvl == "null"), linewidth = 1) +
  geom_density(aes(color = ccv_lvl), size = 1) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  theme(legend.position = "none") +
  xlim(-1, 1) + 
  scale_color_manual(values = RColorBrewer::brewer.pal(9, "Blues")[rev(c(4, 6, 8))]) + 
  xlab("Co-expression (Pearson's r)") + 
  ylab("Density") +
  ggtitle("Correction by co-regulated markers derived MGPs\nfor genes not undergoing co-regulation") + 
  theme(plot.title = element_text(size = 17), plot.subtitle = element_text(size = 14))

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 15, height = 15) 
# =================================================

# =================================================
# OUTPUT S_FILE
xdata <- list()

xdata$fileName <- "sfile_06_cell_type_proportions.tsv" 

xdata$dat <- read_rds(paste0(workspace$outputDir, "exprmats_ccv.rds"))

xdata$dat$ccvs$cv00 <- xdata$dat$ccvsNull$cv01
xdata$dat$ccvs <- xdata$dat$ccvs[sort(names(xdata$dat$ccvs))]

xdata$main <- xdata$dat$ccvs %>% session$collectionUtils$lapplyWithName(function(ccvLvl, ccvMat) { 
  ccvMat %>% session$dataWrangler$setRownameAsColumn("sbj_id") %>% mutate(ccv_lvl = ccvLvl)
}) %>% session$dataWrangler$rbind()

write.table(xdata$main, file = paste0(pdata$sfilesDir, xdata$fileName), 
            sep = "\t", row.names = FALSE)
# =================================================


# =================================================
# OUTPUT S_FILE
xdata <- list()

xdata$fileName <- "sfile_09_ccv_coex_values.tsv" 

xdata$dat <- pdata$coextbls$uncorrected %>% lapply(function(tbl) { tbl %>% mutate(correction = "none") }) %>% 
  c(pdata$coextbls$grdtruth %>% lapply(function(tbl) { tbl %>% mutate(correction = "groundtruth") })) %>% 
  c(list(pdata$coextbls$mgpUncoex$indep %>% mutate(correction = "mgp_indep_mkrs"))) %>% 
  c(list(pdata$coextbls$mgpUncoex$depen %>% mutate(correction = "mgp_coreg_mkrs"))) %>% 
  c(list(pdata$coextbls$mgpCoex$indep %>% mutate(correction = "mgp_indep_mkrs"))) %>% 
  c(list(pdata$coextbls$mgpCoex$depen %>% mutate(correction = "mgp_coreg_mkrs"))) %>% 
  session$dataWrangler$rbind()

write.table(xdata$dat, file = paste0(pdata$sfilesDir, xdata$fileName), 
            sep = "\t", row.names = FALSE)
# =================================================





