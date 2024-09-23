
pdata <- list()

pdata$figuresDir <- paste0(workspace$workspaceDir, "041_ANALYSIS_01_xcell_coex_FIGURES/")

pdata$genes <- readRDS(paste0(workspace$outputDir, "genes_metadata.rds")) # genes

pdata$exprmatFiles <- tibble(outfile = list.files(workspace$outputDir)) %>% 
  filter(grepl("exprmats", outfile)) %>% 
  mutate(outfile_name = gsub(".rds", "", outfile)) %>% 
  filter(!grepl("RAW", outfile_name)) %>% 
  mutate(strs = strsplit(outfile_name, "_")) %>% 
  mutate(level = strs %>% sapply(function(currStr) { currStr[1] })) %>% 
  mutate(dataset = strs %>% sapply(function(currStr) { currStr[3] })) %>% 
  dplyr::select(outfile, level, dataset)

pdata$celltypes <- c("excitatory", "inhibitory", "opc", "oligodendrocyte", "astrocyte", "microglia")
names(pdata$celltypes) <- pdata$celltypes

pdata$ctMkrs <- read_rds(paste0(workspace$outputDir, "mkrs_mancarci.rds"))

pdata$go <- readRDS(paste0(workspace$outputDir, "gene_ontology.rds"))

pdata$files$coexmats <- tibble(outfile = list.files(workspace$outputDir)) %>% 
  filter(grepl("coexmats", outfile)) %>% 
  mutate(outfile_name = gsub(".rds", "", outfile)) %>% 
  mutate(strs = strsplit(outfile_name, "_")) %>% 
  mutate(level = strs %>% sapply(function(currStr) { currStr[1] })) %>% 
  mutate(type = strs %>% sapply(function(currStr) { currStr[2] })) %>% 
  mutate(dataset = strs %>% sapply(function(currStr) { currStr[3] })) %>% 
  mutate(cell_type = strs %>% sapply(function(currStr) { currStr[4] })) %>% 
  dplyr::select(outfile, level, type, dataset, cell_type)


# ++++++++++++++++++++++++++++++=
# let's take a look at the number of subjects vs. samples

fdata <- list()

fdata$files <- pdata$exprmatFiles %>% filter(level == "sc")

fdata$filesFlat <- fdata$files$outfile
names(fdata$filesFlat) <- fdata$files$dataset

(fdata$currFile <- fdata$filesFlat[7])
fdata$currDat <- read_rds(paste0(workspace$outputDir, fdata$currFile))

fdata$currDat$raw$samples$sample_bk %>% unique() %>% sort() %>% length()
fdata$currDat$raw$samples$projid %>% unique() %>% sort() %>% length()


xdata <- list()
xdata$libsizes <- fdata$currDat$raw$ctmat %>% colSums()
xdata$libsizes %>% summary()
xdata$libsizes %>% length()

fdata$currDat$raw$ctmat %>% dim()

# ++++++++++++++++++++++++++++++=
# first analysis, compute the number of samples and cells per dataset, cell type combinations

fdata <- list()

fdata$files <- pdata$exprmatFiles %>% filter(level == "sc")

fdata$filesFlat <- fdata$files$outfile
names(fdata$filesFlat) <- fdata$files$dataset

fdata$nCells <- fdata$filesFlat %>% session$collectionUtils$lapplyWithName(function(currDataset, currFile) {
  currDat <- read_rds(paste0(workspace$outputDir, currFile))
  currDat$samples %>% session$dataWrangler$rbind() %>% 
    group_by(cell_type, sample_bk) %>% 
    summarize(n_cell = n()) %>% 
    ungroup() %>% 
    mutate(dataset = currDataset) %>% 
    dplyr::select(dataset, cell_type, sample_bk, n_cell)
}) %>% session$dataWrangler$rbind()

fdata$nCells <- fdata$nCells %>% filter(cell_type %in% pdata$celltypes)

fdata$nCells <- fdata$nCells %>% filter(n_cell >= 20)

fdata$nCells <- fdata$nCells %>% group_by(dataset, cell_type) %>% mutate(n_sbj = n()) %>% ungroup()

fdata$nCells <- fdata$nCells %>% filter(n_sbj >= 20) # remove the data group without enough sbjs

# =================================  
# COMMIT ==========================
fdata$nCells %>% saveRDS(paste0(workspace$outputDir, "sc_stats_ncells.rds"))
# =================================


# ++++++++++++++++++++++++++++++=
# figures time

fdata <- list()

fdata$quantile_breaks <- function(xs, n = 10) {
  seq(min(xs), max(xs), length.out = 10)
}

fdata$nCells <- read_rds(paste0(workspace$outputDir, "sc_stats_ncells.rds"))

fdata$nCells$n_cell %>% sum()

fdata$nCellsSmry <- fdata$nCells %>% group_by(dataset, cell_type) %>% summarize(n_cell = sum(n_cell)) %>% ungroup()

xdata <- list()

xdata$main <- fdata$nCellsSmry %>% spread(cell_type, n_cell) %>% session$dataWrangler$setColAsRownames("dataset")
xdata$main <- xdata$main %>% session$dataWrangler$fillNa("opc", 0) %>% as.matrix()

xdata$main %>% 
  session$graphingUtils$heatmap(decimalNums = 1, cluster_rows = FALSE, cluster_cols = FALSE)

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_01.eps"

xdata$plot <- fdata$nCellsSmry %>% 
  mutate(dataset = workspace$utils$fmtDataset(dataset), 
         cell_type = workspace$utils$fmtCelltypes(cell_type)) %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = n_cell, color = dataset)) + 
  geom_point(size = 4) + 
  geom_line(aes(group = dataset)) + 
  session$graphingUtils$tiltX(angle = 90) + 
  scale_y_continuous(trans = "log10", labels = label_comma()) + 
  ggtitle("Number of cells per group") + 
  xlab("Cell type") + 
  ylab("Number of cells") + 
  theme(legend.key.size = unit(1.5,"line")) + 
  labs(color = "Dataset")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 20, height = 15, dpi = 1000) 
# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "s_figure_01.eps"

xdata$main <- fdata$nCellsSmry %>% mutate(dataset = workspace$utils$fmtDataset(dataset), 
                                          cell_type = workspace$utils$fmtCelltypes(cell_type))

xdata$celltype <- xdata$main$cell_type %>% unique() %>% sort()

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = n_cell)) + 
  facet_wrap(~dataset, ncol = 1) + 
  geom_bar(stat = "identity", aes(fill = dataset)) + 
  xlab("Cell type") + 
  ylab("Number of cells") + 
  xlim(rev(xdata$celltype)) +
  coord_flip() +
  ggtitle("Breakdown of cells by data groups") + 
  scale_fill_brewer(palette = "Set1") + 
  theme(legend.position = "none")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 35, height = 40, dpi = 1000) 
# =================================================

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "s_figure_02.eps"

xdata$main <- fdata$nCells %>% group_by(dataset, cell_type) %>% summarize(n_subject = n()) %>% ungroup() # fdata$nSbjSmry

xdata$main <- xdata$main %>% mutate(dataset = workspace$utils$fmtDataset(dataset), 
                                          cell_type = workspace$utils$fmtCelltypes(cell_type))

xdata$celltype <- xdata$main$cell_type %>% unique() %>% sort()

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = n_subject)) + 
  facet_wrap(~dataset, ncol = 1) + 
  geom_bar(stat = "identity", aes(fill = dataset)) + 
  xlab("Cell type") + 
  ylab("Number of subjects") + 
  coord_flip() + 
  xlim(rev(xdata$celltype)) + 
  ggtitle("Breakdown of subjects by data groups") + 
  scale_fill_brewer(palette = "Set1") + 
  theme(legend.position = "none")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 35, height = 40, dpi = 1000) 
# =================================================
  

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_02.eps"

xdata$nCells <- fdata$nCells %>% 
  mutate(dataset = workspace$utils$fmtDataset(dataset), 
         cell_type = workspace$utils$fmtCelltypes(cell_type))

xdata$means <- xdata$nCells %>% group_by(dataset, cell_type) %>% 
  summarize(n_cell = mean(n_cell)) %>% 
  ungroup()

xdata$plot <- xdata$nCells %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = n_cell, color = dataset)) + 
  geom_point(position = "jitter", shape = 1, size = 1) + 
  geom_point(data = xdata$means, size = 4) +
  geom_line(data = xdata$means, aes(group = dataset), linewidth = 1) + 
  scale_y_continuous(trans = "log10", labels = label_comma()) + 
  ylab("Number of cells per subject") + 
  xlab("Cell type") + 
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 14), 
        legend.position = "none") + 
  session$graphingUtils$tiltX(angle = 90) + 
  scale_color_brewer(palette = "Set1") + 
  ggtitle("Number of cells per subject by data group")
  
xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 20, height = 15) 

# =================================================

fdata$nCells %>% group_by(dataset, cell_type) %>% summarize(n_cell = mean(n_cell))

# ++++++++++++++
# figures
# now number of genes per data group

fdata <- list()

fdata$main <- read_rds(paste0(workspace$outputDir, "sc_stats_ctprofiles.rds"))

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_03.eps"

xdata$main <- fdata$main %>% 
  group_by(dataset, cell_type) %>% summarize(n_gene = n()) %>% ungroup() 

xdata$main <- xdata$main %>% 
  mutate(dataset = workspace$utils$fmtDataset(dataset), 
         cell_type = workspace$utils$fmtCelltypes(cell_type))

xdata$plot <- xdata$main  %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = n_gene, color = dataset)) + 
  geom_point(size = 4) +
  geom_line(aes(group = dataset)) + 
  session$graphingUtils$tiltX(angle = 90) + 
  scale_y_continuous(labels = label_comma()) + 
  xlab("Cell type") + 
  ylab("Number of genes") + 
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 14), 
        legend.position = "none") + 
  session$graphingUtils$tiltX(angle = 90) + 
  scale_color_brewer(palette = "Set1") + 
  ggtitle("Number of genes detected by data group")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 20, height = 15) 

# =================================================


# this doesn't show anything beyond the main figures

fdata$main %>% 
  group_by(dataset, cell_type) %>% summarize(n_gene = n()) %>% ungroup() %>% 
  session$graphingUtils$ggplot(aes(x = n_gene, y = cell_type)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~dataset, ncol = 1) + 
  ylim(rev(pdata$celltypes)) + 
  ggtitle("Number of genes detected by dataset and cell type")

fdata$main %>% 
  group_by(dataset, cell_type) %>% summarize(n_gene = n()) %>% ungroup() %>% 
  group_by(cell_type) %>% 
  summarize(n_gene = mean(n_gene))

# +++ within dataset, compare cell types

fdata$datasets <- fdata$main$dataset %>% unique() %>% sort()

xdata <- list()

(xdata$dataset <- fdata$datasets[7])

xdata$main <- fdata$main %>% filter(dataset == xdata$dataset)

xdata$main <- xdata$main %>%
  dplyr::select(cell_type, gene_id) %>% 
  mutate(value = 1) %>% 
  spread(cell_type, value) %>% 
  session$dataWrangler$fillNa(colNames = fdata$celltypes, value = 0) %>% 
  session$dataWrangler$setColAsRownames("gene_id") 

xdata$main %>% UpSetR::upset(nsets = 7, text.scale = 2, keep.order = TRUE, nintersects = NA)


# +++ within cell type, compare datasets

xdata <- list()

(xdata$celltype <- pdata$celltypes[5])

xdata$main <- fdata$main %>% filter(cell_type == xdata$celltype)

xdata$main <- xdata$main %>%
  dplyr::select(dataset, gene_id) %>% 
  mutate(value = 1) %>% 
  spread(dataset, value) %>% 
  session$dataWrangler$fillNa(colNames = fdata$datasets, value = 0) %>% 
  session$dataWrangler$setColAsRownames("gene_id") 

xdata$main %>% UpSetR::upset(nsets = 7, text.scale = 2, keep.order = TRUE, nintersects = NA)

# +++++++++++++++++++++++++++++++++++=
# for each cell type, compute intersection vs. union; stack plot

fdata <- list(); gc()

fdata$main <- read_rds(paste0(workspace$outputDir, "sc_stats_ctprofiles.rds"))

fdata$celltypes <- fdata$main$cell_type %>% unique() %>% sort() %>% 
  session$dataWrangler$attachNames()

fdata$geneTally <- fdata$celltypes %>% session$collectionUtils$lapply(function(celltype) {
  
  geneTally <- fdata$main %>% 
    filter(cell_type == celltype) %>% 
    dplyr::select(dataset, gene_id) %>% 
    mutate(value = 1) %>% 
    spread(dataset, value) %>% 
    session$dataWrangler$fillNa(colNames = fdata$datasets, value = 0) %>% 
    session$dataWrangler$setColAsRownames("gene_id") 
  
  geneTally %>% apply(1, function(vals) { sum(vals) }) %>% 
    session$dataWrangler$vectorToTibble() %>% 
    dplyr::select(gene = variable, n_dataset = value) %>% 
    mutate(cell_type = celltype)
  
}) %>% session$dataWrangler$rbind()

fdata$geneTally <- fdata$geneTally %>% 
  group_by(cell_type) %>% mutate(frac_dataset = n_dataset / max(n_dataset)) %>% ungroup()

fdata$geneTally <- fdata$geneTally %>% 
  mutate(class = n_dataset %>% sapply(function(n) { min(n, 2) }) ) %>% 
  mutate(class = class + (frac_dataset %>% sapply(function(n) { if (n == 1) { n } else { 0 } }))) %>% 
  mutate(class = class %>% sapply(function(currClass) { 
    if (currClass >= 3) { "All" }
    else if (currClass >= 2) { "Multiple" }
    else { "Single" }
  }))

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_04_1.eps"

xdata$main <- fdata$geneTally %>% group_by(cell_type, class) %>% summarize(n_gene = n()) %>% ungroup()

xdata$main <- xdata$main %>% mutate(cell_type = workspace$utils$fmtCelltypes(cell_type))

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = n_gene)) + 
  geom_bar(aes(fill = factor(class, levels = c("Single", "Multiple", "All"))), stat = "identity") + 
  session$graphingUtils$tiltX(angle = 90) + 
  scale_fill_manual(values = rev(c("grey10", "grey40", "grey80"))) + 
  xlab("Cell type") + 
  ylab("Number of genes") + 
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 14)) + 
  labs(fill = "N Datasets") + 
  ggtitle("Number of genes intersecting datasets by cell type")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 20, height = 15) 

# =================================================


# ++++++++++++=
# now let's look at the cell type specific sets of genes commonly detected across all datasets

fdata <- list()

fdata$ctprofiles <- read_rds(paste0(workspace$outputDir, "sc_stats_ctprofiles.rds"))

fdata$datasets <- fdata$ctprofiles$dataset %>% unique() %>% sort()

fdata$tally <- pdata$celltypes %>% session$collectionUtils$lapply(function(currCelltype) {
  ctprofiles <- fdata$ctprofiles %>% filter(cell_type == currCelltype)
  
  mat <- ctprofiles %>%
    dplyr::select(dataset, gene_id) %>% 
    mutate(value = 1) %>% 
    spread(dataset, value) %>% 
    session$dataWrangler$fillNa(colNames = fdata$datasets, value = 0) %>% 
    session$dataWrangler$setColAsRownames("gene_id") 
  
  repro <- mat %>% apply(1, sum) %>% session$dataWrangler$vectorToTibble() %>% 
    dplyr::select(gene_id = variable, n_dataset = value)

  
  repro %>% 
    filter(n_dataset == max(n_dataset)) %>% 
    dplyr::select(gene_id) %>% 
    mutate(cell_type = currCelltype)
  
}) %>% session$dataWrangler$rbind()

fdata$tallyFlat <- pdata$celltypes %>% lapply(function(currCelltype) {
  fdata$tally %>% filter(cell_type == currCelltype) %>% session$dataWrangler$extractColumn("gene_id") %>% unique()
})


fdata$main <- pdata$celltypes %>% session$collectionUtils$lapply(function(celltype) {
  
  exprs <- fdata$ctprofiles %>% 
     filter(gene_id %in% fdata$tallyFlat[[celltype]], cell_type == celltype) %>% 
     dplyr::select(gene_id, dataset, expr) %>% 
     spread(dataset, expr)
  
  exprs <- exprs %>% session$dataWrangler$setColAsRownames("gene_id")
  exprs <- exprs %>% workspace$utils$logTrans()
  
  coex <- exprs %>% cor()
  
  tibble(comparison = coex %>% workspace$utils$getPairIds(), 
         cor_coef = coex %>% workspace$utils$vectorize()) %>% 
    mutate(cell_type = celltype) %>% 
    mutate(n_gene = nrow(exprs))
  
}) %>% session$dataWrangler$rbind()

fdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = cor_coef)) + 
  geom_boxplot() +
  geom_point(shape = 1, position = "jitter") + 
  ylim(-1, 1)

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_05.eps"

xdata$main <- fdata$main %>% mutate(cell_type = workspace$utils$fmtCelltypes(cell_type))

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = cor_coef)) + 
  geom_boxplot() +
  geom_point(shape = 1, position = "jitter") + 
  ylim(-1, 1) +
  session$graphingUtils$tiltX(angle = 90) + 
  xlab("Cell type") + 
  ylab("Pearson's r") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 14)) + 
  ggtitle("Correlation of expression among independent datasets")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 20, height = 15) 

# =================================================

# now you can look at correlations?

xdata <- list()

(xdata$celltype <- pdata$celltypes[6])

(xdata$main <- fdata$ctprofiles %>% 
  filter(gene_id %in% fdata$tallyFlat[[xdata$celltype]], cell_type == xdata$celltype) %>% 
  dplyr::select(gene_id, dataset, expr) %>% 
  spread(dataset, expr))

xdata$main %>% 
  session$dataWrangler$setColAsRownames("gene_id") %>% 
  ggpairs()

# ----

xdata <- list()

xdata$main <- fdata$tally %>% mutate(value = 1) %>% 
  spread(cell_type, value) %>% 
  session$dataWrangler$fillNa(colNames = pdata$celltypes, value = 0) %>% 
  session$dataWrangler$setColAsRownames("gene_id") 

xdata$main %>% UpSetR::upset(nsets = 6, text.scale = 2, keep.order = TRUE, nintersects = NA)




# ++++++++++++++
# figures
# Now this is where I'm going to look at reproducibility of the xCell networks

fdata <- list(); gc()

fdata$celltypes <- pdata$files$coexmats$cell_type %>% unique() %>% sort()
names(fdata$celltypes) <- fdata$celltypes

fdata$coexmatFiles <-  pdata$files$coexmats %>% filter(level == "sc")

fdata$celltypes %>% session$collectionUtils$lapply(function(currCelltype) {
  
  print(paste0("working on cell type: ", currCelltype, "..."))
  
  lnktbls <-fdata$coexmatFiles %>% filter(cell_type == currCelltype)
  lnktbls <- lnktbls$outfile %>% session$dataWrangler$attachNames(lnktbls$dataset)
  
  lnktbls <- lnktbls %>% mclapply(function(currFile) {
    
    lnkmat <- read_rds(paste0(workspace$outputDir, currFile))
    lnkmat[which(lnkmat >= 0.99)] <- 1
    lnkmat[which(lnkmat < 0.99)] <- 0
    
    lnktbl <- tibble(pair = lnkmat %>% workspace$utils$getPairIds(), lnk = lnkmat %>% workspace$utils$vectorize())
    
    return(lnktbl)
  }, mc.cores = length(lnktbls))
  
  print("saving to disk...")
  
  # ++++++++++++++++++++++++
  # COMMIT
  lnktbls %>% saveRDS(paste0(workspace$outputDir, paste0("sc_lnktbls_", currCelltype, ".rds")))
  
})


# +++++++++++++++++++++
# compute ORA reproducibility for every cell type, store it as one big tbl

pdata$files$lnktbls <- tibble(outfile = list.files(workspace$outputDir)) %>% 
  filter(grepl("lnktbls", outfile)) %>% 
  mutate(outfile_name = gsub(".rds", "", outfile)) %>% 
  mutate(strs = strsplit(outfile_name, "_")) %>% 
  mutate(level = strs %>% sapply(function(currStr) { currStr[1] })) %>% 
  mutate(type = strs %>% sapply(function(currStr) { currStr[2] })) %>% 
  mutate(cell_type = strs %>% sapply(function(currStr) { currStr[3] })) %>% 
  dplyr::select(outfile, level, type, cell_type)

fdata <- list(); gc()

fdata$celltypes <- pdata$files$lnktbls$cell_type %>% unique()
names(fdata$celltypes) <- fdata$celltypes

fdata$lnktblFiles <- pdata$files$lnktbls %>% filter(level == "sc")

fdata$main <- fdata$celltypes %>% session$collectionUtils$lapply(function(currCelltype) {
  
  print(paste0("working on cell type: ", currCelltype, "..."))
  
  lnkstbls <- (fdata$lnktblFiles %>% filter(cell_type == currCelltype))$outfile
  lnkstbls <- read_rds(paste0(workspace$outputDir, lnkstbls))
  
  datasets <- lnkstbls %>% names()
  names(datasets) <- datasets
  
  print(paste0("preparing lnktbls..."))
  
  lnks_T <- lnkstbls %>% lapply(function(currLnkTbl) { currLnkTbl %>% filter(lnk == 1) })
  lnks_F <- lnkstbls %>% lapply(function(currLnkTbl) { currLnkTbl %>% filter(lnk == 0) })
  
  print(paste0("running fisher's exact tests in parallel..."))
  
  datasets %>% lapply(function(currDat_T) {
    
    set_T <- lnks_T[[currDat_T]]$pair
    
    datasets_S <- datasets[datasets != currDat_T]
    
    datasets_S %>% mclapply(function(currDat_S) {
      
      session$evaluationUtils$fisher(predicted = lnks_T[[currDat_S]]$pair, 
                                     notPredicted = lnks_F[[currDat_S]]$pair, 
                                     trueSet = lnks_T[[currDat_T]]$pair, 
                                     alternative = "greater")
      
    }, mc.cores = 3) # TODO maybe change this to the number of datasets_S
    
  })
  
})

# COMMIT
fdata$main %>% saveRDS(paste0(workspace$outputDir, paste0("sc_repro_fisher.rds")))


# +++++++++++++++++++++
# finally; let's make some figures for ora reproducibility

fdata <- list()

fdata$fishers <- read_rds(paste0(workspace$outputDir, paste0("sc_repro_fisher.rds")))

fdata$main <- fdata$fishers %>% session$collectionUtils$lapplyWithName(function(currCelltype, currCTComps) {
  
  currCTComps %>% session$collectionUtils$lapplyWithName(function(currDatTruth, currComparisons) {
    
    currComparisons %>% session$collectionUtils$lapplyWithName(function(currDatSample, currComparison) {
      
      currComparison$stats %>% as.data.frame() %>% t() %>% as_tibble() %>% 
        mutate(pvalue = currComparison$test$p.value) %>% 
        mutate(dat_truth = currDatTruth, dat_sample = currDatSample, cell_type = currCelltype) %>% 
        dplyr::select(cell_type, dat_truth, dat_sample, everything())
      
    }) %>% session$dataWrangler$rbind()
    
  }) %>% session$dataWrangler$rbind()
  
}) %>% session$dataWrangler$rbind() 


# ++++++++++++++++++++++++
# COMMIT
fdata$output$fishers <- fdata$fishers
fdata$output$tbl <- fdata$main
fdata$output %>% saveRDS(paste0(workspace$outputDir, paste0("sc_repro_fisher.rds")))
# ++++++++++++++++++++++++

fdata <- list()

fdata$main <- read_rds(paste0(workspace$outputDir, paste0("sc_repro_fisher.rds")))$tbl

fdata$main <- fdata$main %>% mutate(qvalue = pvalue %>% p.adjust(method = "fdr"))


# try showing xCell reproducibility for just one dataset? velmeshev?
fdata$main %>% filter(dat_truth == "velmeshev") %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = or)) + 
  geom_boxplot(color = "grey50") +
  geom_point(alpha = 0.4, size = 5) +
  scale_y_continuous(trans = "log2") + 
  xlab("Cell type") + 
  ylab("Odds ratio") +
  ggtitle("Reproducibility of xcell co-expression signals\nin the Velmeshev dataset") +
  session$graphingUtils$tiltX(angle = 90)


# let's show the means per dataset for MSL presentation
xdata <- list()
xdata$main <- fdata$main %>% mutate(cell_type = workspace$utils$fmtCelltypes(cell_type), 
                                    dat_truth = workspace$utils$fmtDataset(dat_truth))
xdata$means <- xdata$main %>% 
  group_by(cell_type, dat_truth) %>% summarize(or = mean(or)) %>% ungroup()
xdata$means %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = or)) + 
  geom_point(aes(color = dat_truth), size = 5) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(trans = "log10") + 
  labs(x = "Cell type", 
       y = "Recovery of edges by other datasets\n(Odds ratio)", 
       color = "Dataset") +
  ggtitle("Reproducibility of cross-cell co-expression signals") +
  session$graphingUtils$tiltX(angle = 90) 
  

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_06.eps"

xdata$main <- fdata$main %>% mutate(cell_type = workspace$utils$fmtCelltypes(cell_type), 
                                    dat_truth = workspace$utils$fmtDataset(dat_truth))

xdata$means <- xdata$main %>% 
  group_by(cell_type, dat_truth) %>% summarize(or = mean(or)) %>% ungroup()

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(y = dat_truth, x = or, color = dat_truth)) + 
  geom_point(shape = 1) +
  geom_point(data = xdata$means, size = 3) +
  scale_x_continuous(trans = "log2") + 
  facet_wrap(~cell_type, ncol = 1) + 
  geom_vline(xintercept = 1, linetype = "dashed") +
  xlab("Odds ratio (log2 scale)") + 
  ggtitle("Reproducibility of xCell co-expression edges") +
  scale_color_brewer(palette = "Set1") +
  labs(color = "Dataset") +
  scale_y_discrete(limits = rev) +
  theme(axis.text = element_text(size = 13))

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 20, height = 25) 

# =================================================


# ++++++++++++++
# figures
# Now, examine the overall topolog (node degrees, etc)
# and focus on those edges that are reproducible across datasets

fdata <- list(); gc()

fdata$celltypes <- pdata$files$coexmats$cell_type %>% unique() %>% sort()
names(fdata$celltypes) <- fdata$celltypes

fdata$celltypes %>% session$collectionUtils$lapply(function(currCelltype) {
  
  print(paste0("working on cell type: ", currCelltype, "..."))
  
  lnkmats <- pdata$files$coexmats %>% filter(cell_type == currCelltype, level == "sc")
  lnkmats <- lnkmats$outfile %>% session$dataWrangler$attachNames(lnkmats$dataset)
  
  lnkmats <- lnkmats %>% mclapply(function(currFile) {
    
    lnkmat <- read_rds(paste0(workspace$outputDir, currFile))
    lnkmat[which(lnkmat >= 0.99)] <- 1
    lnkmat[which(lnkmat < 0.99)] <- 0
  
    return(lnkmat)
  }, mc.cores = length(lnkmats))
  
  print("saving to disk...")
  
  # ++++++++++++++++++++++++
  # COMMIT
  lnkmats %>% saveRDS(paste0(workspace$outputDir, paste0("sc_lnkmats_", currCelltype, ".rds")))
  
})

pdata$files$lnkmats <- tibble(outfile = list.files(workspace$outputDir)) %>% 
  filter(grepl("lnkmats", outfile)) %>% 
  mutate(outfile_name = gsub(".rds", "", outfile)) %>% 
  mutate(strs = strsplit(outfile_name, "_")) %>% 
  mutate(level = strs %>% sapply(function(currStr) { currStr[1] })) %>% 
  mutate(type = strs %>% sapply(function(currStr) { currStr[2] })) %>% 
  mutate(cell_type = strs %>% sapply(function(currStr) { currStr[3] })) %>% 
  dplyr::select(outfile, level, type, cell_type) %>% 
  filter(level == "sc")

# +++++++++++++++++++++
# compute and examine node degrees 

fdata <- list(); gc()

fdata$celltypes <- pdata$files$lnkmats$cell_type %>% unique() %>% sort()

fdata$files <- pdata$files$lnkmats$outfile
names(fdata$files) <- pdata$files$lnkmats$cell_type

fdata$degrees <- fdata$celltypes %>% mclapply(function(currCelltype) {
  
  lnkmats <- read_rds(paste0(workspace$outputDir, fdata$files[currCelltype]))
  lnkmats %>% session$collectionUtils$lapplyWithName(function(currDataset, currlnkmat) {
    currlnkmat %>% workspace$utils$computeNodeDegrees() %>% mutate(dataset = currDataset)
  }) %>% 
    session$dataWrangler$rbind() %>% 
    mutate(cell_type = currCelltype)
  
}, mc.cores = length(fdata$celltypes)) %>% session$dataWrangler$rbind()

# add in log transformed column

fdata$degrees <- fdata$degrees %>% mutate(node_degree_log10 = log10(node_degree + 1))

fdata$degrees %>% group_by(dataset, cell_type) %>% summarize(node_degree_max = max(node_degree)) %>% arrange(desc(node_degree_max)) %>% ungroup()

# in vs. out of network (having least one edge)
fdata$degrees %>% 
  mutate(in_network = (node_degree > 0)) %>% 
  group_by(dataset, cell_type, in_network) %>% summarize(n_gene = n()) %>% ungroup() %>% 
  session$graphingUtils$ggplot(aes(x = in_network, y = n_gene)) + 
  geom_bar(stat = "identity") + 
  geom_text(aes(label = n_gene), size = 3, vjust = -0.5) +
  facet_wrap(~cell_type * dataset) + 
  xlab("Number of genes in network") 

# cumulative distribution of node degrees 
# TODO --- does it matter????
fdata$degrees %>% 
  session$graphingUtils$ggplot(aes(x = node_degree_log10, color = dataset)) + 
  stat_ecdf() + 
  facet_wrap(~cell_type, ncol = 2) + 
  geom_vline(xintercept = 2, linetype = "dashed")

# correlation of node degrees
fdata$degrees %>% 
  filter(cell_type == "inhibitory") %>% 
  dplyr::select(gene_id, dataset, node_degree_log10) %>% unique() %>% 
  spread(dataset, node_degree_log10) %>% 
  na.omit() %>% 
  session$dataWrangler$setColAsRownames("gene_id") %>% 
  GGally::ggpairs()

# +++++++++++++++++++++
# also how many edges are discovered in multiple datasets, etc -- AND NOT necessary all between the same datasets?

fdata <- list(); gc()

fdata$lnktblFiles <- pdata$files$lnktbls %>% filter(level == "sc")

fdata$celltypes <- fdata$lnktblFiles$cell_type %>% unique() %>% sort() %>% session$dataWrangler$attachNames()

fdata$files <- fdata$lnktblFiles$outfile
names(fdata$files) <- fdata$lnktblFiles$cell_type

fdata$lnkSmrys <- fdata$celltypes %>% mclapply(function(celltype) {
  
  lnktbl <- readRDS(paste0(workspace$outputDir, fdata$files[celltype])) %>% 
    session$collectionUtils$lapplyWithName(function(currDataset, currLnktbl) {
      currLnktbl %>% filter(lnk > 0) %>% mutate(dataset = currDataset) 
    }) %>% session$dataWrangler$rbind()
  
  lnkSmry <- lnktbl %>% group_by(pair) %>% summarize(n_dataset = n()) %>% mutate(cell_type = celltype)
  
  return(lnkSmry)
}, mc.cores = length(fdata$celltypes))

fdata$lnkSmrys <- fdata$lnkSmrys %>% session$collectionUtils$lapply(function(tbl) {
  tbl %>% mutate(frac_dataset = n_dataset / max(n_dataset))
})

# ==== COMMIT
fdata$lnkSmrys %>% saveRDS(paste0(workspace$outputDir, "sc_lnk_smrys.rds"))
# =============

# again single vs. multiple. vs all

fdata <- list(); gc()

fdata$lnkSmrys <- read_rds(paste0(workspace$outputDir, "sc_lnk_smrys.rds"))

fdata$lnkSmrys <- fdata$lnkSmrys %>% session$collectionUtils$lapply(function(tbl) {
  tbl %>% 
    mutate(class = n_dataset %>% sapply(function(n) { min(n, 2) }) ) %>% 
    mutate(class = class + (frac_dataset %>% sapply(function(n) { if (n == 1) { n } else { 0 } }))) %>% 
    mutate(class = class %>% sapply(function(currClass) { 
      if (currClass >= 3) { "All" }
      else if (currClass >= 2) { "Multiple" }
      else { "Single" }
    }))
})

fdata$main <- fdata$lnkSmrys %>% session$collectionUtils$lapplyWithName(function(celltype, tbl) {
  tbl %>% group_by(class) %>% summarize(n_pair = n()) %>% mutate(cell_type = celltype)
}) %>% session$dataWrangler$rbind()

fdata$main$n_pair %>% sum() # number of edges in total


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_07.eps"

xdata$main <- fdata$main %>% mutate(cell_type = workspace$utils$fmtCelltypes(cell_type))

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = n_pair)) + 
  geom_bar(aes(fill = factor(class, levels = c("Single", "Multiple", "All"))), stat = "identity") + 
  session$graphingUtils$tiltX(angle = 90) + 
  scale_fill_manual(values = rev(c("grey10", "grey40", "grey80"))) + 
  scale_y_continuous(labels = label_comma()) +
  xlab("Cell type") + 
  ylab("Number of edges") + 
  labs(fill = "Datasets") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 20, height = 7) 

# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_08.eps"

xdata$main <- fdata$main %>% mutate(cell_type = workspace$utils$fmtCelltypes(cell_type))

xdata$plot <- xdata$main %>% 
  filter(class == "All") %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = n_pair)) + 
  geom_bar(color = "grey10", stat = "identity") + 
  geom_text(aes(label = n_pair), vjust = -1, size = 6) +
  session$graphingUtils$tiltX(angle = 90) + 
  scale_fill_manual(values = rev(c("grey10", "grey40", "grey80"))) + 
  scale_y_continuous(labels = label_comma(), limits = c(0, 12000)) +
  xlab("") + 
  ylab("") 

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 18, height = 12) 

# =================================================




fdata$lnktbl <- readRDS(paste0(workspace$outputDir, fdata$files[fdata$celltype])) %>% 
  session$collectionUtils$lapplyWithName(function(currDataset, currLnktbl) {
    currLnktbl %>% filter(lnk > 0) %>% mutate(dataset = currDataset) 
  }) %>% session$dataWrangler$rbind()

fdata$lnkSmry <- fdata$lnktbl %>% group_by(pair) %>% summarize(n_dataset = n())

# look at the distribution of edges by number of datasets

xdata <- list()

xdata$nLnksSmry <- fdata$lnkSmry %>% group_by(n_dataset) %>% summarize(n = n()) %>% 
  mutate(frac = n / sum(n))

xdata$nLnksSmry %>% 
  session$graphingUtils$ggplot(aes(x = n_dataset, y = frac)) + 
  geom_text(aes(label = n), vjust = -1, size = 6) +
  geom_bar(stat = "identity") + 
  ylim(0, 1) + 
  ggtitle(paste0("Distribution of network edges\nby number of independent discoveries:\n", fdata$celltype))


# ++++++++++++++++++++++=
# todo --- compute AUROC????

fdata <- list(); gc()

































