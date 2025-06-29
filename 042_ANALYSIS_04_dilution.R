pdata <- list()

pdata$figuresDir <- paste0(workspace$workspaceDir, "042_ANALYSIS_04_dilution_FIGURES/")

pdata$sfilesDir <- paste0(workspace$workspaceDir, "042_S_FILES/")

# just reuse the simulations from the previous part

pdata$this <- read_rds(paste0(workspace$outputDir, "this_cell_sampling.rds"))

pdata$exprdat <- read_rds(paste0(workspace$outputDir, "exprmats_agg_cell_sampling.rds"))
pdata$exprmats <- pdata$exprdat$exprmats

# +++++++++++++++++++++++++++++++++++++++============
# compute the diluted bulk tissue data

fdata <- list()

fdata$celltypes <- c("excitatory", "inhibitory", "oligodendrocyte", "astrocyte", "microglia", "opc") %>% session$dataWrangler$attachNames()
fdata$nCels <- c(9750, 3250, 8775, 5850, 1950, 2925) %>% session$dataWrangler$attachNames(useNames = fdata$celltypes)

fdata$fracCells <- fdata$nCels / sum(fdata$nCels)

fdata$exprmat$diluted <- Reduce("+", fdata$celltypes %>% lapply(function(currCelltype) {
  pdata$exprmats$sbjv[[currCelltype]] * fdata$fracCells[currCelltype]
}))

fdata$exprmat$sbjv <- pdata$exprmats$sbjv

# COMMIT =====
fdata$exprmat %>% saveRDS(paste0(workspace$outputDir, "exprmats_dilution.rds"))
# ============

# SOFT COMMIT ===
pdata$exprmats$diluted <- read_rds(paste0(workspace$outputDir, "exprmats_dilution.rds"))
# ======================

# +++++++++++++=====================
# show a "schematic" that relates all the cells to one 

# VWA3A ; random genes

fdata <- list(); gc()

# 1 TSPAN7   1186.
# 2 CCDC88A  1127.
# 3 HNRNPDL   847

fdata$gene <- "TSPAN7"

fdata$exprmats$sbjv <- pdata$exprmats$diluted$sbjv
fdata$exprmats$diluted <- pdata$exprmats$diluted$diluted

fdata$exprs <- pdata$exprmats$sbjv %>% session$collectionUtils$lapplyWithName(function(celltype, exprmat) {
  exprmat[fdata$gene, ] %>% session$dataWrangler$vectorToTibble() %>% 
    dplyr::select(subject = variable, expr = value) %>% 
    mutate(cell_type = celltype)
}) %>% session$dataWrangler$rbind()

fdata$exprs %>% 
  session$graphingUtils$ggplot(aes(x = subject, y = expr)) + 
  geom_point() + 
  geom_line(group = 1) +
  facet_wrap(~cell_type, ncol = 1, scales = "free_y") + 
  session$graphingUtils$tiltX(angle = 90) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

fdata$exprmats$diluted[fdata$gene, ] %>% session$dataWrangler$vectorToTibble() %>% 
  dplyr::select(subject = variable, expr = value) %>% 
  session$graphingUtils$ggplot(aes(x = subject, y = expr)) + 
  geom_point() + 
  geom_line(group = 1) +
  session$graphingUtils$tiltX(angle = 90) + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_01.eps" # vector file too big just draw bitmap

xdata$main <- fdata$exprs %>% mutate(cell_type = workspace$utils$fmtCelltypes(cell_type))

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = subject, y = expr)) + 
  geom_point() + 
  geom_line(group = 1) +
  facet_wrap(~cell_type, ncol = 1, scales = "free_y") + 
  session$graphingUtils$tiltX(angle = 90) + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank()) 

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 15, height = 20) 
# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_02.eps" # vector file too big just draw bitmap

xdata$main <- fdata$exprmats$diluted[fdata$gene, ] 

xdata$plot <- xdata$main %>% 
  session$dataWrangler$vectorToTibble() %>% 
  dplyr::select(subject = variable, expr = value) %>% 
  session$graphingUtils$ggplot(aes(x = subject, y = expr)) + 
  geom_point() + 
  geom_line(group = 1) +
  session$graphingUtils$tiltX(angle = 90) + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank()) +
  ggtitle("Bulk")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 15, height = 4) 
# =================================================

# +++++++++++++++++++++++++++++++++++++++============
# compute rsqr for diluted

fdata <- list()

fdata$exprmats <- read_rds(paste0(workspace$outputDir, "exprmats_dilution.rds"))

fdata$rsqrs <- workspace$utils$computeRSqr(exprmatObs = fdata$exprmats$diluted, exprmatCovs = fdata$exprmats$sbjv)

fdata$rsqrs %>% summary()

fdata$rsqrstbl <- fdata$rsqrs %>% gather(cell_type, rsqr, -gene_id) 

# COMMIT =====================
fdata$rsqrstbl %>% saveRDS(paste0(workspace$outputDir, "rsqrs_dilution.rds"))
# ===============================

fdata <- list()

fdata$rsqrstbl <- read_rds(paste0(workspace$outputDir, "rsqrs_dilution.rds"))

fdata$rsqrstbl %>% group_by(gene_id) %>% summarize(rsqr_sum = sum(rsqr)) %>% 
  session$graphingUtils$ggplot(aes(x = rsqr_sum)) + geom_density()


fdata$rsqrGene <- fdata$rsqrstbl %>% filter(gene_id == "DBNDD2")


fdata$rsqrstbl %>% group_by(cell_type) %>% summarize(rsqr = mean(rsqr))

# isolate cases of "winner takes all"
fdata$rsqrstbl %>% filter(rsqr >= 0.8) %>% group_by(cell_type) %>% summarize(n = n()) # show extension to co-expression & elevated expression in these genes


fdata$rsqrstbl %>% filter(rsqr >= 0.8) %>% group_by(cell_type) %>% summarize(n = n()) -> x; x$n %>% sum()

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_03.eps"

xdata$main <- fdata$rsqrstbl %>% 
  mutate(cell_type = strsplit(cell_type, "_")) %>% 
  mutate(cell_type = cell_type %>% sapply(function(currStr) { currStr[2] })) %>% 
  mutate(cell_type = cell_type %>% workspace$utils$fmtCelltypes())

xdata$means <- xdata$main %>% group_by(cell_type) %>% summarize(rsqr = mean(rsqr))

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = rsqr)) + 
  geom_histogram(aes(fill = cell_type)) + 
  geom_vline(data = xdata$means, aes(xintercept = rsqr, color = cell_type)) + 
  geom_vline(xintercept = 0.8, linewidth = 0.1) + 
  facet_wrap(~cell_type, ncol = 1) + 
  scale_fill_brewer(palette = "Set1") + 
  scale_color_brewer(palette = "Set1") + 
  theme(legend.position = "none") + 
  xlab(bquote("PVE by cell type specific expression"~(R^2))) + 
  ylab("Number of genes") + 
  ggtitle("Dilution of cell type specific expression variability") + 
  theme(plot.title = element_text(size = 25), plot.subtitle = element_text(size = 14))  
  
xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 30, height = 30) 
# =================================================

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_04.eps"

xdata$main <- fdata$rsqrstbl %>% 
  mutate(cell_type = strsplit(cell_type, "_")) %>% 
  mutate(cell_type = cell_type %>% sapply(function(currStr) { currStr[2] })) %>% 
  mutate(cell_type = cell_type %>% workspace$utils$fmtCelltypes())

xdata$main <- xdata$main %>% group_by(gene_id) %>% filter(rsqr == max(rsqr)) %>% ungroup()

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = rsqr)) + 
  geom_histogram(aes(fill = cell_type)) + 
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 18)) + 
  ggtitle("Distribution of PVE by the maximal cell type across genes") + 
  xlab(bquote(R^2)) + 
  ylab("Number of genes") +
  xlim(0, 1)

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 18, height = 15) 
# =================================================


# ++++++++++++++++==+++++++++++=======
# show that those with higher rsqr are those with higher expression
# TODO - possibly make this for supplementary

fdata <- list(); gc()

fdata$rsqrtbl <- read_rds(paste0(workspace$outputDir, "rsqrs_dilution.rds"))

fdata$rsqrtbl <- fdata$rsqrtbl %>% 
  group_by(cell_type) %>% 
  mutate(rank = rank(rsqr)) %>%
  mutate(rank = rank / max(rank)) %>% 
  ungroup()

fdata$rsqrtbl %>% filter(rsqr > 0.9) %>% group_by(cell_type) %>% summarize(n_gene = n()) # actually focus on these!!!

fdata$rsqrTop <- fdata$rsqrtbl %>% filter(rank > 0.99) 

fdata$cteprf <- pdata$this$refdat$cteprf[fdata$rsqrTop$gene_id, ] 

fdata$cteprf <- log2(fdata$cteprf + 1)

fdata$cteprf <- fdata$cteprf %>% apply(1, function(vals) { (vals - mean(vals)) / sd(vals) }) %>% t()

fdata$cteprf %>% session$graphingUtils$heatmap(cluster_rows = FALSE, cluster_cols = FALSE, 
                                               show_rownames = FALSE)



# ++++++++++++++++++============
# let's show it clearly why GAD2 is dominated entirely by inhibitory neurons

# PDGFB
# BBS2

fdata <- list(); gc()

fdata$gene <- "BBS2"

fdata$exprmats$sbjv <- pdata$exprmats$diluted$sbjv
fdata$exprmats$diluted <- pdata$exprmats$diluted$diluted

fdata$exprs <- pdata$exprmats$sbjv %>% session$collectionUtils$lapplyWithName(function(celltype, exprmat) {
  exprmat[fdata$gene, ] %>% session$dataWrangler$vectorToTibble() %>% 
    dplyr::select(subject = variable, expr = value) %>% 
    mutate(cell_type = celltype)
}) %>% session$dataWrangler$rbind()

fdata$exprsBulk <- fdata$exprmats$diluted[fdata$gene, ] %>% session$dataWrangler$vectorToTibble() %>% 
  dplyr::select(subject = variable, expr = value) %>% 
  mutate(expr_log = log2(expr + 1))

fdata$exprs %>% 
  mutate(expr_log = log2(expr + 1)) %>% 
  session$graphingUtils$ggplot(aes(x = subject, y = expr_log)) + 
  geom_point(aes(color = cell_type)) + 
  geom_line(aes(color = cell_type, group = cell_type)) +
  geom_point(data = fdata$exprsBulk) +
  geom_line(data = fdata$exprsBulk, group = 1, size = 1) +
  session$graphingUtils$tiltX(angle = 90) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "none") + 
  scale_color_brewer(palette = "Set1") 


# ++++++++++++++++++============
# example of gene dominated by microglia: PDGFB

fdata <- list(); gc()

fdata$gene <- "PDGFB"

fdata$exprmats$sbjv <- pdata$exprmats$diluted$sbjv
fdata$exprmats$diluted <- pdata$exprmats$diluted$diluted

fdata$exprs <- pdata$exprmats$sbjv %>% session$collectionUtils$lapplyWithName(function(celltype, exprmat) {
  exprmat[fdata$gene, ] %>% session$dataWrangler$vectorToTibble() %>% 
    dplyr::select(subject = variable, expr = value) %>% 
    mutate(cell_type = celltype)
}) %>% session$dataWrangler$rbind()

fdata$exprsBulk <- fdata$exprmats$diluted[fdata$gene, ] %>% session$dataWrangler$vectorToTibble() %>% 
  dplyr::select(subject = variable, expr = value) %>% 
  mutate(expr_log = log2(expr + 1))

fdata$exprs %>% filter(cell_type == "microglia") %>% 
  left_join(fdata$exprsBulk %>% dplyr::select(subject, expr_log)) %>% 
  summarize(cor(expr, expr_log))


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_05.eps"

xdata$plot <- fdata$exprs %>% 
  mutate(expr_log = log2(expr + 1)) %>% 
  session$graphingUtils$ggplot(aes(x = subject, y = expr_log)) + 
  geom_point(aes(color = cell_type)) + 
  geom_line(aes(color = cell_type, group = cell_type)) +
  geom_point(data = fdata$exprsBulk, size = 2) +
  geom_line(data = fdata$exprsBulk, group = 1, linewidth = 1) +
  session$graphingUtils$tiltX(angle = 90) + 
  theme(axis.text.x = element_blank(), 
        legend.position = "none") + 
  scale_color_brewer(palette = "Set1")  + 
  ggtitle(paste0("Expression of sim-", fdata$gene)) + 
  ylab("log2(CPM)") + 
  xlab("Subject")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 18, height = 10) 
# =================================================


# ++++++++++++++++++============
# gene where no cell type dominates: BBS2

fdata <- list(); gc()

fdata$gene <- "BBS2"

fdata$exprmats$sbjv <- pdata$exprmats$diluted$sbjv
fdata$exprmats$diluted <- pdata$exprmats$diluted$diluted

fdata$exprs <- pdata$exprmats$sbjv %>% session$collectionUtils$lapplyWithName(function(celltype, exprmat) {
  exprmat[fdata$gene, ] %>% session$dataWrangler$vectorToTibble() %>% 
    dplyr::select(subject = variable, expr = value) %>% 
    mutate(cell_type = celltype)
}) %>% session$dataWrangler$rbind()

fdata$exprsBulk <- fdata$exprmats$diluted[fdata$gene, ] %>% session$dataWrangler$vectorToTibble() %>% 
  dplyr::select(subject = variable, expr = value) %>% 
  mutate(expr_log = log2(expr + 1))

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_06.eps"

xdata$plot <- fdata$exprs %>% 
  mutate(expr_log = log2(expr + 1)) %>% 
  session$graphingUtils$ggplot(aes(x = subject, y = expr_log)) + 
  geom_point(aes(color = cell_type)) + 
  geom_line(aes(color = cell_type, group = cell_type)) +
  geom_point(data = fdata$exprsBulk, size = 2) +
  geom_line(data = fdata$exprsBulk, group = 1, size = 1) +
  session$graphingUtils$tiltX(angle = 90) + 
  theme(axis.text.x = element_blank(), 
        legend.position = "none") + 
  scale_color_brewer(palette = "Set1")  + 
  ggtitle(paste0("Expression of sim-", fdata$gene)) + 
  ylab("log2(CPM)") + 
  xlab("Subject")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 18, height = 10) 
# =================================================














# ++++++++++++++++++++++===============================
# ++++++++++++++++++================================
# let's make simulate a version every all 3k genes are synchronized across 30 subjects, for all cell types 
# and co-expressed at 0.8? 

pdata$utils$simulate <- function(this, nCels, nSbjs) {
  # output a sim object
  
  cellsBsln <- this %>% workspace$api$simBseLnCels(nSubject = nSbjs, nCell = nCels) 
  
  sbjMeans <- this %>% workspace$api$simSbjLvMeans(nSbjs)
  
  celParamsSbjv <- this %>% workspace$api$computeCelParams(sbjMeans$simSbjs) 
  
  cellsSbjv <- cellsBsln$simCells %>% workspace$api$convertCelLvDist(cellsBsln$params$md, celParamsSbjv)
  
  # =================================  
  simObjs <- list()
  simObjs$cellsBsln <- list(simdat = cellsBsln$simCells, params = cellsBsln$params)
  simObjs$cellsSbjv <- list(simdat = cellsSbjv, params = celParamsSbjv)
  simObjs$subjects <- list(simdat = sbjMeans$simSbjs, params = sbjMeans$params)
  # =================================
  
  return(simObjs)  
} 

pdata$utils$xSbjExprmats <- function(simObjs) {
  # return xSbj from 1. bsln, 2. sbjv, 3. sbjmeans
  celltypes <- c("excitatory", "inhibitory", "oligodendrocyte", "astrocyte", "microglia", "opc") %>% session$dataWrangler$attachNames()
  ncels <- c(9750, 3250, 8775, 5850, 1950, 2925) %>% session$dataWrangler$attachNames(useNames = celltypes)
  
  sbjs <- simObjs$cellsBsln$simdat$samples$subject %>% unique() %>% sort() %>% session$dataWrangler$attachNames()
  
  cells <- celltypes %>% session$collectionUtils$lapply(function(celltype) {
    sbjs %>% session$collectionUtils$lapply(function(currSbj) { 
      currSamples <- simObjs$cellsBsln$simdat$samples %>% filter(subject == currSbj, cell_type == celltype)
      bootstrap <- currSamples$sample %>% sample(ncels[celltype], replace = TRUE)
      currSamples <- tibble(sample = bootstrap) %>% left_join(currSamples, by = "sample")
      return(currSamples)
    }, verbose = FALSE) %>% session$dataWrangler$rbind()
  })
  
  exprmats <- list()
  
  exprmats$bsln <- cells %>% session$collectionUtils$lapply(function(currCells) {
    workspace$utils$aggregSbj(simObjs$cellsBsln$simdat$exprmat, samples = currCells)
  })
  
  # exprmats$bsln %>% lapply(function(exprmat) { exprmat %>% t() %>% cor() %>% mean() }) # check coexpression 
  
  exprmats$sbjv <- cells %>% session$collectionUtils$lapply(function(currCells) {
    workspace$utils$aggregSbj(simObjs$cellsSbjv$simdat$exprmat, samples = currCells)
  })
  
  # exprmats$sbjv %>% lapply(function(exprmat) { exprmat %>% t() %>% cor() %>% mean() }) # check coexpression
  
  exprmats$islv <- celltypes %>% session$collectionUtils$lapply(function(celltype) {
    samples <- simObjs$subjects$simdat$samples %>% filter(cell_type == celltype)
    exprmat <- simObjs$subjects$simdat$exprmat[, samples$sample]
    colnames(exprmat) <- samples$subject
    return(exprmat)
  })
  
  # exprmats$islv %>% lapply(function(exprmat) { exprmat %>% t() %>% cor() %>% mean() }) # check coexpression
  
  output <- list()
  output$cells <- cells 
  output$exprmats <- exprmats
  return(output)
}
  
     
# in scenario 4, co-expressed and synchronized
xdata <- list(); gc()
xdata$genes <- rownames(pdata$this$mdParams$sbj$astrocyte)
xdata$coexPrgms$sync <- list(celltypes = c("excitatory", "inhibitory", "astrocyte", "microglia", "oligodendrocyte", "opc"), 
                             genes = xdata$genes, value = 0.8)
xdata$this <- pdata$this %>% workspace$api$setRefDatCoexPrgms(xdata$coexPrgms, level = "cel") 
xdata$this <- xdata$this %>% workspace$api$setRefDatCoexPrgms(xdata$coexPrgms, level = "sbj") 
xdata$this$mdParams$cel <- xdata$this$mdParams$cel %>% lapply(function(mdParam) { mdParam[xdata$genes, ] })
xdata$this$mdParams$sbj <- xdata$this$mdParams$sbj %>% lapply(function(mdParam) { mdParam[xdata$genes, ] })
xdata$simObjs <- xdata$this %>% pdata$utils$simulate(5000, 30)
# COMMIT ==========================
xdata$simObjs %>% saveRDS(paste0(workspace$outputDir, "simObjs_dilution_ct_sync.rds")) 
# =================================     

xdata$simObjs$cellsSbjv$simdat$exprmat %>% dim()

# ++++++++++++++++==============================
# load in the synchrony simbObjs and save into xSubject and baseline bulk data expression matrices 
fdata <- list()

fdata$simObjs <- readRDS(paste0(workspace$outputDir, "simObjs_dilution_ct_sync.rds"))

fdata$exprdat$xsbj <- pdata$utils$xSbjExprmats(fdata$simObjs)

fdata$exprmats$sbjv <- fdata$exprdat$xsbj$exprmats$sbjv

fdata$celltypes <- c("excitatory", "inhibitory", "oligodendrocyte", "astrocyte", "microglia", "opc") %>% session$dataWrangler$attachNames()
fdata$nCels <- c(9750, 3250, 8775, 5850, 1950, 2925) %>% session$dataWrangler$attachNames(useNames = fdata$celltypes)

fdata$fracCells <- fdata$nCels / sum(fdata$nCels)

fdata$exprmat$diluted <- Reduce("+", fdata$celltypes %>% lapply(function(currCelltype) {
  fdata$exprmats$sbjv[[currCelltype]] * fdata$fracCells[currCelltype]
}))

fdata$exprmat$sbjv <- fdata$exprmats$sbjv

# COMMIT =====
fdata$exprmat %>% saveRDS(paste0(workspace$outputDir, "exprmats_dilution_sync.rds"))
# ============

   
# +++++++++++++++++++++++++++++++++++++++============
# compute rsqr for diluted

fdata <- list()

fdata$exprmats <- read_rds(paste0(workspace$outputDir, "exprmats_dilution_sync.rds"))

fdata$rsqrs <- workspace$utils$computeRSqr(exprmatObs = fdata$exprmats$diluted, exprmatCovs = fdata$exprmats$sbjv)

fdata$rsqrs %>% summary()

fdata$rsqrstbl <- fdata$rsqrs %>% gather(cell_type, rsqr, -gene_id) 

# COMMIT =====================
fdata$rsqrstbl %>% saveRDS(paste0(workspace$outputDir, "rsqrs_dilution_sync.rds"))
# ===============================

fdata <- list()

fdata$rsqrstbl <- read_rds(paste0(workspace$outputDir, "rsqrs_dilution_sync.rds"))

fdata$rsqrstbl %>% group_by(gene_id) %>% summarize(rsqr_sum = sum(rsqr)) %>% 
  session$graphingUtils$ggplot(aes(x = rsqr_sum)) + geom_density()


fdata$rsqrGene <- fdata$rsqrstbl %>% filter(gene_id == "DBNDD2")


fdata$rsqrstbl %>% group_by(cell_type) %>% summarize(rsqr = mean(rsqr))

# isolate cases of "winner takes all"
fdata$rsqrstbl %>% filter(rsqr > 0.8) %>% group_by(cell_type) %>% summarize(n = n()) # show extension to co-expression & elevated expression in these genes


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_07.eps"

xdata$main <- fdata$rsqrstbl %>% 
  mutate(cell_type = strsplit(cell_type, "_")) %>% 
  mutate(cell_type = cell_type %>% sapply(function(currStr) { currStr[2] })) %>% 
  mutate(cell_type = cell_type %>% workspace$utils$fmtCelltypes())

xdata$means <- xdata$main %>% group_by(cell_type) %>% summarize(rsqr = mean(rsqr))

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = rsqr)) + 
  geom_histogram(aes(fill = cell_type)) + 
  geom_vline(data = xdata$means, aes(xintercept = rsqr, color = cell_type)) + 
  geom_vline(xintercept = 0.8, linewidth = 0.1) + 
  facet_wrap(~cell_type, ncol = 1) + 
  scale_fill_brewer(palette = "Set1") + 
  scale_color_brewer(palette = "Set1") + 
  theme(legend.position = "none") + 
  xlab(bquote("PVE by cell type specific expression"~(R^2))) + 
  ylab("Number of genes") + 
  ggtitle("Dilution of cell type specific expression variability\nfor genes with inter-cell type synchrony") + 
  theme(plot.title = element_text(size = 25), plot.subtitle = element_text(size = 14))  

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 30, height = 30) 
# =================================================

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_08.eps"

xdata$main <- fdata$rsqrstbl %>% 
  mutate(cell_type = strsplit(cell_type, "_")) %>% 
  mutate(cell_type = cell_type %>% sapply(function(currStr) { currStr[2] })) %>% 
  mutate(cell_type = cell_type %>% workspace$utils$fmtCelltypes())

xdata$main <- xdata$main %>% group_by(gene_id) %>% filter(rsqr == max(rsqr)) %>% ungroup()

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = rsqr)) + 
  geom_histogram(aes(fill = cell_type)) + 
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 18)) + 
  ggtitle("Distribution of PVE by the maximal cell type across genes\nwith inter-cell type synchrony") + 
  xlab(bquote(R^2)) + 
  ylab("Number of genes") +
  xlim(0, 1)

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 18, height = 15) 
# =================================================


# ++++++++++++++++++============
# example of gene dominated by microglia: PDGFB

fdata <- list(); gc()

fdata$exprdat <- read_rds(paste0(workspace$outputDir, "exprmats_dilution_sync.rds"))

fdata$gene <- "PDGFB"

fdata$exprmats$sbjv <- fdata$exprdat$sbjv
fdata$exprmats$diluted <- fdata$exprdat$diluted

fdata$exprs <- fdata$exprmats$sbjv %>% session$collectionUtils$lapplyWithName(function(celltype, exprmat) {
  exprmat[fdata$gene, ] %>% session$dataWrangler$vectorToTibble() %>% 
    dplyr::select(subject = variable, expr = value) %>% 
    mutate(cell_type = celltype)
}) %>% session$dataWrangler$rbind()

fdata$exprsBulk <- fdata$exprmats$diluted[fdata$gene, ] %>% session$dataWrangler$vectorToTibble() %>%
  dplyr::select(subject = variable, expr = value) %>%
  mutate(expr_log = log2(expr + 1))

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_09.eps"

xdata$plot <- fdata$exprs %>% 
  mutate(expr_log = log2(expr + 1)) %>% 
  session$graphingUtils$ggplot(aes(x = subject, y = expr_log)) + 
  geom_point(aes(color = cell_type)) + 
  geom_line(aes(color = cell_type, group = cell_type)) +
  geom_point(data = fdata$exprsBulk, size = 2) +
  geom_line(data = fdata$exprsBulk, group = 1, size = 1) +
  session$graphingUtils$tiltX(angle = 90) + 
  theme(axis.text.x = element_blank(), 
        legend.position = "none") + 
  scale_color_brewer(palette = "Set1")  + 
  ggtitle(paste0("Expression of sim-", fdata$gene)) + 
  ylab("log2(CPM)") + 
  xlab("Subject")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 18, height = 10) 
# =================================================

# ++++++++++++++++++============
# gene where no cell type dominates: BBS2

fdata <- list(); gc()

fdata$exprdat <- read_rds(paste0(workspace$outputDir, "exprmats_dilution_sync.rds"))

fdata$gene <- "BBS2"

fdata$exprmats$sbjv <- fdata$exprdat$sbjv
fdata$exprmats$diluted <- fdata$exprdat$diluted

fdata$exprs <- fdata$exprmats$sbjv %>% session$collectionUtils$lapplyWithName(function(celltype, exprmat) {
  exprmat[fdata$gene, ] %>% session$dataWrangler$vectorToTibble() %>% 
    dplyr::select(subject = variable, expr = value) %>% 
    mutate(cell_type = celltype)
}) %>% session$dataWrangler$rbind()

fdata$exprsBulk <- fdata$exprmats$diluted[fdata$gene, ] %>% session$dataWrangler$vectorToTibble() %>%
  dplyr::select(subject = variable, expr = value) %>%
  mutate(expr_log = log2(expr + 1))

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_10.eps"

xdata$plot <- fdata$exprs %>% 
  mutate(expr_log = log2(expr + 1)) %>% 
  session$graphingUtils$ggplot(aes(x = subject, y = expr_log)) + 
  geom_point(aes(color = cell_type)) + 
  geom_line(aes(color = cell_type, group = cell_type)) +
  geom_point(data = fdata$exprsBulk, size = 2) +
  geom_line(data = fdata$exprsBulk, group = 1, size = 1) +
  session$graphingUtils$tiltX(angle = 90) + 
  theme(axis.text.x = element_blank(), 
        legend.position = "none") + 
  scale_color_brewer(palette = "Set1")  + 
  ggtitle(paste0("Expression of sim-", fdata$gene)) + 
  ylab("log2(CPM)") + 
  xlab("Subject")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 18, height = 10) 
# =================================================


# =================================================
# OUTPUT S_FILE
xdata <- list()

xdata$fileName <- "sfile_05_dilution_effects.tsv" 

xdata$rsqr <- read_rds(paste0(workspace$outputDir, "rsqrs_dilution.rds"))
xdata$rsqrSync <- read_rds(paste0(workspace$outputDir, "rsqrs_dilution_sync.rds"))

xdata$main <- xdata$rsqr %>% mutate(synchronized = "sync_not") %>% 
  rbind(xdata$rsqrSync %>% mutate(synchronized = "sync")) 

xdata$main <- xdata$main %>% spread(synchronized, rsqr)

write.table(xdata$main, file = paste0(pdata$sfilesDir, xdata$fileName), 
            sep = "\t", row.names = FALSE) 
# =================================================















