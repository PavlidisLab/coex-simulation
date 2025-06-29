pdata <- list()

pdata$figuresDir <- paste0(workspace$workspaceDir, "042_ANALYSIS_03_cell_sampling_FIGURES/")
pdata$sfilesDir <- paste0(workspace$workspaceDir, "042_S_FILES/")

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++===========
# use the demo this and just change the coex programs to none
# no co-expressions in this dataset

fdata <- list(); gc()

fdata$this <- read_rds(paste0(workspace$outputDir, "this_demo.rds"))

# no coex program defined for the first demo version
fdata$this <- fdata$this %>% workspace$api$setRefDatCoexPrgms(list(), level = "sbj") 
fdata$this <- fdata$this %>% workspace$api$setRefDatCoexPrgms(list(), level = "cel") 

# =================================  
# COMMIT ==========================
fdata$this %>% saveRDS(paste0(workspace$outputDir, "this_cell_sampling.rds"))
# =================================  

# ++++++++++++++++++
# now simulate cells
# in this case you actually need to generate enough cells per cell type 
# go with a large superset to work with, maybe 5k cells per cell type then sample with replacement later
# just 15 samples? 

fdata <- list(); gc()

fdata$this <- read_rds(paste0(workspace$outputDir, "this_cell_sampling.rds"))

# "MATK"         "DLX6-AS1"     "LOC101927745" "PHLDB2"       "TRIM54"       "EMBP1"       

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

fdata$simObjs %>% saveRDS(paste0(workspace$outputDir, "simObjs_cell_sampling.rds"))
# ================================= 

# SOFT COMMIT ++++++++++++++++++++++++++++++++++============
pdata$this <- read_rds(paste0(workspace$outputDir, "this_cell_sampling.rds"))
pdata$simObjs <- read_rds(paste0(workspace$outputDir, "simObjs_cell_sampling.rds"))
# ++++++++++++++++++++++++++++++++++++++=
 
# +++++++++++++++++++++++++++++++++++++++============
# start with excitatory neurons: 100 cells

fdata <- list(); gc()

fdata$simObj <- pdata$simObjs$cellsBsln
fdata$samples <- fdata$simObj$simdat$samples

fdata$celltypes <- c("excitatory") %>% session$dataWrangler$attachNames()
fdata$nCels <- c(100) %>% session$dataWrangler$attachNames(useNames = fdata$celltypes)

fdata$sbjs <- pdata$simObjs$cellsBsln$simdat$samples$subject %>% unique() %>% sort() %>% session$dataWrangler$attachNames()

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
  workspace$utils$aggregSbj(pdata$simObjs$cellsBsln$simdat$exprmat, samples = currCells)
})

fdata$exprmats$sbjv <- fdata$cells %>% session$collectionUtils$lapply(function(currCells) {
  workspace$utils$aggregSbj(pdata$simObjs$cellsSbjv$simdat$exprmat, samples = currCells)
})

fdata$exprmats$islv <- fdata$celltypes %>% session$collectionUtils$lapply(function(celltype) {
  samples <- pdata$simObjs$subjects$simdat$samples %>% filter(cell_type == celltype)
  exprmat <- pdata$simObjs$subjects$simdat$exprmat[, samples$sample]
  colnames(exprmat) <- samples$subject
  return(exprmat)
})

# COMMIT =====
fdata$output <- list()
fdata$output$cells <- fdata$cells 
fdata$output$exprmats <- fdata$exprmats
fdata$output %>% saveRDS(paste0(workspace$outputDir, "exprmats_agg_cell_sampling_100_cells.rds"))
# ============


# +++++++++++++++++++++++++++++++++++++++============
# start with excitatory neurons: 1000 cells

fdata <- list(); gc()

fdata$simObj <- pdata$simObjs$cellsBsln
fdata$samples <- fdata$simObj$simdat$samples

fdata$celltypes <- c("excitatory") %>% session$dataWrangler$attachNames()
fdata$nCels <- c(1000) %>% session$dataWrangler$attachNames(useNames = fdata$celltypes)

fdata$sbjs <- pdata$simObjs$cellsBsln$simdat$samples$subject %>% unique() %>% sort() %>% session$dataWrangler$attachNames()

set.seed(0); fdata$cells <-  fdata$celltypes %>% session$collectionUtils$lapply(function(celltype) {
  fdata$sbjs %>% session$collectionUtils$lapply(function(currSbj) { 
    currSamples <- fdata$samples %>% filter(subject == currSbj, cell_type == celltype)
    bootstrap <- currSamples$sample %>% sample(fdata$nCels[celltype], replace = TRUE)
    currSamples <- tibble(sample = bootstrap) %>% left_join(currSamples, by = "sample")
    return(currSamples)
  }, verbose = FALSE) %>% session$dataWrangler$rbind()
})

# now compute the sample level expressions

fdata$exprmats$bsln <- fdata$cells %>% session$collectionUtils$lapply(function(currCells) {
  workspace$utils$aggregSbj(pdata$simObjs$cellsBsln$simdat$exprmat, samples = currCells)
})

fdata$exprmats$sbjv <- fdata$cells %>% session$collectionUtils$lapply(function(currCells) {
  workspace$utils$aggregSbj(pdata$simObjs$cellsSbjv$simdat$exprmat, samples = currCells)
})

fdata$exprmats$islv <- fdata$celltypes %>% session$collectionUtils$lapply(function(celltype) {
  samples <- pdata$simObjs$subjects$simdat$samples %>% filter(cell_type == celltype)
  exprmat <- pdata$simObjs$subjects$simdat$exprmat[, samples$sample]
  colnames(exprmat) <- samples$subject
  return(exprmat)
})

# COMMIT =====
fdata$output <- list()
fdata$output$cells <- fdata$cells 
fdata$output$exprmats <- fdata$exprmats
fdata$output %>% saveRDS(paste0(workspace$outputDir, "exprmats_agg_cell_sampling_1000_cells.rds"))
# ============


# +++++++++++++++++++++++++=+++==========
# plot the basic things with 100 cells to show the ratios between bsln variance and isv variance

fdata <- list()

fdata$main <- read_rds(paste0(workspace$outputDir, "exprmats_agg_cell_sampling_100_cells.rds"))

fdata$exprmats <- fdata$main$exprmats
fdata$exprmats <- fdata$exprmats %>% lapply(function(currExprmats) { currExprmats[["excitatory"]] })

fdata$rsqrs <- workspace$utils$computeRSqr(exprmatObs = fdata$exprmats$sbjv, exprmatCovs = fdata$exprmats[c("bsln", "islv")])

p <- fdata$rsqrs %>% session$graphingUtils$ggplot(aes(x = rsqr_islv , y = rsqr_bsln)) + 
  geom_point(shape = 1) + 
  geom_point(data = fdata$rsqrs %>% filter(gene_id == "GINM1"), color = "red", size = 3) +
  xlim(0, 1) + 
  ylim(0, 1) 

# p %>% ggMarginal(type = "histogram", size = 3)

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_01.png" # vector file too big just draw bitmap

xdata$means <- fdata$rsqrs %>% summarize(rsqr_bsln = mean(rsqr_bsln), rsqr_islv = mean(rsqr_islv))
# > xdata$means 
# # A tibble: 1 × 2
# rsqr_bsln rsqr_islv
# <dbl>     <dbl>
#   1     0.437     0.562

xdata$plot <- fdata$rsqrs %>% session$graphingUtils$ggplot(aes(x = rsqr_islv , y = rsqr_bsln)) + 
  geom_point(shape = 1) + 
  geom_point(data = fdata$rsqrs %>% filter(gene_id == "GINM1"), color = "red", size = 3) +
  geom_vline(data = xdata$means, aes(xintercept = rsqr_islv), linetype = "dashed", color = "red") +
  geom_hline(data = xdata$means, aes(yintercept = rsqr_bsln), linetype = "dashed", color = "red") +
  xlim(0, 1) + 
  ylim(0, 1) + 
  xlab(bquote("PVE by ISV"~(R^2))) +
  ylab(bquote("PVE by CSIV"~(R^2))) + 
  ggtitle("CSIV vs. ISV for 100 cells per sample") +
  theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 16))

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 15, height = 15) 
# =================================================


# +++++++++++++++++++++++++=+++==========
# now lets explore the dynamics of these genes expressions 
# first at baseline (only sampling)

fdata <- list(); gc() 

fdata$main <- read_rds(paste0(workspace$outputDir, "exprmats_agg_cell_sampling_100_cells.rds"))

fdata$exprmats <- fdata$main$exprmats
fdata$exprmats <- fdata$exprmats %>% lapply(function(currExprmats) { currExprmats[["excitatory"]] })
fdata$cells <- fdata$main$cells$excitatory

fdata$gene <- "GINM1" 

# first let's look at the expression of this gene at baseline
fdata$exprGene <- pdata$simObjs$cellsBsln$simdat$exprmat[fdata$gene, fdata$cells$sample]
fdata$exprGene <- fdata$exprGene %>% 
  session$dataWrangler$vectorToTibble() %>% 
  dplyr::select(sample = variable, expr = value) %>% 
  unique()

fdata$exprGene <- fdata$cells %>% left_join(fdata$exprGene, by = "sample")

# put together the bsln stuff
fdata$mm$bsln <- pdata$simObjs$cellsBsln$params$md$excitatory %>% 
  workspace$utils$computeGamMM() %>% 
  session$dataWrangler$setRownameAsColumn("gene")

fdata$mm$bsln <- fdata$mm$bsln %>%
  mutate(sd = sqrt(variance)) %>% 
  mutate(se_100 = (sd / sqrt(100)))

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_02.png" # vector file too big just draw bitmap

xdata$estimates <- fdata$exprGene %>% group_by(subject) %>% summarize(mu_est = mean(expr))
xdata$mm <- fdata$mm$bsln %>% filter(gene == fdata$gene)

xdata$main <- xdata$estimates %>% 
  mutate(mu_true = xdata$mm$mean, 
         sd_true = xdata$mm$sd, 
         se_100_truth = xdata$mm$se_100)

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = subject)) + 
  geom_errorbar(aes(ymin = mu_true - se_100_truth, ymax = mu_true + se_100_truth), width = 0, size = 2, color = "grey70") +
  geom_line(aes(y = mu_true), group = 1, color = "red") +
  geom_point(aes(y = mu_true), color = "red", size = 2) +
  geom_line(aes(y = mu_est), group = 1, color = "blue") +
  geom_point(aes(y = mu_est), color = "blue", size = 2) +
  session$graphingUtils$tiltX(angle = 90) + 
  ylim(8, 60) + 
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank()) + 
  ylab("Expression (CPM)") + 
  xlab("Subject") + 
  ggtitle("xSubject variability in absence of ISV", "Gene: sim-GINM1") +
  theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 16))
  
xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 17, height = 15) 
# =================================================


# +++++++++++++++++++++++++=+++==========
# now, let's look at the version with SBJV imposed

fdata <- list(); gc() 

fdata$main <- read_rds(paste0(workspace$outputDir, "exprmats_agg_cell_sampling_100_cells.rds"))

fdata$exprmats <- fdata$main$exprmats
fdata$exprmats <- fdata$exprmats %>% lapply(function(currExprmats) { currExprmats[["excitatory"]] })
fdata$cells <- fdata$main$cells$excitatory

fdata$gene <- "GINM1"

# first let's look at the expression of this gene at baseline
fdata$exprGene <- pdata$simObjs$cellsSbjv$simdat$exprmat[fdata$gene, fdata$cells$sample]
fdata$exprGene <- fdata$exprGene %>% 
  session$dataWrangler$vectorToTibble() %>% 
  dplyr::select(sample = variable, expr = value) %>% 
  unique()

fdata$exprGene <- fdata$cells %>% left_join(fdata$exprGene, by = "sample")

# put together the adjusted distribution mean & se_100
fdata$mm <- pdata$simObjs$cellsSbjv$params$excitatory %>% lapply(function(currParams) {
  currParams %>% 
    workspace$utils$computeGamMM() %>% 
    session$dataWrangler$setRownameAsColumn("gene")
})

fdata$mm <- fdata$mm %>% session$collectionUtils$lapplyWithName(function(id, currParams) {
  currParams %>%
    mutate(sd = sqrt(variance)) %>% 
    mutate(se_100 = (sd / sqrt(100))) %>% 
    mutate(subject = id)
}) %>% session$dataWrangler$rbind()


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_03.png" # vector file too big just draw bitmap

xdata$estimates <- fdata$exprGene %>% group_by(subject) %>% summarize(mu_est = mean(expr))
xdata$mm <- fdata$mm %>% filter(gene == fdata$gene)

xdata$main <- xdata$estimates %>% left_join(xdata$mm, by = "subject")

xdata$main <- xdata$main %>% 
  dplyr::select(subject, mu_est, mu_true = mean, sd_truth = sd, se_100_true = se_100)

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = subject)) + 
  geom_errorbar(aes(ymin = mu_true - se_100_true, ymax = mu_true + se_100_true), width = 0, size = 2, color = "grey70") +
  geom_line(aes(y = mu_true), group = 1, color = "red") +
  geom_point(aes(y = mu_true), color = "red", size = 2) +
  geom_line(aes(y = mu_est), group = 1, color = "blue") +
  geom_point(aes(y = mu_est), color = "blue", size = 2) +
  session$graphingUtils$tiltX(angle = 90)  + 
  ylim(8, 60) + 
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank()) + 
  ylab("Expression (CPM)") + 
  xlab("Subject") + 
  ggtitle("xSubject variability in presence of ISV", "Gene: sim-GINM1") +
  theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 16))

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 17, height = 15) 
# =================================================


# +++++++++++++++++++++++++=+++==========
# plot the same thing but this time with 5000 cells!!
 
fdata <- list()

fdata$main <- read_rds(paste0(workspace$outputDir, "exprmats_agg_cell_sampling_1000_cells.rds"))

fdata$exprmats <- fdata$main$exprmats
fdata$exprmats <- fdata$exprmats %>% lapply(function(currExprmats) { currExprmats[["excitatory"]] })

fdata$rsqrs <- workspace$utils$computeRSqr(exprmatObs = fdata$exprmats$sbjv, exprmatCovs = fdata$exprmats[c("bsln", "islv")])

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_04.png" # vector file too big just draw bitmap

xdata$means <- fdata$rsqrs %>% summarize(rsqr_bsln = mean(rsqr_bsln), rsqr_islv = mean(rsqr_islv))

# > xdata$means 
# # A tibble: 1 × 2
# rsqr_bsln rsqr_islv
# <dbl>     <dbl>
#   1     0.143     0.877

xdata$plot <- fdata$rsqrs %>% session$graphingUtils$ggplot(aes(x = rsqr_islv , y = rsqr_bsln)) + 
  geom_point(shape = 1) + 
  geom_point(data = fdata$rsqrs %>% filter(gene_id == "GINM1"), color = "red", size = 3) +
  geom_vline(data = xdata$means, aes(xintercept = rsqr_islv), linetype = "dashed", color = "red") +
  geom_hline(data = xdata$means, aes(yintercept = rsqr_bsln), linetype = "dashed", color = "red") +
  xlim(0, 1) + 
  ylim(0, 1) + 
  xlab(bquote("PVE by ISV"~(R^2))) +
  ylab(bquote("PVE by CSIV"~(R^2))) + 
  ggtitle("CSIV vs. ISV for 1000 cells per sample") +
  theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 16))

# xdata$plot <- xdata$plot %>% ggMarginal(type = "histogram", size = 3)

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 15, height = 15) 
# =================================================


# +++++++++++++++++++++++++=+++==========
# now lets explore the dynamics of these genes expressions 
# first at baseline (only sampling)

fdata <- list(); gc() 

fdata$main <- read_rds(paste0(workspace$outputDir, "exprmats_agg_cell_sampling_1000_cells.rds"))

fdata$exprmats <- fdata$main$exprmats
fdata$exprmats <- fdata$exprmats %>% lapply(function(currExprmats) { currExprmats[["excitatory"]] })
fdata$cells <- fdata$main$cells$excitatory

fdata$gene <- "GINM1" 

# first let's look at the expression of this gene at baseline
fdata$exprGene <- pdata$simObjs$cellsBsln$simdat$exprmat[fdata$gene, fdata$cells$sample]
fdata$exprGene <- fdata$exprGene %>% 
  session$dataWrangler$vectorToTibble() %>% 
  dplyr::select(sample = variable, expr = value) %>% 
  unique()

fdata$exprGene <- fdata$cells %>% left_join(fdata$exprGene, by = "sample")

# put together the bsln stuff
fdata$mm$bsln <- pdata$simObjs$cellsBsln$params$md$excitatory %>% 
  workspace$utils$computeGamMM() %>% 
  session$dataWrangler$setRownameAsColumn("gene")

fdata$mm$bsln <- fdata$mm$bsln %>%
  mutate(sd = sqrt(variance)) %>% 
  mutate(se_100 = (sd / sqrt(1000)))


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_05.png" # vector file too big just draw bitmap

xdata$estimates <- fdata$exprGene %>% group_by(subject) %>% summarize(mu_est = mean(expr))
xdata$mm <- fdata$mm$bsln %>% filter(gene == fdata$gene)

xdata$main <- xdata$estimates %>% 
  mutate(mu_true = xdata$mm$mean, 
         sd_true = xdata$mm$sd, 
         se_100_truth = xdata$mm$se_100)

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = subject)) + 
  geom_errorbar(aes(ymin = mu_true - se_100_truth, ymax = mu_true + se_100_truth), width = 0, size = 2, color = "grey70") +
  geom_line(aes(y = mu_true), group = 1, color = "red") +
  geom_point(aes(y = mu_true), color = "red", size = 2) +
  geom_line(aes(y = mu_est), group = 1, color = "blue") +
  geom_point(aes(y = mu_est), color = "blue", size = 2) +
  session$graphingUtils$tiltX(angle = 90) + 
  ylim(8, 60) + 
  theme(axis.ticks.x = element_blank(), 
      axis.text.x = element_blank()) + 
  ylab("Expression (CPM)") + 
  xlab("Subject") + 
  ggtitle("xSubject variability in absence of ISV", "Gene: sim-GINM1") +
  theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 16))

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 17, height = 15) 
# =================================================


# +++++++++++++++++++++++++=+++==========
# now, let's look at the version with SBJV imposed

fdata <- list(); gc() 

fdata$main <- read_rds(paste0(workspace$outputDir, "exprmats_agg_cell_sampling_1000_cells.rds"))

fdata$exprmats <- fdata$main$exprmats
fdata$exprmats <- fdata$exprmats %>% lapply(function(currExprmats) { currExprmats[["excitatory"]] })
fdata$cells <- fdata$main$cells$excitatory

fdata$gene <- "GINM1"

# first let's look at the expression of this gene at baseline
fdata$exprGene <- pdata$simObjs$cellsSbjv$simdat$exprmat[fdata$gene, fdata$cells$sample]
fdata$exprGene <- fdata$exprGene %>% 
  session$dataWrangler$vectorToTibble() %>% 
  dplyr::select(sample = variable, expr = value) %>% 
  unique()

fdata$exprGene <- fdata$cells %>% left_join(fdata$exprGene, by = "sample")

# put together the adjusted distribution mean & se_100
fdata$mm <- pdata$simObjs$cellsSbjv$params$excitatory %>% lapply(function(currParams) {
  currParams %>% 
    workspace$utils$computeGamMM() %>% 
    session$dataWrangler$setRownameAsColumn("gene")
})

fdata$mm <- fdata$mm %>% session$collectionUtils$lapplyWithName(function(id, currParams) {
  currParams %>%
    mutate(sd = sqrt(variance)) %>% 
    mutate(se_100 = (sd / sqrt(1000))) %>% 
    mutate(subject = id)
}) %>% session$dataWrangler$rbind()


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_06.png" # vector file too big just draw bitmap

xdata$estimates <- fdata$exprGene %>% group_by(subject) %>% summarize(mu_est = mean(expr))
xdata$mm <- fdata$mm %>% filter(gene == fdata$gene)

xdata$main <- xdata$estimates %>% left_join(xdata$mm, by = "subject")

xdata$main <- xdata$main %>% 
  dplyr::select(subject, mu_est, mu_true = mean, sd_truth = sd, se_100_true = se_100)

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = subject)) + 
  geom_errorbar(aes(ymin = mu_true - se_100_true, ymax = mu_true + se_100_true), width = 0, size = 2, color = "grey70") +
  geom_line(aes(y = mu_true), group = 1, color = "red") +
  geom_point(aes(y = mu_true), color = "red", size = 2) +
  geom_line(aes(y = mu_est), group = 1, color = "blue") +
  geom_point(aes(y = mu_est), color = "blue", size = 2) +
  session$graphingUtils$tiltX(angle = 90) + 
  ylim(8, 60) + 
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank()) + 
  ylab("Expression (CPM)") + 
  xlab("Subject") + 
  ggtitle("xSubject variability in absence of ISV", "Gene: sim-GINM1") +
  theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 16))

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 17, height = 15) 
# =================================================

# +++++++++++++++++++++++++++++++++++++++============
# compute for all the cell types
 
fdata <- list(); gc()

fdata$simObj <- pdata$simObjs$cellsBsln
fdata$samples <- fdata$simObj$simdat$samples

fdata$celltypes <- c("excitatory", "inhibitory", "oligodendrocyte", "astrocyte", "microglia", "opc") %>% session$dataWrangler$attachNames()
fdata$nCels <- c(9750, 3250, 8775, 5850, 1950, 2925) %>% session$dataWrangler$attachNames(useNames = fdata$celltypes)

fdata$sbjs <- pdata$simObjs$cellsBsln$simdat$samples$subject %>% unique() %>% sort() %>% session$dataWrangler$attachNames()

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
  workspace$utils$aggregSbj(pdata$simObjs$cellsBsln$simdat$exprmat, samples = currCells)
})

fdata$exprmats$sbjv <- fdata$cells %>% session$collectionUtils$lapply(function(currCells) {
  workspace$utils$aggregSbj(pdata$simObjs$cellsSbjv$simdat$exprmat, samples = currCells)
})

fdata$exprmats$islv <- fdata$celltypes %>% session$collectionUtils$lapply(function(celltype) {
  samples <- pdata$simObjs$subjects$simdat$samples %>% filter(cell_type == celltype)
  exprmat <- pdata$simObjs$subjects$simdat$exprmat[, samples$sample]
  colnames(exprmat) <- samples$subject
  return(exprmat)
})

# COMMIT =====
fdata$output <- list()
fdata$output$cells <- fdata$cells 
fdata$output$exprmats <- fdata$exprmats
fdata$output %>% saveRDS(paste0(workspace$outputDir, "exprmats_agg_cell_sampling.rds"))
# ============

# +++++++++++++++++++++++++++++++++++++++============
# OUTPUT FIGURE

fdata <- list(); gc()

xdata$figName <- "figure_07.eps" # vector file too big just draw bitmap

fdata$main <- read_rds(paste0(workspace$outputDir, "exprmats_agg_cell_sampling.rds"))

fdata$sbjs <- fdata$main$cells$excitatory$subject %>% unique() %>% sort()
fdata$cells <- fdata$main$cells %>% 
  session$collectionUtils$lapplyWithName(function(celltype, cells) { tibble(cell_type = celltype, n_cells = nrow(cells) / length(fdata$sbjs) ) }) %>% 
  session$dataWrangler$rbind()

fdata$cells <- fdata$cells %>% mutate(cell_type = cell_type %>% workspace$utils$fmtCelltypes())

xdata$plot <- fdata$cells %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = n_cells)) + 
  geom_bar(aes(fill = cell_type), stat = "identity") + 
  geom_text(aes(label = n_cells), size = 5, vjust = -1) + 
  session$graphingUtils$tiltX(angle = 90) + 
  theme(legend.position = "none") + 
  ylim(0, 10500) + 
  ggtitle("Estimated number of cells\nin 1 mg of cortical tissue") +
  theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 16)) +
  xlab("Cell type") + 
  ylab("Number of cells") + 
  scale_fill_brewer(palette = "Set1")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, units = "cm", width = 13, height = 15) 
# =================================================


fdata$exprmats <- fdata$main$exprmats

fdata$celltypes <- fdata$exprmats$bsln %>% names() %>% sort() %>% session$dataWrangler$attachNames()

fdata$rsqrs <- fdata$celltypes %>% session$collectionUtils$lapply(function(celltype) {
  
  exprmats <- fdata$exprmats %>% lapply(function(currExprmats) { currExprmats[[celltype]] })
  
  rsqrs <- workspace$utils$computeRSqr(exprmatObs = exprmats$sbjv, exprmatCovs = exprmats[c("bsln", "islv")])
  
  vars <- exprmats %>% 
    sapply(function(exprmat) { exprmat %>% apply(1, var) })
  colnames(vars) <- paste0("var_", colnames(vars))
  vars <- vars %>% session$dataWrangler$setRownameAsColumn("gene_id")
  
  means <- exprmats %>% 
    sapply(function(exprmat) { exprmat %>% apply(1, mean) })
  colnames(means) <- paste0("mean_", colnames(means))
  means <- means %>% session$dataWrangler$setRownameAsColumn("gene_id")
  
  rsqrs %>% left_join(vars, by = "gene_id") %>% left_join(means, by = "gene_id") %>%  mutate(cell_type = celltype)
  
}) %>% session$dataWrangler$rbind()

fdata$rsqrs %>% group_by(cell_type) %>% summarize(cor(rsqr_bsln, rsqr_islv))

fdata$nCels <- c("excitatory", "inhibitory", "oligodendrocyte", "astrocyte", "microglia", "opc") %>% session$dataWrangler$attachNames()
fdata$nCels <- c(9750, 3250, 8775, 5850, 1950, 2925) %>% session$dataWrangler$attachNames(useNames = fdata$nCels)

fdata$rsqrs$n_cel <- fdata$nCels[fdata$rsqrs$cell_type]

# =================================================
# OUTPUT S_FILE
xdata <- list()

xdata$fileName <- "sfile_04_csiv_vs_isv.tsv" 

xdata$colnames <- c(gene_id = "gene_d", 
                    rsqr_bsln = "rsqr_csiv", 
                    rsqr_islv = "rsqr_isv", 
                    var_bsln = "var_csiv", 
                    var_sbjv = "var_combined", 
                    var_islv = "var_isv", 
                    mean_bsln = "mean_csiv", 
                    mean_sbjv = "mean_combined", 
                    mean_islv = "mean_isv", 
                    cell_type = "cell_type", 
                    n_cel = "n_cel")

xdata$dat <- fdata$rsqrs
colnames(xdata$dat) <- xdata$colnames[colnames(xdata$dat)]

write.table(xdata$dat, file = paste0(pdata$sfilesDir, xdata$fileName), 
            sep = "\t", row.names = FALSE) 
# =================================================



# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_08.eps"

xdata$means <- fdata$rsqrs %>% group_by(cell_type) %>% 
  summarize(rsqr_bsln = mean(rsqr_bsln))

xdata$plot <- fdata$rsqrs %>% 
  session$graphingUtils$ggplot(aes(x = rsqr_bsln)) + 
  geom_density(aes(color = cell_type)) + 
  xlim(0, 1) + 
  ylab("Density") + 
  xlab(bquote("PVE by CSIV"~(R^2))) +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "none") +
  ggtitle("CSIV importance by cell type")
  
xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 14, height = 15) 


# 
# xdata$means %>% 
#   mutate(cell_type = workspace$utils$fmtCelltypes(cell_type)) %>% 
#   session$graphingUtils$ggplot(aes(x = cell_type, y = rsqr_bsln)) + 
#   geom_bar(stat = "identity", aes(fill = cell_type)) + 
#   session$graphingUtils$tiltX(angle = 90) + 
#   ylim(0, 1) +
#   xlab("Cell type") + 
#   ylab("") + 
#   scale_fill_brewer(palette = "Set1") + 
#   theme(legend.position = "none")


# =================================================

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_09.eps"

xdata$plot <- fdata$rsqrs %>% 
  session$graphingUtils$ggplot(aes(x = rsqr_islv)) + 
  geom_density(aes(color = cell_type)) + 
  xlim(0, 1) + 
  ylab("Density") + 
  xlab(bquote("PVE by ISV"~(R^2))) +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "none") +
  ggtitle("ISV importance by cell type") 

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 14, height = 15) 
# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_10.eps" # vector file too big just draw bitmap

xdata$plot <- fdata$rsqrs %>% 
  mutate(rv = var_bsln / mean_bsln) %>% 
  session$graphingUtils$ggplot(aes(x = (rv))) + 
  geom_density(aes(color = cell_type)) +
  xlim(0, 1) + 
  ylab("Density") + 
  xlab("Relative variance of CSIV (v/m)") + 
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Set1") +
  ggtitle("CSIV magnitude by cell type")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 14, height = 15) 
# =================================================








# this is kind of interesting - move it somewhere else to a sandbox file or whatever
# TODO - sandbox, get rid of later

fdata <- list(); gc()

fdata$sbjs <- paste0("s_", formatC(1:1000, digits = 5, flag = "0"))

# simulated the 300 correlated genes
fdata$genes$coex <- paste0("gc_", formatC(1:300, digits = 5, flag = "0"))
fdata$sigma <- matrix(0.2, nrow = 300, ncol = 300)
rownames(fdata$sigma) <- fdata$genes$coex
colnames(fdata$sigma) <- fdata$genes$coex
diag(fdata$sigma) <- 1

set.seed(0); fdata$exprmat$coex <- rmvnorm(n = 1000, sigma = fdata$sigma, checkSymmetry = FALSE) %>% t()
rownames(fdata$exprmat$coex) <- fdata$genes$coex
colnames(fdata$exprmat$coex) <- fdata$sbjs

# simulated the 2700 null correlation genes
fdata$genes$coex <- paste0("gn_", formatC(1:2700, digits = 5, flag = "0"))
fdata$exprmat$null <- matrix(rnorm(2700 * 1000), nrow = 2700) # have to impose a correlation here
rownames(fdata$exprmat$null) <- fdata$genes$coex 
colnames(fdata$exprmat$null) <- fdata$sbjs

# put them together
fdata$exprmat$all <- fdata$exprmat %>% session$dataWrangler$rbind()

fdata$exprmat$all %>% dim()

fdata$coexmat$n1000 <- fdata$exprmat$all %>% workspace$utils$computeCoexmat()
set.seed(0); fdata$coexmat$n500 <- fdata$exprmat$all[, sample(colnames(fdata$exprmat$all), 500)] %>% workspace$utils$computeCoexmat()
set.seed(0); fdata$coexmat$n100 <- fdata$exprmat$all[, sample(colnames(fdata$exprmat$all), 100)] %>% workspace$utils$computeCoexmat()
set.seed(0); fdata$coexmat$n50 <- fdata$exprmat$all[, sample(colnames(fdata$exprmat$all), 50)] %>% workspace$utils$computeCoexmat()
set.seed(0); fdata$coexmat$n40 <- fdata$exprmat$all[, sample(colnames(fdata$exprmat$all), 40)] %>% workspace$utils$computeCoexmat()
set.seed(0); fdata$coexmat$n30 <- fdata$exprmat$all[, sample(colnames(fdata$exprmat$all), 30)] %>% workspace$utils$computeCoexmat()
set.seed(0); fdata$coexmat$n20 <- fdata$exprmat$all[, sample(colnames(fdata$exprmat$all), 20)] %>% workspace$utils$computeCoexmat()
set.seed(0); fdata$coexmat$n10 <- fdata$exprmat$all[, sample(colnames(fdata$exprmat$all), 10)] %>% workspace$utils$computeCoexmat()
set.seed(0); fdata$coexmat$n5 <- fdata$exprmat$all[, sample(colnames(fdata$exprmat$all), 5)] %>% workspace$utils$computeCoexmat()

fdata$coextbl <- fdata$coexmat %>% session$collectionUtils$lapply(function(currCoexmat) {
  tibble(pair = currCoexmat %>% workspace$utils$getPairIds(), 
         coex = currCoexmat %>% workspace$utils$vectorize())
})

fdata$coextop <- fdata$coextbl %>% session$collectionUtils$lapply(function(currTbl) {
  currTbl %>% mutate(rank = rank(coex)) %>% mutate(perc = rank / max(rank)) %>% filter(perc >= 0.99)
})

fdata$main <- fdata$coextop %>% session$collectionUtils$lapplyWithName(function(id, edges) {
  tibble(n_cel = id, n_overlap = length(intersect(edges$pair, fdata$coextop$n1000$pair)), n_total = nrow(edges))
}) %>% session$dataWrangler$rbind() %>% mutate(frac = n_overlap / n_total)

fdata$main






