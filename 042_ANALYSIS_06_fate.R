pdata <- list()

pdata$figuresDir <- paste0(workspace$workspaceDir, "042_ANALYSIS_06_fate_FIGURES/")

pdata$sfilesDir <- paste0(workspace$workspaceDir, "042_S_FILES/")

# just reuse the simulations from the previous part

pdata$cteprf <- read_rds(paste0(workspace$outputDir, "seed_cteprf.rds"))

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


pdata$utils$xbulkExprmats <- function(sbjvExprmats, cv) {
  
  celltypes <- c("excitatory", "inhibitory", "oligodendrocyte", "astrocyte", "microglia", "opc") %>% session$dataWrangler$attachNames()
  nCels <- c(9750, 3250, 8775, 5850, 1950, 2925) %>% session$dataWrangler$attachNames(useNames = celltypes)
  fracCells <- nCels / sum(nCels)
  sdCells <- nCels * cv
  
  nSbjs <- ncol(sbjvExprmats[[1]])
  sbjs <- colnames(sbjvExprmats[[1]])
  
  ccv <- celltypes %>% sapply(function(celltype) {
    rand <- rnorm(nSbjs, mean = nCels[celltype], sd = sdCells[celltype]) * (fracCells[celltype] / nCels[celltype])
    abs(rand) # take the absolute fraction to prevent negative numbers; shouldn't have a big impact on SD
  })
  rownames(ccv) <- sbjs
  
  exprmats <- list()
  
  exprmats$ccv <- Reduce("+", celltypes %>% lapply(function(celltype) {
      ccvCelltype <- ccv[, celltype]
      exprmat <- sbjvExprmats[[celltype]] 
      exprmat %>% apply(1, function(vals) { vals * ccvCelltype }) %>% t() 
  }))
  
  exprmats$null <- Reduce("+", celltypes %>% lapply(function(celltype) {
    sbjvExprmats[[celltype]] * fracCells[celltype]
  }))
  
  return(exprmats)
}


pdata$utils$extracCelLvlExprmat <- function(exprmat, samples, celltype, sbj) {
  exprmat[, (samples %>% filter(cell_type == celltype, subject == sbj))$sample ]
}


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++===========
# segregate 3000  genes into 5 different scenarios, each has 600 genes

fdata <- list(); gc()

fdata$celltypes <- pdata$cteprf$mkrs$mkr_ct %>% unique() %>% sort() %>% session$dataWrangler$attachNames()

set.seed(0); fdata$genes <- rownames(pdata$cteprf$cteprf) %>% sample(3000, replace = FALSE) # shuffle it

# separate into 5 groups of 600 genes
fdata$genes <- tibble(gene_id = fdata$genes, scenario_i = ceiling((seq_along(fdata$genes) / 600)))
fdata$genes <- fdata$genes %>% mutate(scenario = paste0("s_", formatC(scenario_i, digits = 1, flag = "0"))) 

fdata$genes %>% group_by(scenario) %>% summarize(n_gene = n()) # double check that genes are binned correctly

fdata$genesFlat <- fdata$genes$scenario %>% unique() %>% sort() %>% session$dataWrangler$attachNames() %>% 
  lapply(function(currScenario) { (fdata$genes %>% filter(scenario == currScenario))$gene_id })

# put together the co-expression programs

fdata$this <- read_rds(paste0(workspace$outputDir, "this_cell_sampling.rds"))

# define coexpression programs
fdata$coexPrgms <- list() 

# IDEA - run a separate simulation for each scenario
# for each scenario, define coex program, run simualtion, save simulation to disk

# in scenario 1, only xcell co-expression is present
xdata <- list(); gc()
xdata$genes <- fdata$genesFlat$s_01
xdata$coexPrgms$s_01_00 <- list(celltypes = "excitatory", genes = xdata$genes, value = 0.8)
xdata$this <- fdata$this %>% workspace$api$setRefDatCoexPrgms(xdata$coexPrgms, level = "cel") 
xdata$this$mdParams$cel <- xdata$this$mdParams$cel %>% lapply(function(mdParam) { mdParam[xdata$genes, ] })
xdata$this$mdParams$sbj <- xdata$this$mdParams$sbj %>% lapply(function(mdParam) { mdParam[xdata$genes, ] })
xdata$simObjs <- xdata$this %>% pdata$utils$simulate(5000, 30)
# COMMIT ==========================
xdata$simObjs %>% saveRDS(paste0(workspace$outputDir, "simObjs_fate_scenario_1.rds"))
# =================================  

# in scenario 2, only xcell co-expression are propagated to the xsbj level
xdata <- list(); gc()
xdata$genes <- fdata$genesFlat$s_02
xdata$coexPrgms$s_02_00 <- list(celltypes = "excitatory", genes = xdata$genes, value = 0.8)
xdata$this <- fdata$this %>% workspace$api$setRefDatCoexPrgms(xdata$coexPrgms, level = "cel") 
xdata$this <- xdata$this %>% workspace$api$setRefDatCoexPrgms(xdata$coexPrgms, level = "sbj") 
xdata$this$mdParams$cel <- xdata$this$mdParams$cel %>% lapply(function(mdParam) { mdParam[xdata$genes, ] })
xdata$this$mdParams$sbj <- xdata$this$mdParams$sbj %>% lapply(function(mdParam) { mdParam[xdata$genes, ] })
xdata$simObjs <- xdata$this %>% pdata$utils$simulate(5000, 30)
# COMMIT ==========================
xdata$simObjs %>% saveRDS(paste0(workspace$outputDir, "simObjs_fate_scenario_2.rds"))
# =================================  

# in scenario 3, only xcell co-expression are propagated to the xsbj level, and independently co-regulated witihn each cell type
# at both levels
xdata <- list(); gc()
xdata$genes <- fdata$genesFlat$s_03
xdata$coexPrgms$s_03_exc <- list(celltypes = "excitatory", genes = xdata$genes, value = 0.8)
xdata$coexPrgms$s_03_ihn <- list(celltypes = "inhibitory", genes = xdata$genes, value = 0.8)
xdata$coexPrgms$s_03_ast <- list(celltypes = "astrocyte", genes = xdata$genes, value = 0.8)
xdata$coexPrgms$s_03_mic <- list(celltypes = "microglia", genes = xdata$genes, value = 0.8)
xdata$coexPrgms$s_03_olig <- list(celltypes = "oligodendrocyte", genes = xdata$genes, value = 0.8)
xdata$coexPrgms$s_03_opc <- list(celltypes = "opc", genes = xdata$genes, value = 0.8)
xdata$this <- fdata$this %>% workspace$api$setRefDatCoexPrgms(xdata$coexPrgms, level = "cel") 
xdata$this <- xdata$this %>% workspace$api$setRefDatCoexPrgms(xdata$coexPrgms, level = "sbj") 
xdata$this$mdParams$cel <- xdata$this$mdParams$cel %>% lapply(function(mdParam) { mdParam[xdata$genes, ] })
xdata$this$mdParams$sbj <- xdata$this$mdParams$sbj %>% lapply(function(mdParam) { mdParam[xdata$genes, ] })
xdata$simObjs <- xdata$this %>% pdata$utils$simulate(5000, 30)
# COMMIT ==========================
xdata$simObjs %>% saveRDS(paste0(workspace$outputDir, "simObjs_fate_scenario_3.rds"))
# =================================  

# in scenario 4, co-expressed and synchronized  
xdata <- list(); gc()
xdata$genes <- fdata$genesFlat$s_04
xdata$coexPrgms$s_04 <- list(celltypes = c("excitatory", "inhibitory", "astrocyte", "microglia", "oligodendrocyte", "opc"), 
                             genes = xdata$genes, value = 0.8)
xdata$this <- fdata$this %>% workspace$api$setRefDatCoexPrgms(xdata$coexPrgms, level = "cel") 
xdata$this <- xdata$this %>% workspace$api$setRefDatCoexPrgms(xdata$coexPrgms, level = "sbj") 
xdata$this$mdParams$cel <- xdata$this$mdParams$cel %>% lapply(function(mdParam) { mdParam[xdata$genes, ] })
xdata$this$mdParams$sbj <- xdata$this$mdParams$sbj %>% lapply(function(mdParam) { mdParam[xdata$genes, ] })
xdata$simObjs <- xdata$this %>% pdata$utils$simulate(5000, 30)
# COMMIT ==========================
xdata$simObjs %>% saveRDS(paste0(workspace$outputDir, "simObjs_fate_scenario_4.rds"))
# =================================    
   
# in scenario 5, no co-expression   
xdata <- list(); gc() 
xdata$genes <- fdata$genesFlat$s_04 
xdata$coexPrgms$s_04 <- list()
xdata$this <- fdata$this %>% workspace$api$setRefDatCoexPrgms(xdata$coexPrgms, level = "cel") 
xdata$this <- xdata$this %>% workspace$api$setRefDatCoexPrgms(xdata$coexPrgms, level = "sbj") 
xdata$this$mdParams$cel <- xdata$this$mdParams$cel %>% lapply(function(mdParam) { mdParam[xdata$genes, ] })
xdata$this$mdParams$sbj <- xdata$this$mdParams$sbj %>% lapply(function(mdParam) { mdParam[xdata$genes, ] })
xdata$simObjs <- xdata$this %>% pdata$utils$simulate(5000, 30)
# COMMIT ==========================
xdata$simObjs %>% saveRDS(paste0(workspace$outputDir, "simObjs_fate_scenario_5.rds"))
# =================================  

# +++++++++++++++++++++++=====================
# load in each simulation scenario
# 1. check their co-expression patterns
# 2. construct expression matrix at the different levels
# 3. compute coexmat & coextbl for each
# 4. save everything

# Scenario 1 - co-expression in excitatory neurons at the xcell level only

fdata <- list(); gc()

fdata$simObjs <- read_rds(paste0(workspace$outputDir, "simObjs_fate_scenario_1.rds"))
fdata$exprmats$xsbj <- pdata$utils$xSbjExprmats(fdata$simObjs)
fdata$exprmats$xbulk <- pdata$utils$xbulkExprmats(fdata$exprmats$xsbj$exprmats$sbjv, cv = 0.3)

# compute co-expressions
fdata$coexmats <- fdata$simObjs$cellsBsln$simdat$samples

# exprmat at different levesl for a single cell types

fdata$celltype <- "excitatory"

fdata$exprmats2$xcel <- fdata$simObjs$cellsBsln$simdat$exprmat %>% 
  pdata$utils$extracCelLvlExprmat(fdata$simObjs$cellsBsln$simdat$samples, fdata$celltype, "sbj_001")
fdata$exprmats2$bsln <- fdata$exprmats$xsbj$exprmats$bsln[[fdata$celltype]]
fdata$exprmats2$xsbj <- fdata$exprmats$xsbj$exprmats$sbjv[[fdata$celltype]]
fdata$exprmats2$diluted <- fdata$exprmats$xbulk$null
fdata$exprmats2$ccv <- fdata$exprmats$xbulk$ccv

fdata$coexmats <- fdata$exprmats2 %>% session$collectionUtils$lapply(workspace$utils$computeCoexmat)

# fdata$coexmats <- fdata$exprmats2 %>% session$collectionUtils$lapply(workspace$utils$computeCoexmat)

fdata$coextbls <- fdata$coexmats %>% session$collectionUtils$lapplyWithName(function(lvl, coexmat) {
  tibble(pair_id = workspace$utils$getPairIds(coexmat), coex = workspace$utils$vectorize(coexmat)) %>% 
    mutate(level = lvl)
}) %>% session$dataWrangler$rbind()

fdata$coextbls %>% session$graphingUtils$ggplot(aes(x = coex)) + geom_density(aes(color = level))

fdata$coextbls %>% group_by(level) %>% summarize(coex = mean(coex))
fdata$coextbls %>% group_by(level) %>% summarize(coex_max = max(coex))

# +++++++++=========== write to disk
fdata$coextbls %>% saveRDS(paste0(workspace$outputDir, "coextbls_fate_scenario_1.rds"))
# ===================================

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_01.png" 

xdata$means <- fdata$coextbls %>% 
  filter(scenario == "scenario_1") %>% 
  group_by(level) %>% summarize(coex = mean(coex)) 

set.seed(0); xdata$pairsSample <- fdata$coextbls$pair_id %>% unique() %>% sample(1000)

xdata$plot <- fdata$coextbls %>% 
  filter(pair_id %in% xdata$pairsSample) %>% 
  session$graphingUtils$ggplot(aes(x = level, y = coex)) + 
  geom_line(aes(group = pair_id), size = 0.05) + 
  geom_line(data = xdata$means, aes(group = 1), color = "red", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  xlim(c("xcel", "bsln", "xsbj", "diluted", "ccv")) + 
  session$graphingUtils$tiltX(angle = 90) + 
  ylim(-1, 1) +
  ylab("Co-expression (Pearson's r)") +
  ggtitle("Scenario 2: genes co-regulated at the xCell level only") + 
  theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 14), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 20, height = 15) 
# =================================================


# +++++++++++++++++++++++====================
# Scenario 2

fdata <- list(); gc()

fdata$simObjs <- read_rds(paste0(workspace$outputDir, "simObjs_fate_scenario_2.rds"))
fdata$exprmats$xsbj <- pdata$utils$xSbjExprmats(fdata$simObjs)
fdata$exprmats$xbulk <- pdata$utils$xbulkExprmats(fdata$exprmats$xsbj$exprmats$sbjv, cv = 0.3)

# compute co-expressions
fdata$coexmats <- fdata$simObjs$cellsBsln$simdat$samples

# exprmat at different levesl for a single cell types

fdata$celltype <- "excitatory"

fdata$exprmats2$xcel <- fdata$simObjs$cellsBsln$simdat$exprmat %>% 
  pdata$utils$extracCelLvlExprmat(fdata$simObjs$cellsBsln$simdat$samples, fdata$celltype, "sbj_001")
fdata$exprmats2$bsln <- fdata$exprmats$xsbj$exprmats$bsln[[fdata$celltype]]
fdata$exprmats2$xsbj <- fdata$exprmats$xsbj$exprmats$sbjv[[fdata$celltype]]
fdata$exprmats2$diluted <- fdata$exprmats$xbulk$null
fdata$exprmats2$ccv <- fdata$exprmats$xbulk$ccv

fdata$coexmats <- fdata$exprmats2 %>% session$collectionUtils$lapply(workspace$utils$computeCoexmat)

fdata$coexmats <- fdata$exprmats2 %>% session$collectionUtils$lapply(workspace$utils$computeCoexmat)

fdata$coextbls <- fdata$coexmats %>% session$collectionUtils$lapplyWithName(function(lvl, coexmat) {
  tibble(pair_id = workspace$utils$getPairIds(coexmat), coex = workspace$utils$vectorize(coexmat)) %>% 
    mutate(level = lvl)
}) %>% session$dataWrangler$rbind()

fdata$coextbls %>% session$graphingUtils$ggplot(aes(x = coex)) + geom_density(aes(color = level))

fdata$coextbls %>% group_by(level) %>% summarize(coex = mean(coex))
fdata$coextbls %>% group_by(level) %>% summarize(coex_max = max(coex))

# +++++++++=========== write to disk
fdata$coextbls %>% saveRDS(paste0(workspace$outputDir, "coextbls_fate_scenario_2.rds"))
# ===================================

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_02.png" 

xdata$means <- fdata$coextbls %>% group_by(level) %>% summarize(coex = mean(coex))

set.seed(0); xdata$pairsSample <- fdata$coextbls$pair_id %>% unique() %>% sample(1000)

xdata$plot <- fdata$coextbls %>% 
  filter(pair_id %in% xdata$pairsSample) %>% 
  session$graphingUtils$ggplot(aes(x = level, y = coex)) + 
  geom_line(aes(group = pair_id), size = 0.05) + 
  geom_line(data = xdata$means, aes(group = 1), color = "red", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  xlim(c("xcel", "bsln", "xsbj", "diluted", "ccv")) + 
  session$graphingUtils$tiltX(angle = 90) + 
  ylim(-1, 1) +
  ylab("Co-expression (Pearson's r)") +
  ggtitle("Scenario 3: genes co-regulated at both the xCell and\nxSubject levels") + 
  theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 14), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 20, height = 15) 
# =================================================

# +++++++++++++++++++++++====================
# Scenario 3

fdata <- list(); gc()

fdata$simObjs <- read_rds(paste0(workspace$outputDir, "simObjs_fate_scenario_3.rds"))
fdata$exprmats$xsbj <- pdata$utils$xSbjExprmats(fdata$simObjs)
fdata$exprmats$xbulk <- pdata$utils$xbulkExprmats(fdata$exprmats$xsbj$exprmats$sbjv, cv = 0.3)

# compute co-expressions
fdata$coexmats <- fdata$simObjs$cellsBsln$simdat$samples

# exprmat at different levesl for a single cell types

fdata$celltype <- "excitatory"

fdata$exprmats2$xcel <- fdata$simObjs$cellsBsln$simdat$exprmat %>% 
  pdata$utils$extracCelLvlExprmat(fdata$simObjs$cellsBsln$simdat$samples, fdata$celltype, "sbj_001")
fdata$exprmats2$bsln <- fdata$exprmats$xsbj$exprmats$bsln[[fdata$celltype]]
fdata$exprmats2$xsbj <- fdata$exprmats$xsbj$exprmats$sbjv[[fdata$celltype]]
fdata$exprmats2$diluted <- fdata$exprmats$xbulk$null
fdata$exprmats2$ccv <- fdata$exprmats$xbulk$ccv

fdata$coexmats <- fdata$exprmats2 %>% session$collectionUtils$lapply(workspace$utils$computeCoexmat)

fdata$coexmats <- fdata$exprmats2 %>% session$collectionUtils$lapply(workspace$utils$computeCoexmat)

fdata$coextbls <- fdata$coexmats %>% session$collectionUtils$lapplyWithName(function(lvl, coexmat) {
  tibble(pair_id = workspace$utils$getPairIds(coexmat), coex = workspace$utils$vectorize(coexmat)) %>% 
    mutate(level = lvl)
}) %>% session$dataWrangler$rbind()

fdata$coextbls %>% session$graphingUtils$ggplot(aes(x = coex)) + geom_density(aes(color = level))

fdata$coextbls %>% group_by(level) %>% summarize(coex = mean(coex))
fdata$coextbls %>% group_by(level) %>% summarize(coex_max = max(coex))

# +++++++++=========== write to disk
fdata$coextbls %>% saveRDS(paste0(workspace$outputDir, "coextbls_fate_scenario_3.rds"))
# ===================================

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_03.png" 

xdata$means <- fdata$coextbls %>% group_by(level) %>% summarize(coex = mean(coex))

set.seed(1); xdata$pairsSample <- fdata$coextbls$pair_id %>% unique() %>% sample(1000)

xdata$plot <- fdata$coextbls %>% 
  filter(pair_id %in% xdata$pairsSample) %>% 
  session$graphingUtils$ggplot(aes(x = level, y = coex)) + 
  geom_line(aes(group = pair_id), size = 0.05) + 
  geom_line(data = xdata$means, aes(group = 1), color = "red", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  xlim(c("xcel", "bsln", "xsbj", "diluted", "ccv")) + 
  session$graphingUtils$tiltX(angle = 90) + 
  ylim(-1, 1) +
  ylab("Co-expression (Pearson's r)") +
  ggtitle("Scenario 4: genes co-regulated in all cell types") + 
  theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 14), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 20, height = 15) 
# =================================================

# +++++++++++++++++++++++====================
# Scenario 4

fdata <- list(); gc()

fdata$simObjs <- read_rds(paste0(workspace$outputDir, "simObjs_fate_scenario_4.rds"))
fdata$exprmats$xsbj <- pdata$utils$xSbjExprmats(fdata$simObjs)
fdata$exprmats$xbulk <- pdata$utils$xbulkExprmats(fdata$exprmats$xsbj$exprmats$sbjv, cv = 0.3)

# compute co-expressions
fdata$coexmats <- fdata$simObjs$cellsBsln$simdat$samples

# exprmat at different levesl for a single cell types

fdata$celltype <- "excitatory"

fdata$exprmats2$xcel <- fdata$simObjs$cellsBsln$simdat$exprmat %>% 
  pdata$utils$extracCelLvlExprmat(fdata$simObjs$cellsBsln$simdat$samples, fdata$celltype, "sbj_001")
fdata$exprmats2$bsln <- fdata$exprmats$xsbj$exprmats$bsln[[fdata$celltype]]
fdata$exprmats2$xsbj <- fdata$exprmats$xsbj$exprmats$sbjv[[fdata$celltype]]
fdata$exprmats2$diluted <- fdata$exprmats$xbulk$null
fdata$exprmats2$ccv <- fdata$exprmats$xbulk$ccv

fdata$coexmats <- fdata$exprmats2 %>% session$collectionUtils$lapply(workspace$utils$computeCoexmat)

fdata$coexmats <- fdata$exprmats2 %>% session$collectionUtils$lapply(workspace$utils$computeCoexmat)

fdata$coextbls <- fdata$coexmats %>% session$collectionUtils$lapplyWithName(function(lvl, coexmat) {
  tibble(pair_id = workspace$utils$getPairIds(coexmat), coex = workspace$utils$vectorize(coexmat)) %>% 
    mutate(level = lvl)
}) %>% session$dataWrangler$rbind()

fdata$coextbls %>% session$graphingUtils$ggplot(aes(x = coex)) + geom_density(aes(color = level))

fdata$coextbls %>% group_by(level) %>% summarize(coex = mean(coex))
fdata$coextbls %>% group_by(level) %>% summarize(coex_max = max(coex))

# +++++++++=========== write to disk
fdata$coextbls %>% saveRDS(paste0(workspace$outputDir, "coextbls_fate_scenario_4.rds"))
# ===================================

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_04.png" 

xdata$means <- fdata$coextbls %>% group_by(level) %>% summarize(coex = mean(coex))

set.seed(1); xdata$pairsSample <- fdata$coextbls$pair_id %>% unique() %>% sample(1000)

xdata$plot <- fdata$coextbls %>% 
  filter(pair_id %in% xdata$pairsSample) %>% 
  session$graphingUtils$ggplot(aes(x = level, y = coex)) + 
  geom_line(aes(group = pair_id), size = 0.05) + 
  geom_line(data = xdata$means, aes(group = 1), color = "red", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  xlim(c("xcel", "bsln", "xsbj", "diluted", "ccv")) + 
  session$graphingUtils$tiltX(angle = 90) + 
  ylim(-1, 1) +
  ylab("Co-expression (Pearson's r)") +
  ggtitle("Scenario 5: genes co-regulated and synchronized\namong all cell types") + 
  theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 14), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 20, height = 15) 
# =================================================


# +++++++++++++++++++++++====================
# Scenario 5

fdata <- list(); gc()

fdata$simObjs <- read_rds(paste0(workspace$outputDir, "simObjs_fate_scenario_5.rds"))
fdata$exprmats$xsbj <- pdata$utils$xSbjExprmats(fdata$simObjs)
fdata$exprmats$xbulk <- pdata$utils$xbulkExprmats(fdata$exprmats$xsbj$exprmats$sbjv, cv = 0.3)

# compute co-expressions
fdata$coexmats <- fdata$simObjs$cellsBsln$simdat$samples

# exprmat at different levesl for a single cell types

fdata$celltype <- "excitatory"

fdata$exprmats2$xcel <- fdata$simObjs$cellsBsln$simdat$exprmat %>% 
  pdata$utils$extracCelLvlExprmat(fdata$simObjs$cellsBsln$simdat$samples, fdata$celltype, "sbj_001")
fdata$exprmats2$bsln <- fdata$exprmats$xsbj$exprmats$bsln[[fdata$celltype]]
fdata$exprmats2$xsbj <- fdata$exprmats$xsbj$exprmats$sbjv[[fdata$celltype]]
fdata$exprmats2$diluted <- fdata$exprmats$xbulk$null
fdata$exprmats2$ccv <- fdata$exprmats$xbulk$ccv

fdata$coexmats <- fdata$exprmats2 %>% session$collectionUtils$lapply(workspace$utils$computeCoexmat)

fdata$coexmats <- fdata$exprmats2 %>% session$collectionUtils$lapply(workspace$utils$computeCoexmat)

fdata$coextbls <- fdata$coexmats %>% session$collectionUtils$lapplyWithName(function(lvl, coexmat) {
  tibble(pair_id = workspace$utils$getPairIds(coexmat), coex = workspace$utils$vectorize(coexmat)) %>% 
    mutate(level = lvl)
}) %>% session$dataWrangler$rbind()

fdata$coextbls %>% session$graphingUtils$ggplot(aes(x = coex)) + geom_density(aes(color = level))

fdata$coextbls %>% group_by(level) %>% summarize(coex = mean(coex))
fdata$coextbls %>% group_by(level) %>% summarize(coex_max = max(coex))

# +++++++++=========== write to disk
fdata$coextbls %>% saveRDS(paste0(workspace$outputDir, "coextbls_fate_scenario_5.rds"))
# ===================================

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_05.png" 

xdata$means <- fdata$coextbls %>% group_by(level) %>% summarize(coex = mean(coex))

set.seed(1); xdata$pairsSample <- fdata$coextbls$pair_id %>% unique() %>% sample(1000)

xdata$plot <- fdata$coextbls %>% 
  filter(pair_id %in% xdata$pairsSample) %>% 
  session$graphingUtils$ggplot(aes(x = level, y = coex)) + 
  geom_line(aes(group = pair_id), size = 0.05) + 
  geom_line(data = xdata$means, aes(group = 1), color = "red", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  xlim(c("xcel", "bsln", "xsbj", "diluted", "ccv")) + 
  session$graphingUtils$tiltX(angle = 90) + 
  ylim(-1, 1) +
  ylab("Co-expression (Pearson's r)") +
  ggtitle("Scenario 1: genes not co-regulated") + 
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 14), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 30, height = 15) 
# =================================================


# =================================================
# OUTPUT S_FILE
xdata <- list()

xdata$fileName <- "sfile_10_fate_coex_values.tsv" 

xdata$files <- paste0("coextbls_fate_scenario_", 1:5, ".rds") 
xdata$files <- xdata$files %>% session$dataWrangler$attachNames(paste0("scenario_", c(2:5, 1))) # match the scenario # in the manuscript

xdata$dat <- xdata$files %>% session$collectionUtils$lapplyWithName(function(scenario, file) {
  tbl <- read_rds(paste0(workspace$outputDir, file))
  tbl %>% mutate(scenario = scenario)
}) %>% session$dataWrangler$rbind()

# change "bsln" to "csiv"
xdata$dat <- xdata$dat %>% 
  mutate(level = level %>% sapply(function(lvl) { if (lvl == "bsln") { "csiv" } else { lvl } }))

xdata$dat <- xdata$dat %>% arrange(scenario)

# some code for sanity check - that the scenarios demonstrate the right observations
# xdata$dat %>% group_by(scenario, level) %>% summarize(mean(coex)) -> x
# x %>% filter(scenario == "scenario_5")

write.table(xdata$dat, file = paste0(pdata$sfilesDir, xdata$fileName), 
            sep = "\t", row.names = FALSE) 
# =================================================






