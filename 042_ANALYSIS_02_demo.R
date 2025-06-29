pdata <- list()

pdata$figuresDir <- paste0(workspace$workspaceDir, "042_ANALYSIS_02_demo_FIGURES/")
pdata$sfilesDir <- paste0(workspace$workspaceDir, "042_S_FILES/")

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

# no coex program defined for the first demo version
fdata$this <- fdata$this %>% fdata$simulator$api$setRefDatCoexPrgms(list(), level = "sbj") 
fdata$this <- fdata$this %>% fdata$simulator$api$setRefDatCoexPrgms(list(), level = "cel") 

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
fdata$this %>% saveRDS(paste0(workspace$outputDir, "this_demo.rds"))
# ================================= 

# # "MATK"         "DLX6-AS1"     "LOC101927745" "PHLDB2"       "TRIM54"       "EMBP1"       
fdata$this$mdParams$cel$astrocyte["MATK", ]
fdata$this$mdParams$cel$excitatory["DLX6-AS1", ]


# in this script, demonstrate that I'm able to 
# 1. reproduce the mean-variance relationship at both the cell and the subject levels
# 2. impose correlations at the both the cell level and the subject level; which is absent in the real data
# 3. the distributions are integrated (for 1 gene?)

# do this all in just 1 cell type: excitatory neurons?
# 1000 genes; 1000 cells; 100 subjects; 200 at the cell level and 200 at the subject level

# +++++++++++++++++++++++++============
# first, add in the coex programs

# define coex programs (300 of 1k at the cell level vs. 300 at the subject level)

fdata <- list(); gc()

fdata$this <- read_rds(paste0(workspace$outputDir, "this_demo.rds"))
fdata$cteprf <- read_rds(paste0(workspace$outputDir, "seed_cteprf.rds"))$cteprf

fdata$celltypes <- colnames(fdata$cteprf)
names(fdata$celltypes) <- fdata$celltypes

fdata$genes <- rownames(fdata$cteprf)

set.seed(0); fdata$genesCoexCell <- fdata$genes %>% sample(500, replace = FALSE)
set.seed(2); fdata$genesCoexSbj <- fdata$genes %>% sample(500, replace = FALSE)

fdata$prgmCell <- list(program = list(
  celltypes = c("excitatory"), 
  genes = fdata$genesCoexCell, 
  value = 0.9
))

fdata$prgmSbj <- list(program = list(
  celltypes = c("excitatory"), 
  genes = fdata$genesCoexSbj, 
  value = 0.9
))

fdata$this <- fdata$this %>% workspace$api$setRefDatCoexPrgms(fdata$prgmCell, level = "cel")
fdata$this <- fdata$this %>% workspace$api$setRefDatCoexPrgms(fdata$prgmSbj, level = "sbj")

# =================================  
# COMMIT ==========================
fdata$this %>% saveRDS(paste0(workspace$outputDir, "this_demo.rds"))
# =================================  

# ++++++++++++++++++
# now simulate cells

fdata <- list()

fdata$this <- read_rds(paste0(workspace$outputDir, "this_demo.rds"))

fdata$cellsBsln <- fdata$this %>% workspace$api$simBseLnCels(nSubject = 100, nCell = 1000) 

fdata$sbjMeans <- fdata$this %>% workspace$api$simSbjLvMeans(100)

fdata$celParamsSbjv <- fdata$this %>% workspace$api$computeCelParams(fdata$sbjMeans$simSbjs) 

fdata$cellsSbjv <- fdata$cellsBsln$simCells %>% workspace$api$convertCelLvDist(fdata$cellsBsln$params$md, fdata$celParamsSbjv)

# =================================  
# COMMIT ==========================
fdata$simObjs <- list()
fdata$simObjs$cellsBsln <- list(simdat = fdata$cellsBsln$simCells, params = fdata$cellsBsln$params)
fdata$simObjs$cellsSbjv <- list(simdat = fdata$cellsSbjv, params = fdata$celParamsSbjv)
fdata$simObjs$subjects <- list(simdat = fdata$sbjMeans$simSbjs, params = fdata$sbjMeans$params)

fdata$simObjs %>% saveRDS(paste0(workspace$outputDir, "sim_demo_simObjs.rds"))
# =================================  

# SOFT COMMIT ++++++++++++++++++++++++++++++++++============
pdata$this <- read_rds(paste0(workspace$outputDir, "this_demo.rds"))
pdata$simObjs <- read_rds(paste0(workspace$outputDir, "sim_demo_simObjs.rds"))
# ++++++++++++++++++++++++++++++++++++++=

# +++++++++++++++++++++++++++++++++++++==================
# let's write out a few supplementary material files here

# =================================================
# OUTPUT S_FILE
xdata <- list()

xdata$fileName <- "sfile_01_cell_type_profiles_seed.tsv" 

xdata$dat <- pdata$this$refdat$cteprf %>% session$dataWrangler$setRownameAsColumn("gene_id")

write.table(xdata$dat, file = paste0(pdata$sfilesDir, xdata$fileName), 
            sep = "\t", row.names = FALSE) 
# =================================================

# =================================================
# OUTPUT S_FILE
xdata <- list()

xdata$fileName <- "sfile_02_gamma_model_parameters.tsv" 

xdata$dat <- pdata$this$mdParams %>% session$collectionUtils$lapplyWithName(function(currLvl, params) {
  params %>% session$collectionUtils$lapplyWithName(function(currCelltype, param) {
    param %>% session$dataWrangler$setRownameAsColumn("gene_id") %>% mutate(cell_type = currCelltype)
  }) %>% session$dataWrangler$rbind() %>% mutate(level = currLvl)
}) %>% session$dataWrangler$rbind()

write.table(xdata$dat, file = paste0(pdata$sfilesDir, xdata$fileName), 
            sep = "\t", row.names = FALSE)
# =================================================

# ++++++++++++++++++++===++++++++++++==========
# first examine the baseline cells & sbjv cells as well
# 1. mean variance relationship is maintained
# 2. mean expression for each gene is preserved # this is true with variation introduced by quantile normalization
# 3. co-expression is imposed successfully 

# baseline cell level mean & variance

fdata <- list()

# get mean / variance for simulated cell level data

xdata <- list()
xdata$celltype <- "excitatory"
xdata$subject <- "sbj_006"
xdata$samples <- pdata$simObjs$cellsBsln$simdat$samples %>% filter(cell_type == xdata$celltype, subject == xdata$subject)
xdata$exprmat <- pdata$simObjs$cellsBsln$simdat$exprmat[, xdata$samples$sample]

fdata$mm$celSim <- xdata$exprmat %>% 
  apply(1, function(vals) { c(m = mean(vals), v = var(vals)) }) %>% 
  t() %>% 
  session$dataWrangler$setRownameAsColumn("gene")

fdata$mm$celSim <- fdata$mm$celSim %>% mutate(level = "cel", type = "simulated")

# get mean / variance for real cell level data

fdata$mm$celReal <- pdata$this$refdat$exprmat$cel %>% 
  apply(1, function(vals) { c(m = mean(vals), v = var(vals)) }) %>% 
  t() %>% 
  session$dataWrangler$setRownameAsColumn("gene")

fdata$mm$celReal <- fdata$mm$celReal %>% mutate(level = "cel", type = "real")

set.seed(0); fdata$mm$celRealSample <- fdata$mm$celReal %>% filter(gene %in% sample(gene, 3000))

# subject level simulated

xdata <- list()
xdata$celltype <- "excitatory"
xdata$samples <- pdata$simObjs$subjects$simdat$samples %>% filter(cell_type == xdata$celltype)
xdata$exprmat <- pdata$simObjs$subjects$simdat$exprmat[, xdata$samples$sample]

fdata$mm$sbjSim <- xdata$exprmat %>% 
  apply(1, function(vals) { c(m = mean(vals), v = var(vals)) }) %>% 
  t() %>% 
  session$dataWrangler$setRownameAsColumn("gene")

fdata$mm$sbjSim <- fdata$mm$sbjSim %>% mutate(level = "sbj", type = "simulated")

# subject level real

fdata$mm$sbjReal <- pdata$this$refdat$exprmat$sbj %>% 
  apply(1, function(vals) { c(m = mean(vals), v = var(vals)) }) %>% 
  t() %>% 
  session$dataWrangler$setRownameAsColumn("gene")

fdata$mm$sbjReal <- fdata$mm$sbjReal %>% mutate(level = "sbj", type = "real")

set.seed(0); fdata$mm$sbjRealSample <- fdata$mm$sbjReal %>% filter(gene %in% sample(gene, 3000))

# ==========

fdata$main <- fdata$mm$celSim %>% 
  rbind(fdata$mm$celRealSample) %>% 
  rbind(fdata$mm$sbjSim) %>% 
  rbind(fdata$mm$sbjRealSample)


# COMMIT =============
fdata$main %>% saveRDS(paste0(workspace$outputDir, "02_demo_mv.rds")) 
# ======================


# =================================================
# OUTPUT S_FILE
xdata <- list(); xdata$fileName <- "sfile_03_mv_trend.tsv" 

xdata$dat <- read_rds(paste0(workspace$outputDir, "02_demo_mv.rds"))

write.table(xdata$dat, file = paste0(pdata$sfilesDir, xdata$fileName), 
            sep = "\t", row.names = FALSE) 
# =================================================


# ++++++++++============================
# work on the plots showing the hierarchical distribution nature between cells vs. subjects

fdata <- list(); gc()

fdata$main <- read_rds(paste0(workspace$outputDir, "02_demo_mv.rds"))

fdata$main %>% session$graphingUtils$ggplot(aes(x = m)) + geom_density()

# let's try this: PCLO

fdata$gene <- "TSPAN7"
fdata$cells <- pdata$simObjs$cellsSbjv$simdat$samples %>% filter(cell_type == "excitatory") %>% group_by(subject) %>% ungroup()
fdata$exprCels <- pdata$simObjs$cellsSbjv$simdat$exprmat[fdata$gene, fdata$cells$sample]
fdata$exprCels <- fdata$exprCels %>% session$dataWrangler$vectorToTibble() %>% dplyr::select(sample = variable, expr = value)
fdata$exprCels <- fdata$cells %>% left_join(fdata$exprCels, by = "sample")

fdata$main <- fdata$exprCels

fdata$sbjMeans <- fdata$main %>% group_by(subject) %>% summarize(m_sbjv = mean(expr))


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_09.png" # vector file too big just draw bitmap

xdata$plot <- fdata$main %>% 
  session$graphingUtils$ggplot(aes(x = expr)) + 
  geom_density(aes(group = subject), size = 0.2, color = "#bf8502") + 
  geom_rug(data = fdata$sbjMeans, aes(x = m_sbjv), size = 0.7, color = "#0576b5") +
  xlim(0, 6500) + 
  xlab("Expression (CPM)") + 
  ylab("Density") + 
  ggtitle("Distributions of xCell level expression of sim-TSPAN7")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 30, height = 15) 
# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_10.png" # vector file too big just draw bitmap

xdata$plot <- fdata$sbjMeans %>% 
  session$graphingUtils$ggplot(aes(x = m_sbjv)) + 
  geom_density(color = "#0576b5") + 
  geom_rug(data = fdata$sbjMeans, aes(x = m_sbjv), size = 0.7, color = "#0576b5") +
  xlim(0, 6500) + 
  xlab("Expression (CPM)") + 
  ylab("Density") + 
  ggtitle("Distribution of xSubjct level expression of sim-TSPAN7")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 30, height = 15) 
# =================================================






# ++++++++++============================
# plot the m-v relationship scatter plots, at both the cel and sbj levels

fdata <- list(); gc()

fdata$main <- read_rds(paste0(workspace$outputDir, "02_demo_mv.rds"))

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_01.png" # vector file too big just draw bitmap

xdata$plot <- fdata$main %>% 
  filter(type == "real") %>% 
  session$graphingUtils$ggplot(aes(x = log2(m + 1), y = log2(v + 1))) + 
  geom_point(aes(color = level), shape = 1) + 
  geom_smooth(aes(group = level), method = "lm", se = FALSE, fill = NA,
              formula = y ~ poly(x, 3, raw = TRUE), colour = "black", size = 0.5) +
  ylim(0, 22) +
  ylab("Variance (log2)") + 
  theme(panel.spacing = unit(2, "lines")) + 
  ggtitle("Real data") + 
  scale_color_manual(values = c("#fac243", "#56B4E9")) + 
  theme(legend.position = "none") + 
  scale_x_continuous(name = "Mean (log2)", limits = c(0, 13)) 
  
xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 15, height = 15) 
# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_02.png" # vector file too big just draw bitmap

xdata$plot <- fdata$main %>% 
  filter(type == "simulated") %>% 
  session$graphingUtils$ggplot(aes(x = log2(m + 1), y = log2(v + 1))) + 
  geom_point(aes(color = level), shape = 1) + 
  geom_smooth(aes(group = level), method = "lm", se = FALSE, fill = NA,
              formula = y ~ poly(x, 3, raw = TRUE), colour = "black", size = 0.5) +
  ylim(0, 22) +
  ylab("Variance (log2)") + 
  theme(panel.spacing = unit(2, "lines")) + 
  ggtitle("Simulated data") + 
  scale_color_manual(values = c("#bf8502", "#0576b5")) + 
  theme(legend.position = "none") +
  scale_x_continuous(name = "Mean (log2)", limits = c(0, 13)) 

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 15, height = 15) 
# =================================================


# ++++++++++++++++++++++============
# let's check out some marginal distributions

# comparable m and v in real and simulated:
# 4 EIF2S3        50.2 3687. cel   simulated
# 5 RPL41         51.8 3845. cel   real  

# 8 MED19         51.3 4353. cel   simulated
# 9 LINGO1        50.9 4443. cel   real  

fdata <- list()

fdata$main <- read_rds(paste0(workspace$outputDir, "02_demo_mv.rds"))

# print the means > 50 to pick comparable genes
fdata$main %>% filter(level == "cel", m > 50) %>% arrange(v)

fdata$geneSim <- "MED19"
fdata$geneReal <- "LINGO1"

# expression & ranking for the cell level simulated & real data

xdata <- list()
xdata$celltype <- "excitatory"
xdata$subject <- "sbj_006"
xdata$samples <- pdata$simObjs$cellsBsln$simdat$samples %>% filter(cell_type == xdata$celltype, subject == xdata$subject)
xdata$exprmat <- pdata$simObjs$cellsBsln$simdat$exprmat[, xdata$samples$sample]

xdata$exprGene <- xdata$exprmat[fdata$geneSim, ] %>% # right, only 1k cells here!!!
  session$dataWrangler$vectorToTibble() %>% 
  mutate(rank = rank(value, ties.method = "random")) %>% 
  mutate(rank = rank / max(rank)) %>% 
  mutate(type = "simulated")

# FDATA COMMIT ====
fdata$exprGene <- xdata$exprGene
# =================

# now the real version

xdata <- list()

xdata$exprGene <- pdata$this$refdat$exprmat$cel[fdata$geneReal, ] %>% 
  session$dataWrangler$vectorToTibble() %>% 
  filter(variable %in% sample(variable, nrow(fdata$exprGene))) # also use 1k genes here so its fair comparison

xdata$exprGene <- xdata$exprGene %>% 
  mutate(rank = rank(value, ties.method = "random")) %>% 
  mutate(rank = rank / max(rank)) %>% 
  mutate(type = "real")

# FDATA COMMIT ====
fdata$exprGene <- fdata$exprGene %>% rbind(xdata$exprGene)
# =================

# cumulative probability distribution
fdata$exprGene %>%
  session$graphingUtils$ggplot(aes(x = log2(value + 1), y = rank)) + 
  geom_point(aes(color = type), size = 0.5)

# probability desntiy distribution
fdata$exprGene %>% 
  session$graphingUtils$ggplot(aes(x = log2(value + 1))) +
  geom_density(aes(color = type), size = 1)

fdata$mvSim <- fdata$main %>% filter(gene == fdata$geneSim, level == "cel", type == "simulated")
fdata$mvReal <- fdata$main %>% filter(gene == fdata$geneReal, level == "cel", type == "real")

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_03.png" # vector file too big just draw bitmap

xdata$plot <- fdata$exprGene %>%
  session$graphingUtils$ggplot(aes(x = log2(value + 1), y = rank)) + 
  geom_point(aes(color = type), size = 0.5) + 
  ylab("Rank percentile") + 
  xlab("Expression (log2 - CPM)") + 
  ggtitle(paste0("CDF of xCell level expression\n", 
                 "Simulated: sim-MED19 - m = ", round(fdata$mvSim$m, 0), "; v = ", round(fdata$mvSim$v, 0), 
                 "\nReal: LINGO1 - m = ", round(fdata$mvReal$m, 0), "; v = ", round(fdata$mvReal$v, 0))) + 
  theme(plot.title = element_text(size = 15)) + 
  scale_color_manual(values = c("#fac243", "#bf8502")) + 
  theme(legend.position = "none")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 15, height = 15) 
# =================================================

# ++++++++++++++++++++++============
# let's check out some marginal distributions at the subject level this time

# comparable m and v in real and simulated:
# ENSG00000120008.10
# MOB1B

fdata <- list()

fdata$main <- read_rds(paste0(workspace$outputDir, "02_demo_mv.rds"))

# print the means > 50 to pick comparable genes
fdata$main %>% filter(level == "sbj", m > 50) %>% arrange(m)

fdata$geneSim <- "MOB1B"
fdata$geneReal <- "ENSG00000120008.10"

# expression & ranking for the cell level simulated & real data

xdata <- list()
xdata$celltype <- "excitatory"
xdata$samples <- pdata$simObjs$subjects$simdat$samples %>% filter(cell_type == xdata$celltype)
xdata$exprmat <- pdata$simObjs$subjects$simdat$exprmat[, xdata$samples$sample]

xdata$exprGene <- xdata$exprmat[fdata$geneSim, ] %>% # right, only 1k cells here!!!
  session$dataWrangler$vectorToTibble() %>% 
  mutate(rank = rank(value, ties.method = "random")) %>% 
  mutate(rank = rank / max(rank)) %>% 
  mutate(type = "simulated")

# FDATA COMMIT ====
fdata$exprGene <- xdata$exprGene
# =================

# now the real version

xdata <- list()

xdata$exprGene <- pdata$this$refdat$exprmat$sbj[fdata$geneReal, ] %>% 
  session$dataWrangler$vectorToTibble() %>% 
  filter(variable %in% sample(variable, nrow(fdata$exprGene))) # also use 1k genes here so its fair comparison

xdata$exprGene <- xdata$exprGene %>% 
  mutate(rank = rank(value, ties.method = "random")) %>% 
  mutate(rank = rank / max(rank)) %>% 
  mutate(type = "real")

# FDATA COMMIT ====
fdata$exprGene <- fdata$exprGene %>% rbind(xdata$exprGene)
# =================

# cumulative probability distribution
fdata$exprGene %>%
  session$graphingUtils$ggplot(aes(x = log2(value + 1), y = rank)) + 
  geom_point(aes(color = type), size = 0.5)

# probability desntiy distribution
fdata$exprGene %>% 
  session$graphingUtils$ggplot(aes(x = log2(value + 1))) +
  geom_density(aes(color = type), size = 1)

fdata$mvSim <- fdata$main %>% filter(gene == fdata$geneSim, level == "sbj", type == "simulated")
fdata$mvReal <- fdata$main %>% filter(gene == fdata$geneReal, level == "sbj", type == "real")

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_04.png" # vector file too big just draw bitmap

xdata$plot <- fdata$exprGene %>%
  session$graphingUtils$ggplot(aes(x = log2(value + 1), y = rank)) + 
  geom_point(aes(color = type), size = 1) + 
  ylab("Rank percentile") + 
  xlab("Expression (log2 - CPM)") + 
  ggtitle(paste0("CDF of xSubject level expression\n", 
                 "Simulated: sim-MOB1B - m = ", round(fdata$mvSim$m, 0), "; v = ", round(fdata$mvSim$v, 0), 
                 "\nReal: WDR11 - m = ", round(fdata$mvReal$m, 0), "; v = ", round(fdata$mvReal$v, 0))) + 
  theme(plot.title = element_text(size = 15)) + 
  scale_color_manual(values = c("#56B4E9", "#0576b5")) + 
  theme(legend.position = "none")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 15, height = 15) 
# =================================================

# check the co-expressions of the intended genes 

fdata <- list()

fdata$cel$geneA <- "TMX3" 
fdata$cel$geneB <- "KCNB1" 

fdata$sbj$geneA <- "ACTR6" 
fdata$sbj$geneB <- "ADD1" 

# work out cross cell co-expresison for these genes
xdata <- list()
xdata$celltype <- "excitatory"
xdata$subject <- "sbj_001"
xdata$samples <- pdata$simObjs$cellsBsln$simdat$samples %>% filter(cell_type == xdata$celltype, subject == xdata$subject)
xdata$exprmat <- pdata$simObjs$cellsBsln$simdat$exprmat[, xdata$samples$sample]
xdata$exprGenes$celPair <- xdata$exprmat %>% workspace$utils$getExprGenes(fdata$cel$geneA , fdata$cel$geneB)
xdata$exprGenes$sbjPair <- xdata$exprmat %>% workspace$utils$getExprGenes(fdata$sbj$geneA , fdata$sbj$geneB)

# FDATA COMMIT =====================
fdata$exprGenes$xCell <- xdata$exprGenes
# ================================


# work out cross cell co-expression for these genes
xdata <- list()
xdata$celltype <- "excitatory"
xdata$samples <- pdata$simObjs$subjects$simdat$samples %>% filter(cell_type == xdata$celltype)
xdata$exprmat <- pdata$simObjs$subjects$simdat$exprmat[, xdata$samples$sample]
xdata$exprGenes$celPair <- xdata$exprmat %>% workspace$utils$getExprGenes(fdata$cel$geneA , fdata$cel$geneB)
xdata$exprGenes$sbjPair <- xdata$exprmat %>% workspace$utils$getExprGenes(fdata$sbj$geneA , fdata$sbj$geneB)

# FDATA COMMIT =====================
fdata$exprGenes$xSubject <- xdata$exprGenes
# ================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_05.png" # vector file too big just draw bitmap

xdata$main <- fdata$exprGenes$xCell$celPair

xdata$coex <- cor(xdata$main$gene_a, xdata$main$gene_b)

xdata$plot <- xdata$main %>%
  session$graphingUtils$ggplot(aes(x = gene_a, y = gene_b)) + 
  geom_point(color = "#bf8502") + 
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) + 
  ylab(paste0(paste0("sim-", fdata$cel$geneA), " (CPM)")) + 
  xlab(paste0(paste0("sim-", fdata$cel$geneB), " (CPM)")) + 
  ggtitle(paste0("xCell co-expression"), paste0("Pearson's r = ", round(xdata$coex, 2))) + 
  theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 14))  

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 15, height = 15) 
# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_06.png" # vector file too big just draw bitmap

xdata$main <- fdata$exprGenes$xCell$sbjPair

xdata$coex <- cor(xdata$main$gene_a, xdata$main$gene_b)

xdata$plot <- xdata$main %>%
  session$graphingUtils$ggplot(aes(x = gene_a, y = gene_b)) + 
  geom_point(color = "#bf8502") + 
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) + 
  ylab(paste0("sim-", fdata$sbj$geneA, " (CPM)")) + 
  xlab(paste0("sim-", fdata$sbj$geneB, " (CPM)")) + 
  ggtitle(paste0("xCell co-expression"), paste0("Pearson's r = ", round(xdata$coex, 2))) + 
  theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 14))  

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 15, height = 15) 
# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_07.png" # vector file too big just draw bitmap

xdata$main <- fdata$exprGenes$xSubject$celPair

xdata$coex <- cor(xdata$main$gene_a, xdata$main$gene_b)

xdata$plot <- xdata$main %>%
  session$graphingUtils$ggplot(aes(x = gene_a, y = gene_b)) + 
  geom_point(color = "#0576b5") + 
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) + 
  ylab(paste0(paste0("sim-", fdata$cel$geneA), " (CPM)")) + 
  xlab(paste0(paste0("sim-", fdata$cel$geneB), " (CPM)")) + 
  ggtitle(paste0("xSubject co-expression"), paste0("Pearson's r = ", round(xdata$coex, 2))) + 
  theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 14))  

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 15, height = 15) 
# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_08.png" # vector file too big just draw bitmap

xdata$main <- fdata$exprGenes$xSubject$sbjPair

xdata$coex <- cor(xdata$main$gene_a, xdata$main$gene_b)

xdata$plot <- xdata$main %>%
  session$graphingUtils$ggplot(aes(x = gene_a, y = gene_b)) + 
  geom_point(color = "#0576b5") + 
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) + 
  ylab(paste0("sim-", fdata$sbj$geneA, " (CPM)")) + 
  xlab(paste0("sim-", fdata$sbj$geneB, " (CPM)")) + 
  ggtitle(paste0("xSubject co-expression"), paste0("Pearson's r = ", round(xdata$coex, 2))) + 
  theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 14))  

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 15, height = 15) 
# =================================================

# =====
# below is the code for showing co-expression distribution, which at this point is probably not needed


fdata$coexProgram$cel <- pdata$simObjs$cellsBsln$params$coexPrograms$program
fdata$coexProgram$sbj <- pdata$simObjs$subjects$params$coexPrograms$program

fdata$celltype <- "excitatory"
fdata$subject <- "sbj_001"
fdata$samples <- pdata$simObjs$cellsBsln$simdat$samples %>% filter(cell_type == fdata$celltype, subject == fdata$subject)
fdata$exprmat <- pdata$simObjs$cellsBsln$simdat$exprmat[, fdata$samples$sample]
fdata$coexmat <- fdata$exprmat %>% workspace$utils$computeCoexmat()
fdata$coextbl <- fdata$coexmat %>% 
  workspace$utils$getCoords() %>% 
  mutate(coex = fdata$coexmat %>% workspace$utils$vectorize())

fdata$coextbl <- fdata$coextbl %>% 
  mutate(gene_a_coex = (gene_a %in% fdata$coexProgram$cel$genes), 
         gene_b_coex = (gene_b %in% fdata$coexProgram$cel$genes)) %>% 
  mutate(coex_cel = (as.numeric(gene_a_coex) + as.numeric(gene_b_coex)) == 2) %>% 
  dplyr::select(-gene_a_coex, -gene_b_coex) %>% 
  mutate(gene_a_coex = (gene_a %in% fdata$coexProgram$sbj$genes), 
         gene_b_coex = (gene_b %in% fdata$coexProgram$sbj$genes)) %>% 
  mutate(coex_sbj = (as.numeric(gene_a_coex) + as.numeric(gene_b_coex)) == 2) %>% 
  dplyr::select(-gene_a_coex, -gene_b_coex) %>% 
  mutate(coex_prgm = paste0(coex_cel, ".", coex_sbj)) %>% 
  mutate(coex_prgm = coex_prgm %>% sapply(function(str) {
    if (str == "TRUE.TRUE") { "both" }
    else if (str == "TRUE.FALSE") { "cell" }
    else if (str == "FALSE.TRUE") { "sbj" }
    else { "none" }
  }))

fdata$coextbl %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_density(aes(color = coex_prgm)) + 
  ggtitle("", "xCell coexpression distirbution") +
  xlim(-1, 1)


fdata$exprmat %>% workspace$utils$plotScatter("TMX3", "KCNB1")

# check the co-expressions of the intended genes at the correct levels

fdata <- list()

fdata$coexProgram$cel <- pdata$simObjs$cellsBsln$params$coexPrograms$program
fdata$coexProgram$sbj <- pdata$simObjs$subjects$params$coexPrograms$program

fdata$celltype <- "excitatory"
fdata$samples <- pdata$simObjs$subjects$simdat$samples %>% filter(cell_type == fdata$celltype)
fdata$exprmat <- pdata$simObjs$subjects$simdat$exprmat[, fdata$samples$sample]
fdata$coexmat <- fdata$exprmat %>% workspace$utils$computeCoexmat()
fdata$coextbl <- fdata$coexmat %>% 
  workspace$utils$getCoords() %>% 
  mutate(coex = fdata$coexmat %>% workspace$utils$vectorize())

fdata$coextbl <- fdata$coextbl %>% 
  mutate(gene_a_coex = (gene_a %in% fdata$coexProgram$cel$genes), 
         gene_b_coex = (gene_b %in% fdata$coexProgram$cel$genes)) %>% 
  mutate(coex_cel = (as.numeric(gene_a_coex) + as.numeric(gene_b_coex)) == 2) %>% 
  dplyr::select(-gene_a_coex, -gene_b_coex) %>% 
  mutate(gene_a_coex = (gene_a %in% fdata$coexProgram$sbj$genes), 
         gene_b_coex = (gene_b %in% fdata$coexProgram$sbj$genes)) %>% 
  mutate(coex_sbj = (as.numeric(gene_a_coex) + as.numeric(gene_b_coex)) == 2) %>% 
  dplyr::select(-gene_a_coex, -gene_b_coex) %>% 
  mutate(coex_prgm = paste0(coex_cel, ".", coex_sbj)) %>% 
  mutate(coex_prgm = coex_prgm %>% sapply(function(str) {
    if (str == "TRUE.TRUE") { "both" }
    else if (str == "TRUE.FALSE") { "cell" }
    else if (str == "FALSE.TRUE") { "sbj" }
    else { "none" }
  }))

fdata$coextbl %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_density(aes(color = coex_prgm)) + 
  ggtitle("", "xSbj coexpression distirbution") + 
  xlim(-1, 1)


fdata$exprmat %>% workspace$utils$plotScatter("ACTR6", "ADD1")
fdata$exprmat %>% workspace$utils$plotScatter("TMX3", "KCNB1")

