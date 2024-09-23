# +++++++++++++++++++++
# In this script, I'm performing the PSDBK-BULK analysis
# First, show that I can predict it singificantly, and everything passes sanity checks, etc
# Show that highly cell type specific things could travel through (otherwise things are spread thin)
# Show that synchrony is the way things could pass through - use ribosomal protein genes as a case study
# generalize the synchrony finding

pdata <- list()

pdata$figuresDir <- paste0(workspace$workspaceDir, "041_ANALYSIS_06_dilution_FIGURES/")

pdata$clusters <- read_rds(paste0(workspace$outputDir, paste0("sc_clusters_objects.rds")))
pdata$clusters <- pdata$clusters$clusterFlats %>% unname() %>% unlist(recursive = FALSE)

pdata$genes <- readRDS(paste0(workspace$outputDir, "genes_metadata.rds")) # genes

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


pdata$files$lnktbls <- tibble(outfile = list.files(workspace$outputDir)) %>% 
  filter(grepl("lnktbls", outfile)) %>% 
  mutate(outfile_name = gsub(".rds", "", outfile)) %>% 
  mutate(strs = strsplit(outfile_name, "_")) %>% 
  mutate(level = strs %>% sapply(function(currStr) { currStr[1] })) %>% 
  mutate(type = strs %>% sapply(function(currStr) { currStr[2] })) %>% 
  mutate(cell_type = strs %>% sapply(function(currStr) { currStr[3] })) %>% 
  dplyr::select(outfile, level, type, cell_type)

pdata$files$exprmats <- tibble(outfile = list.files(workspace$outputDir)) %>% 
  filter(grepl("exprmats", outfile)) %>% 
  mutate(outfile_name = gsub(".rds", "", outfile)) %>% 
  mutate(strs = strsplit(outfile_name, "_")) %>% 
  mutate(level = strs %>% sapply(function(currStr) { currStr[1] })) %>% 
  mutate(type = strs %>% sapply(function(currStr) { currStr[2] })) %>% 
  mutate(dataset = strs %>% sapply(function(currStr) { currStr[3] })) %>% 
  dplyr::select(outfile, level, type, dataset)

pdata$files$uexprdats <- tibble(outfile = list.files(workspace$outputDir)) %>% 
  filter(grepl("uexprdats", outfile)) %>% 
  mutate(outfile_name = gsub(".rds", "", outfile)) %>% 
  mutate(strs = strsplit(outfile_name, "_")) %>% 
  mutate(type = strs %>% sapply(function(currStr) { currStr[1] })) %>% 
  mutate(dataset = strs %>% sapply(function(currStr) { currStr[2] })) %>% 
  dplyr::select(outfile, type, dataset)

pdata$ctprofiles <- read_rds(paste0(workspace$outputDir, "sc_stats_ctprofiles.rds")) 




# +++++++++++++++++++++
# generate null distirbutiosn for comparison!!

# for Velmeshev

fdata <- list(); gc()

fdata$psdbkres <- c("velmeshev", "rosmap") %>% session$dataWrangler$attachNames()
  
fdata$psdbkres <- fdata$psdbkres %>% session$collectionUtils$lapply(function(dataset) { 
  read_rds(paste0(workspace$outputDir, paste0("bk_exprmats_", dataset,"-psdbkres.rds"))) 
})

fdata$lmstats <- fdata$psdbkres %>% lapply(function(obj) { obj$lmModel$lmStats })
fdata$coefstats <- fdata$psdbkres %>% lapply(function(obj) { obj$lmModel$coefStats })
fdata$sigGenes <- fdata$lmstats %>% lapply(function(lmstats) { lmstats %>% filter(qvalue < 0.1) })

fdata$models <- workspace$utils$fitModels(fdata$exprdat$bk$exprmat, fdata$exprmats$sbj) 


# generate null distributions

workspace$utils$fitModelsSimple(fdata$psdbkres$velmeshev$exprmat$bkOrig, 
                                fdata$psdbkres$velmeshev$exprmat$sbj) 


fdata$null <- list()

# generate null distribution # do this later and for rosmap as well
fdata$null$velmeshev <- paste0("null_", formatC(1:100, digits = 2, flag = "0")) %>% session$collectionUtils$lapply(function(i) {
  currExprmatBk <- fdata$psdbkres$velmeshev$exprmat$bkOrig
  colnames(currExprmatBk) <- sample(colnames(currExprmatBk), replace = FALSE)
  currLms <- workspace$utils$fitModelsSimple(currExprmatBk, fdata$psdbkres$velmeshev$exprmat$sbj)
  currLms <- currLms %>% mutate(iteration = i)
  return(currLms)
}) %>% session$dataWrangler$rbind()

fdata$null$rosmap <- paste0("null_", formatC(1:100, digits = 2, flag = "0")) %>% session$collectionUtils$lapply(function(i) {
  currExprmatBk <- fdata$psdbkres$rosmap$exprmat$bkOrig
  colnames(currExprmatBk) <- sample(colnames(currExprmatBk), replace = FALSE)
  currLms <- workspace$utils$fitModelsSimple(currExprmatBk, fdata$psdbkres$rosmap$exprmat$sbj)
  currLms <- currLms %>% mutate(iteration = i)
  return(currLms)
}) %>% session$dataWrangler$rbind()

# ============COMMIT
fdata$null %>% saveRDS(paste0(workspace$outputDir, "psdbk_nulls.rds"))
# ++++++++++++++=

# +++++++++++++++++++++
# First, prove that things actually do work

# for Velmeshev

fdata <- list(); gc()

fdata$psdbkres <- c("velmeshev", "rosmap") %>% session$dataWrangler$attachNames()

fdata$psdbkres <- fdata$psdbkres %>% session$collectionUtils$lapply(function(dataset) { 
  read_rds(paste0(workspace$outputDir, paste0("bk_exprmats_", dataset,"-psdbkres.rds"))) 
})

fdata$lmstats <- fdata$psdbkres %>% lapply(function(obj) { obj$lmModel$lmStats })
fdata$coefstats <- fdata$psdbkres %>% lapply(function(obj) { obj$lmModel$coefStats })
fdata$sigGenes <- fdata$lmstats %>% lapply(function(lmstats) { lmstats %>% filter(qvalue < 0.1) })

fdata$nulls <- read_rds(paste0(workspace$outputDir, "psdbk_nulls.rds"))

# plot the rsqr & pvalue distributions to show that the signal is REAL


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_01.png"

(xdata$mean <- fdata$lmstats$velmeshev$rsqr %>% max())

xdata$plot <- fdata$nulls$velmeshev %>% 
  session$graphingUtils$ggplot(aes(x = rsqr)) + 
  geom_density(aes(group = iteration), linewidth = 0.02) + 
  geom_density(data = fdata$lmstats$velmeshev, color = "#07519C", linewidth = 1) + 
  geom_vline(xintercept = xdata$mean, color = "#07519C", linetype = "dashed") +
  xlab(bquote(R^2)) + 
  ylab("Density") +
  ggtitle("PVE distribution across genes", "Velmeshev")
  
xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 19, height = 19) 
# =================================================

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_02.png"

# fdata$nulls$rosmap %>% group_by(iteration) %>% summarize(rsqr = max(rsqr)) %>% 
#   summarize(rsqr = mean(rsqr))

xdata$plot <- fdata$nulls$velmeshev %>% 
  session$graphingUtils$ggplot(aes(x = pvalue)) + 
  geom_density(aes(group = iteration), linewidth = 0.03) + 
  geom_density(data = fdata$lmstats$velmeshev, color = "#07519C", linewidth = 1) + 
  xlab("P-value") + 
  ylab("Density") + 
  ggtitle("P-value distribution", "Velmeshev")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 19, height = 19) 
# =================================================

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_03.png"

(xdata$mean <- fdata$lmstats$rosmap$rsqr %>% mean())

xdata$plot <- fdata$nulls$rosmap %>% 
  session$graphingUtils$ggplot(aes(x = rsqr)) + 
  geom_density(aes(group = iteration), linewidth = 0.02) + 
  geom_density(data = fdata$lmstats$rosmap, color = "#07519C", linewidth = 1) + 
  geom_vline(xintercept = xdata$mean, color = "#07519C", linetype = "dashed") +
  xlab(bquote(R^2)) + 
  ylab("Density") +
  ggtitle("PVE distribution across genes", "ROSMAP")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 19, height = 19) 
# =================================================

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_04.png"

xdata$plot <- fdata$nulls$rosmap %>% 
  session$graphingUtils$ggplot(aes(x = pvalue)) + 
  geom_density(aes(group = iteration), linewidth = 0.03) + 
  geom_density(data = fdata$lmstats$rosmap, color = "#07519C", linewidth = 1) + 
  xlab("P-value") + 
  ylab("Density") + 
  xlim(0, 1) + 
  ggtitle("P-value distribution", "ROSMAP")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 19, height = 19) 
# =================================================


fdata$nulls$rosmap %>% group_by(gene_id) %>% summarize(rsqr_null_max = max(rsqr)) -> x
x %>% left_join(fdata$lmstats$rosmap %>% dplyr::select(gene_id, rsqr)) %>% 
  mutate(sig_em = rsqr - rsqr_null_max) %>% filter(sig_em > 0)

fdata$lmstats$rosmap
fdata$lmstats$velmeshev


# next, plot rsqr_inde vs rsqr_full

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_05.png"
xdata$figName2 <- "figure_05_1.png"

xdata$statsFull <- fdata$lmstats$velmeshev %>% dplyr::select(gene_id, rsqr_full = rsqr, aic_full = aic) 

xdata$statsIndep <- fdata$coefstats$velmeshev %>% 
  group_by(gene_id) %>% filter(rsqr_indep == max(rsqr_indep)) %>% ungroup() %>% 
  dplyr::select(gene_id, cell_type, rsqr_indep, aic_indep)

xdata$main <- xdata$statsFull %>% left_join(xdata$statsIndep, by = "gene_id")

xdata$plot <- xdata$main %>% mutate(cell_type = workspace$utils$fmtCelltypes(cell_type)) %>% 
  session$graphingUtils$ggplot(aes(x = rsqr_indep, y = rsqr_full)) + 
  geom_point(aes(color = cell_type), shape = 1, size = 3) + 
  geom_abline(slope = 1, linetype = "dashed") + 
  scale_color_brewer(palette = "Set1") +
  xlim(0, 1) + 
  ylim(0, 1) + 
  xlab(bquote("Maximal PVE by one cell type"~(R^2))) + 
  ylab(bquote("PVE by the all cell types"~(R^2))) +
  ggtitle("PVE by models using one vs.\nall cell types", "Velmeshev") + 
  theme(legend.position = "none")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 19, height = 19) 

xdata$plot2 <- xdata$main %>% 
  mutate(aic_diff = aic_full - aic_indep) %>% 
  session$graphingUtils$ggplot(aes(x = aic_diff)) + 
  geom_histogram() + geom_vline(xintercept = 0, linetype = "dashed") + 
  xlab("AIC (all cell types) - AIC (best cell type) ") + 
  ylab("Number of genes") +
  ggtitle("AIC difference between models\nusing all vs. one cell type", "Velmeshev")

xdata$plot2

ggsave(filename = paste0(pdata$figuresDir, xdata$figName2), 
       plot = xdata$plot2, device = "png", units = "cm", width = 19, height = 19) 


# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_06.png"
xdata$figName2 <- "figure_06_1.png"

xdata$statsFull <- fdata$lmstats$rosmap %>% dplyr::select(gene_id, rsqr_full = rsqr, aic_full = aic) 

xdata$statsIndep <- fdata$coefstats$rosmap %>% 
  group_by(gene_id) %>% filter(rsqr_indep == max(rsqr_indep)) %>% ungroup() %>% 
  dplyr::select(gene_id, cell_type, rsqr_indep, aic_indep)

xdata$main <- xdata$statsFull %>% left_join(xdata$statsIndep, by = "gene_id")

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = rsqr_indep, y = rsqr_full)) + 
  geom_point(aes(color = cell_type), shape = 1, size = 3) + 
  geom_abline(slope = 1, linetype = "dashed") + 
  scale_color_brewer(palette = "Set1") +
  xlim(0, 1) + 
  ylim(0, 1) + 
  xlab(bquote("Maximal PVE by one cell type"~(R^2))) + 
  ylab(bquote("PVE by the all cell types"~(R^2))) +
  ggtitle("PVE by models using one vs.\nall cell types", "ROSMAP") + 
  theme(legend.position = "none")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 19, height = 19) 

xdata$plot2 <- xdata$main %>% 
  mutate(aic_diff = aic_full - aic_indep) %>% 
  session$graphingUtils$ggplot(aes(x = aic_diff)) + 
  geom_histogram() + geom_vline(xintercept = 0, linetype = "dashed") + 
  xlab("AIC (all cell types) - AIC (best cell type) ") + 
  ylab("Number of genes") +
  ggtitle("AIC difference between models \nusing all vs. one cell type", "ROSMAP")

xdata$plot2

ggsave(filename = paste0(pdata$figuresDir, xdata$figName2), 
       plot = xdata$plot2, device = "png", units = "cm", width = 19, height = 19) 

# =================================================

# show that genes explained specifically by a given cell type will also have elevated expression in that cell type

fdata <- list(); gc()

fdata$psdbkres <- c("velmeshev", "rosmap") %>% session$dataWrangler$attachNames()

fdata$psdbkres <- fdata$psdbkres %>% session$collectionUtils$lapply(function(dataset) { 
  read_rds(paste0(workspace$outputDir, paste0("bk_exprmats_", dataset,"-psdbkres.rds"))) 
})

fdata$lmstats <- fdata$psdbkres %>% lapply(function(obj) { obj$lmModel$lmStats })
fdata$coefstats <- fdata$psdbkres %>% lapply(function(obj) { obj$lmModel$coefStats })
fdata$sigGenes <- fdata$lmstats %>% lapply(function(lmstats) { lmstats %>% filter(qvalue < 0.1) })

fdata$celltypes <- pdata$celltypes %>% sort()

# velmeshev

fdata$rsqrIndeps$velmeshev <- fdata$coefstats$velmeshev %>% 
  dplyr::select(gene_id, cell_type, rsqr_indep) %>% 
  spread(cell_type, rsqr_indep) %>% 
  session$dataWrangler$setColAsRownames("gene_id")

fdata$rsqrIndeps$velmeshev <- fdata$celltypes %>% lapply(function(celltype) {
  
  targetExpr <- fdata$rsqrIndeps$velmeshev[, celltype]
  othrExprMax <- fdata$rsqrIndeps$velmeshev[, colnames(fdata$rsqrIndeps$velmeshev) != celltype] %>% apply(1, max) 
  
  tibble(gene_id = rownames(fdata$rsqrIndeps$velmeshev)) %>% 
    mutate(rsqr = targetExpr, 
           rsqr_max = othrExprMax[gene_id]) %>% 
    mutate(min_diff = rsqr - rsqr_max) %>% 
    mutate(cell_type = celltype)
  
}) %>% session$dataWrangler$rbind()

fdata$rsqrIndeps$velmeshev <- fdata$rsqrIndeps$velmeshev %>% arrange(desc(min_diff))

# rosmap

fdata$rsqrIndeps$rosmap <- fdata$coefstats$rosmap%>% 
  dplyr::select(gene_id, cell_type, rsqr_indep) %>% 
  spread(cell_type, rsqr_indep) %>% 
  session$dataWrangler$setColAsRownames("gene_id")

fdata$rsqrIndeps$rosmap <- fdata$celltypes %>% lapply(function(celltype) {
  
  targetExpr <- fdata$rsqrIndeps$rosmap[, celltype]
  othrExprMax <- fdata$rsqrIndeps$rosmap[, colnames(fdata$rsqrIndeps$rosmap) != celltype] %>% apply(1, max) 
  
  tibble(gene_id = rownames(fdata$rsqrIndeps$rosmap)) %>% 
    mutate(rsqr = targetExpr, 
           rsqr_max = othrExprMax[gene_id]) %>% 
    mutate(min_diff = rsqr - rsqr_max) %>% 
    mutate(cell_type = celltype)
  
}) %>% session$dataWrangler$rbind()

fdata$rsqrIndeps$rosmap <- fdata$rsqrIndeps$rosmap %>% arrange(desc(min_diff))

# now work out the expression profiles

fdata$umkrs <- read_rds(paste0(workspace$outputDir, "umkrs.rds"))

fdata$umkrs$minfc %>% group_by(dataset, cell_type) %>% summarize(n_gene = n())

fdata$ctprofiles <- fdata$umkrs$minfc %>% dplyr::select(gene_id, cell_type, expr, dataset) %>%  mutate(expr = log2(expr + 1))


# =================================================
# OUTPUT FIGURE

# velmeshev first

xdata <- list()

xdata$figName <- "figure_07.eps"
xdata$figName2 <- "figure_07_1.eps"

xdata$topGenes <- fdata$rsqrIndeps$velmeshev %>% 
  group_by(gene_id) %>% 
  filter(rsqr == max(rsqr)) %>% 
  filter(rsqr > 0.3) %>% dplyr::select(gene_id, max_ct = cell_type) %>% 
  ungroup()

xdata$main <- fdata$ctprofiles %>% 
  filter(dataset == "velmeshev") %>% 
  inner_join(xdata$topGenes) %>% 
  group_by(max_ct, cell_type) %>% 
  summarize(expr = mean(expr)) %>% 
  ungroup() 

xdata$main <- xdata$main %>% group_by(max_ct) %>% mutate(expr = expr - mean(expr)) %>% ungroup()
xdata$main <- xdata$main %>% mutate(max_ct = workspace$utils$fmtCelltypes(max_ct), 
                                    cell_type = workspace$utils$fmtCelltypes(cell_type))

xdata$max <- xdata$main %>% group_by(max_ct) %>% filter(expr == max(expr)) %>% ungroup()

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = max_ct)) + 
  geom_tile(aes(fill = expr)) + 
  scale_fill_gradient(low = "white", high = "red", na.value = "grey60") + 
  geom_tile(data = xdata$max, aes(fill = expr), colour = "black", linewidth = 1) +
  session$graphingUtils$tiltX(angle = 90) + 
  ggtitle("Cell type specific genes are\nmore resistant to dilution", "Velmeshev") + 
  labs(fill = "Relative expr.")
  
xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 19, height = 19) 

xdata$plot2 <- xdata$topGenes %>% group_by(max_ct) %>% summarize(n_gene = n()) %>% 
  session$graphingUtils$ggplot(aes(x = n_gene, y = max_ct)) + 
  geom_bar(stat = "identity") + 
  scale_x_continuous(trans = "log10") +
  geom_text(aes(label = n_gene), hjust = -0.1, size = 8) +
  xlim(0, 660) +
  xlab("Number of genes") + 
  ggtitle("", "")

xdata$plot2

ggsave(filename = paste0(pdata$figuresDir, xdata$figName2), 
       plot = xdata$plot2, device = "eps", units = "cm", width = 19, height = 19) 

# =================================================


# =================================================
# OUTPUT FIGURE

# ROSMAP first

xdata <- list()

xdata$figName <- "figure_08.eps"
xdata$figName2 <- "figure_08_1.eps"

xdata$topGenes <- fdata$rsqrIndeps$rosmap %>% 
  group_by(gene_id) %>% 
  filter(rsqr == max(rsqr)) %>% 
  filter(rsqr > 0.3) %>% dplyr::select(gene_id, max_ct = cell_type) %>% 
  ungroup()

xdata$main <- fdata$ctprofiles %>% 
  filter(dataset == "rosmap") %>% 
  inner_join(xdata$topGenes) %>% 
  group_by(max_ct, cell_type) %>% 
  summarize(expr = mean(expr)) %>% 
  ungroup() 

xdata$main <- xdata$main %>% group_by(max_ct) %>% mutate(expr = expr - mean(expr)) %>% ungroup()
xdata$main <- xdata$main %>% mutate(max_ct = workspace$utils$fmtCelltypes(max_ct), 
                                    cell_type = workspace$utils$fmtCelltypes(cell_type))

xdata$max <- xdata$main %>% group_by(max_ct) %>% filter(expr == max(expr)) %>% ungroup()

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = max_ct)) + 
  geom_tile(aes(fill = expr)) + 
  scale_fill_gradient(low = "white", high = "red", na.value = "grey60") + 
  geom_tile(data = xdata$max, aes(fill = expr), colour = "black", linewidth = 1) +
  session$graphingUtils$tiltX(angle = 90) + 
  ggtitle("Cell type specific genes are\nmore resistant to dilution", "ROSMAP") + 
  labs(fill = "Relative expr.")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 19, height = 19) 

xdata$plot2 <- xdata$topGenes %>% group_by(max_ct) %>% summarize(n_gene = n()) %>% 
  session$graphingUtils$ggplot(aes(x = n_gene, y = max_ct)) + 
  geom_bar(stat = "identity") + 
  scale_x_continuous(trans = "log10") +
  geom_text(aes(label = n_gene), hjust = -0.1, size = 8) +
  xlim(0, 270) +
  xlab("Number of genes") + 
  ggtitle("", "")

xdata$plot2

ggsave(filename = paste0(pdata$figuresDir, xdata$figName2), 
       plot = xdata$plot2, device = "eps", units = "cm", width = 19, height = 19) 

# =================================================


# +++++++++++++++++++++====================================
# now work in inter cell type synchrony and generate the plots here

fdata <- list(); gc()

fdata$psdbkres <- read_rds(paste0(workspace$outputDir, paste0("bk_exprmats_rosmap-psdbkres.rds"))) # OK BY ACCIDENT I USED A DIFFERENT DATASET FOR TESTING HA
fdata$exprmats$sbj <- fdata$psdbkres$exprmat$sbj
fdata$exprmats$bk <- fdata$psdbkres$exprmat$bkOrig

fdata$samples <- fdata$exprmats$sbj$microglia %>% colnames() %>% sort() # use microglia as it is the limiting factor
fdata$genes <- fdata$exprmats$sbj$microglia %>% rownames() %>% unique() %>% sort()
names(fdata$genes) <- fdata$genes 

fdata$psdbkexprmats <- fdata$exprmats$sbj %>% 
  lapply(function(currExprmat) { currExprmat[fdata$genes, fdata$samples] })

fdata$ctSyncTbl <- do.call(rbind, fdata$genes %>% session$collectionUtils$lapply(function(currGene) {
  coexmat <- fdata$psdbkexprmats %>% sapply(function(currExprmat) {
    currExprmat[currGene, ]
  }) %>% t() %>% workspace$utils$computeCoex() 
  workspace$utils$getCoords(coexmat) %>% 
    mutate(cor_coef = workspace$utils$vectorize(coexmat)) %>% 
    mutate(gene_id = currGene) %>% 
    dplyr::select(gene_id, everything())
}))

fdata$ctSyncSmry <- fdata$ctSyncTbl %>%
  group_by(gene_id) %>% 
  summarize(cor_coef = mean(cor_coef, na.rm = TRUE))

fdata$ctSyncSmry <- fdata$ctSyncSmry %>% 
  left_join(pdata$genes %>% dplyr::select(gene_id = ensembl, gene), by = "gene_id")

fdata$ctSyncSmry <- fdata$ctSyncSmry %>% 
  mutate(cor_perc = rank(cor_coef) / length(cor_coef))

fdata$ctSyncSmry <- fdata$ctSyncSmry %>% dplyr::select(gene_id, gene, cor_coef, cor_perc)

fdata$bins <- fdata$ctSyncSmry %>% mutate(bin = cut(cor_perc, breaks = 10))
fdata$bins %>% group_by(bin) %>% summarize(n = n())

fdata$syncImpact <- fdata$bins$bin %>% unique() %>% session$collectionUtils$lapply(function(currBin) {
  genes <- fdata$bins %>% filter(bin == currBin)
  
  pdata$celltypes %>% lapply(function(celltype) {
    
    exprmatSbj <- fdata$exprmats$sbj[[celltype]]
    exprmatBk <- fdata$exprmats$bk
    
    geneIds <- genes$gene_id %>% intersect(rownames(exprmatSbj)) %>% intersect(rownames(exprmatBk))
    sbj <- exprmatSbj[geneIds,] %>% workspace$utils$computeCoexmat() %>% workspace$utils$coexMatToTbl()
    bk <- exprmatBk[geneIds,] %>% workspace$utils$computeCoexmat() %>% workspace$utils$coexMatToTbl()
    main <- sbj %>% dplyr::select(pair_id, sbj = cor_coef) %>% 
      left_join(bk %>% dplyr::select(pair_id, bk = cor_coef), by = "pair_id")
    tibble(bin = currBin, 
           pearson = cor(main$sbj, main$bk, use = 'pairwise.complete.obs')[1, 1], 
           spearman = cor(main$sbj, main$bk, use = 'pairwise.complete.obs', method = "spearman")[1, 1], 
           n_genes = length(geneIds), 
           cell_type = celltype)
    
  }) %>% session$dataWrangler$rbind()
  
}) %>% session$dataWrangler$rbind()


fdata$output <- list()
fdata$output$psdbkexprmats <- fdata$psdbkexprmats 
fdata$output$ctSyncTbl <- fdata$ctSyncTbl 
fdata$output$ctSyncSmry <- fdata$ctSyncSmry 
fdata$output$syncImpact <- fdata$syncImpact 
fdata$output %>% saveRDS(paste0(workspace$outputDir, "ctsync_smry_rosmap.rds"))


# ++++++++++++++++++++++++++++++++++++++++++++++++
# figures time again... 

fdata <- list()

fdata$ctSyncDat$rosmap <- read_rds(paste0(workspace$outputDir, "ctsync_smry_rosmap.rds"))
fdata$ctSyncDat$velmeshev <- read_rds(paste0(workspace$outputDir, "ctsync_smry_velmeshev.rds"))

fdata$rpgs <- read.table(paste0(workspace$miscDir, "ribosomal_protein_genes.txt"), header = TRUE)$rpg


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "sfigure_01.png"

xdata$main <- read_rds(paste0(workspace$outputDir, "ctsync_smry_velmeshev.rds"))$syncImpact

xdata$main <- xdata$main %>% mutate(cell_type = workspace$utils$fmtCelltypes(cell_type))

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(y = bin, x = pearson)) + geom_bar(stat = "identity") + 
  facet_wrap(~cell_type) + 
  xlab("Correlation of co-expression (Pearson's r)") +
  ylab("Inter-cell type synchrony (normalized correlation)") +
  ggtitle("Inter-cell type synchrony enables preservation of xSubject co-expression\npatterns at the xBulk level", 
          "Velmeshev")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 37, height = 28) 
# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "sfigure_02.png"

xdata$main <- read_rds(paste0(workspace$outputDir, "ctsync_smry_rosmap.rds"))$syncImpact

xdata$main <- xdata$main %>% mutate(cell_type = workspace$utils$fmtCelltypes(cell_type))

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(y = bin, x = pearson)) + geom_bar(stat = "identity") + 
  facet_wrap(~cell_type) + 
  xlab("Correlation of co-expression (Pearson's r)") +
  ylab("Inter-cell type synchrony (normalized correlation)") +
  ggtitle("Inter-cell type synchrony enables preservation of xSubject co-expression\npatterns at the xBulk level", 
          "ROSMAP")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 37, height = 28) 
# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "sfigure_03.png"


xdata$plot <- fdata$ctSyncDat$rosmap$ctSyncSmry %>% mutate(dataset = "ROSMAP") %>% 
  rbind(fdata$ctSyncDat$velmeshev$ctSyncSmry %>% mutate(dataset = "Velmeshev")) %>% 
  filter(gene %in% fdata$rpgs) %>% 
  session$graphingUtils$ggplot(aes(x = cor_perc)) + 
  geom_histogram() + 
  facet_wrap(~dataset, ncol = 1) + 
  ggtitle("Distribution of inter-cell type synchrony across RPGs") +
  xlab("Inter-cell type synchrony (normalized correlation)") + 
  ylab("Count")


xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 30, height = 30) 
# =================================================



















# TODO --- I don't think there is anything useful after this point




# sanity check that the ones that are significant have elevated beta values

fdata$coefstats$velmeshev %>% 
  mutate(significant = gene_id %in% fdata$sigGenes$gene_id) %>% 
  session$graphingUtils$ggplot(aes(x = beta, y = significant)) + 
  geom_boxplot() +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  facet_wrap(~cell_type, ncol = 1)

fdata$coefstats$rosmap %>% 
  mutate(significant = gene_id %in% fdata$sigGenes$rosmap$gene_id) %>% 
  session$graphingUtils$ggplot(aes(x = beta, y = significant)) + 
  geom_boxplot() +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  facet_wrap(~cell_type, ncol = 1)

# no look at the relationship between model r & independent r

# there are 2 types - ones with single cell type contribution, others with multiple 

xdata <- list()

fdata$coefstats$velmeshev %>% filter(qvalue < 0.1) %>% group_by(gene_id) %>% summarize(n_cell_type = n()) %>% arrange(desc(n_cell_type)) %>% 
  session$graphingUtils$ggplot(aes(x = n_cell_type)) + geom_histogram()

xdata$main <- fdata$coefstats$rosmap %>% filter(gene_id %in% fdata$sigGenes$rosmap$gene_id)

xdata$main <- xdata$main %>% group_by(gene_id) %>% filter(rsqr_inde == max(rsqr_inde)) %>% ungroup()

xdata$main %>% session$graphingUtils$ggplot(aes(x = rsqr_inde, y = rsqr_full)) + geom_point(aes(color = cell_type), shape = 1)

xdata$main %>% group_by(cell_type) %>% summarize(n_gene = n()) %>% session$graphingUtils$ggplot(aes(x = cell_type, y = n_gene)) + geom_bar(stat = "identity")

xdata$main <- xdata$main %>% mutate(fc_ovr_max = rsqr_full / rsqr_inde)   

# first look for genes that are "spread thin"
xdata$main %>% arrange(desc(fc_ovr_max)) %>% filter(qvalue < 0.1) %>% arrange(gene_id)
  
xdata$main %>% session$graphingUtils$ggplot(aes(x = fc_ovr_max)) + geom_histogram()

xdata$main %>% filter(rsqr_full > 0.6) %>% arrange(desc(fc_ovr_max)) %>% View()

fdata$coefstats$velmeshev %>% filter(gene_id %in% fdata$clusterFlat$microglia.id_0001) %>% left_join(pdata$genes)
fdata$lmstats$velmeshev %>% filter(gene_id == "ENSG00000008394") %>% left_join(pdata$genes)

# examples where the information is spread out across multiple cell types
# ENSG00000013573 (excitatory & inhibitory) ; ENSG00000182985 (inhibitory & astrocyte) 

# FOR VELMESHEV -- try plotting these scatter plots
# ENSG00000138378 STAT4 (excitatory)
# ENSG00000008394 MGST1 (astrocyte)

# ++ now let's plot some scatter plots for velmeshev

xdata <- list()

xdata$gene <- "ENSG00000014138"

xdata$expr$sbj <- fdata$psdbkres$velmeshev$exprmat$sbj %>% session$collectionUtils$lapplyWithName(function(celltype, exprmat) { 
  exprmat[xdata$gene, ] %>% 
    workspace$utils$stdVec() %>% 
    session$dataWrangler$vectorToTibble() %>% dplyr::select(sample = variable, expr = value) %>% mutate(cell_type = celltype)
}) %>% session$dataWrangler$rbind()

xdata$expr$bk <- fdata$psdbkres$velmeshev$exprmat$bkOrig[xdata$gene, ] %>% 
  workspace$utils$stdVec() %>% 
  session$dataWrangler$vectorToTibble() %>% dplyr::select(sample = variable, expr = value)

xdata$expr$sbj %>% 
  left_join(xdata$expr$bk, by = "sample") %>% 
  session$graphingUtils$ggplot(aes(x = expr.y, y = expr.x)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  facet_wrap(~cell_type)


fdata$coefstats$velmeshev %>% filter(cell_type == "oligodendrocyte") %>% filter(rsqr_boost > 0.1) %>% arrange(desc(rsqr_inde))



# 1. those genes significantly predictable by a particular cell type would be especially highly expressed in that cell type
# 2. those genes with high expression in a particular cell type would also have high co-expression preservation; No you cannot do this because it is heavily influeced by CCV
# 3. those genes with particularly high synchrony would also have high co-expression preservation


fdata$ctprofiles <- read_rds(paste0(workspace$outputDir, "sc_stats_ctprofiles.rds"))

fdata$uexprdat <- read_rds(paste0(workspace$outputDir, "uexprdats_velmeshev.rds"))

fdata$umkrs <- read_rds(paste0(workspace$outputDir, "umkrs.rds"))


xdata <- list()

(xdata$genes <- fdata$coefstats$velmeshev %>% filter(cell_type == "inhibitory") %>% filter(rsqr_inde > 0.1, qvalue < 0.1, beta > 0) %>% arrange(desc(rsqr_inde))) 

y %>% filter(gene_id %in% fdata$clusterFlat$excitatory.id_0010) %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = expr_mean)) + 
  geom_boxplot() + 
  session$graphingUtils$tiltX(angle = 90)

y %>% filter(gene_id %in% fdata$clusterFlat$microglia.id_0006) %>% 
  spread(cell_type, expr_mean) %>% 
  session$dataWrangler$setColAsRownames("gene_id") %>% 
  apply(1, workspace$utils$stdVec) %>% t() %>% 
  session$graphingUtils$heatmap()



fdata$ctprofiles %>% filter(gene_id %in% xdata$genes$gene_id) %>% 
  session$graphingUtils$ggplot(aes(y = cell_type, x = expr)) + 
  geom_boxplot() +
  facet_wrap(~dataset, ncol = 1) +
  scale_x_continuous(trans = "log2")

fdata$ctprofiles %>% filter(gene_id %in% xdata$genes$gene_id) %>% 
  group_by(dataset, cell_type) %>% 
  summarize(expr = mean(expr)) %>% 
  ungroup() %>% 
  spread(cell_type, expr) %>% 
  session$dataWrangler$setColAsRownames("dataset") %>% 
  session$graphingUtils$heatmap(cluster_row = FALSE, cluster_col = FALSE)


fdata$coexmats$sbj <- pdata$celltypes %>% session$collectionUtils$lapply(function(celltype) { read_rds(paste0(workspace$outputDir, paste0("sbj_coexmats_velmeshev_", celltype, ".rds"))) })
fdata$coexmats$bk <- read_rds(paste0(workspace$outputDir, paste0("bk_coexmats_velmeshev.rds")))


fdata$exprmats$sbj <- read_rds(paste0(workspace$outputDir, paste0("sbj_exprmats_velmeshev.rds")))
fdata$exprmats$bk <- read_rds(paste0(workspace$outputDir, paste0("bk_exprmats_velmeshev.rds")))

xdata <- list()

# xdata$genes <- fdata$umkrs$minfc %>% filter(dataset == "velmeshev", cell_type == "microglia") %>% filter(min_fc_rank < 100)

# just compute the co-expressions here
# TODO BELOW THIS IS THE CODE FOR COMPUTING SYNCHRONY & CORRELATION OF CO-EXPRESSIONS

xdata$exprmat$sbj <- fdata$exprmats$sbj$exprmats$excitatory
xdata$exprmat$bk <- fdata$exprmats$bk$exprmat

xdata$bins <- fdata$ctSyncSmry %>% mutate(bin = cut(cor_coef, breaks = 5))
xdata$bins %>% group_by(bin) %>% summarize(n = n())

xdata$syncImpact <- xdata$bins$bin %>% unique() %>% session$collectionUtils$lapply(function(currBin) {
  genes <- xdata$bins %>% filter(bin == currBin)
  
  pdata$celltypes %>% lapply(function(celltype) {
    
    exprmatSbj <- fdata$exprmats$sbj$exprmats[[celltype]]
    exprmatBk <- fdata$exprmats$bk$exprmat
    
    geneIds <- genes$gene_id %>% intersect(rownames(exprmatSbj)) %>% intersect(rownames(exprmatBk))
    sbj <- exprmatSbj[geneIds,] %>% workspace$utils$computeCoexmat() %>% workspace$utils$coexMatToTbl()
    bk <- exprmatBk[geneIds,] %>% workspace$utils$computeCoexmat() %>% workspace$utils$coexMatToTbl()
    main <- sbj %>% dplyr::select(pair_id, sbj = cor_coef) %>% 
      left_join(bk %>% dplyr::select(pair_id, bk = cor_coef), by = "pair_id")
    tibble(bin = currBin, 
           pearson = cor(main$sbj, main$bk, use = 'pairwise.complete.obs')[1, 1], 
           spearman = cor(main$sbj, main$bk, use = 'pairwise.complete.obs', method = "spearman")[1, 1], 
           n_genes = length(geneIds), 
           cell_type = celltype)
    
  }) %>% session$dataWrangler$rbind()
  
}) %>% session$dataWrangler$rbind()

xdata$syncImpact %>% session$graphingUtils$ggplot(aes(y = bin, x = pearson)) + geom_bar(stat = "identity") + 
  geom_text(aes(label = round(pearson, 2)), hjust = -0.5, size = 5) +
  facet_wrap(~cell_type, ncol = 1) + 
  xlim(0, 1) + 
  xlab("Correlation of COEX: xSubject vs. xBulk") +
  ylab("Level of synchrony among cell types")
  





xdata$genes <- x$gene_id %>% intersect(rownames(xdata$exprmat$sbj)) %>% intersect(rownames(xdata$exprmat$bk))

fdata$coextbl$sbj <- xdata$exprmat$sbj[xdata$genes,] %>% workspace$utils$computeCoexmat() %>% workspace$utils$coexMatToTbl()
fdata$coextbl$bk <- xdata$exprmat$bk[xdata$genes,] %>% workspace$utils$computeCoexmat() %>% workspace$utils$coexMatToTbl()

fdata$coextbl$main <- fdata$coextbl$sbj %>% dplyr::select(pair_id, sbj = cor_coef) %>% 
  left_join(fdata$coextbl$bk %>% dplyr::select(pair_id, bk = cor_coef), by = "pair_id")

fdata$coextbl$main <- fdata$coextbl$main %>% filter(bk < 0.9, bk > -0.9)

fdata$coextbl$main %>% session$graphingUtils$ggplot(aes(x = sbj, y = bk)) + geom_point()
cor.test(fdata$coextbl$main$sbj, fdata$coextbl$main$bk)


pdata$genes %>% filter(gene_id %in% xdata$genes) %>% View()


# next, let's try the synchrony hypothesis


# fdata <- list()

fdata$psdbkres <- read_rds(paste0(workspace$outputDir, paste0("bk_exprmats_velmeshev-psdbkres.rds"))) # OK BY ACCIDENT I USED A DIFFERENT DATASET FOR TESTING HA
fdata$exprmats <- fdata$psdbkres$exprmat$sbj

fdata$samples <- fdata$exprmats$microglia %>% colnames() %>% sort() # use microglia as it is the limiting factor
fdata$genes <- fdata$exprmats$astrocyte %>% rownames() %>% unique() %>% sort()
names(fdata$genes) <- fdata$genes 

fdata$psdbkexprmats <- fdata$exprmats %>% 
  lapply(function(currExprmat) { currExprmat[fdata$genes, fdata$samples] })

fdata$ctSyncTbl <- do.call(rbind, fdata$genes %>% session$collectionUtils$lapply(function(currGene) {
  coexmat <- fdata$psdbkexprmats %>% sapply(function(currExprmat) {
    currExprmat[currGene, ]
  }) %>% t() %>% workspace$utils$computeCoex() 
  workspace$utils$getCoords(coexmat) %>% 
    mutate(cor_coef = workspace$utils$vectorize(coexmat)) %>% 
    mutate(gene_id = currGene) %>% 
    dplyr::select(gene_id, everything())
}))

fdata$ctSyncSmry <- fdata$ctSyncTbl %>%
  group_by(gene_id) %>% 
  summarize(cor_coef = mean(cor_coef, na.rm = TRUE)) # ENSG00000004848

fdata$ctSyncSmry <- fdata$ctSyncSmry %>% 
  left_join(pdata$genes %>% dplyr::select(gene_id = ensembl, gene), by = "gene_id")

fdata$ctSyncSmry <- fdata$ctSyncSmry %>% 
  mutate(cor_perc = rank(cor_coef) / length(cor_coef))

fdata$ctSyncSmry <- fdata$ctSyncSmry %>% dplyr::select(gene_id, gene, cor_coef, cor_perc)

fdata$ctSyncSmry %>% filter(cor_perc > 0, cor_perc < 0.1) -> x
fdata$ctSyncSmry %>% filter(cor_perc > 0.99, cor_perc < 1) -> x


yy <- fdata$clusterFlat[x$cluster] %>% session$collectionUtils$lapplyWithName(function(clusterId, cluster) {
  fdata$ctSyncSmry %>% filter(gene_id %in% cluster) %>% mutate(cluster = clusterId)
}) %>% session$dataWrangler$rbind()



yy %>% 
  group_by(cluster) %>% summarize(cor_coef = mean(cor_coef)) %>% 
  session$graphingUtils$ggplot(aes(y = cluster, x = cor_coef)) + 
  geom_bar(stat = "identity") + 
  ylim(rev(x$cluster))

fdata$diffs %>% filter(cluster %in% xdata$ylim) %>% 
  filter(procedure == "mgp", dataset == "velmeshev") %>% dplyr::select(cluster, diff) %>% 
  left_join(yy %>% 
              group_by(cluster) %>% summarize(cor_coef = mean(cor_coef)), by = "cluster") -> zz

zz %>% session$graphingUtils$ggplot(aes(x = cor_coef, y = diff)) + geom_point() + geom_smooth(method = "lm", se = FALSE)+ 
  ggtitle("Correlation between cell type synchrony and\nCCV (un)importance; Spearman = 0.42")

zz$diff %>% cor.test(zz$cor_coef, method = "spearman") 





ydata <- list()

ydata <- fdata$uexprdat$sbj$exprmats %>% lapply(function(exprmat) { 
  exprmat %>% apply(1, mean) -> x
  # (x - mean(x)) / sd(x)
  rank(-x) -> x
  x <- max(x) - x
})

ydata <- ydata %>% 
  session$collectionUtils$lapplyWithName(function(celltype, vals) { vals %>% session$dataWrangler$vectorToTibble() %>% dplyr::select(gene_id = variable, expr_mean = value) %>% mutate(cell_type = celltype) }) %>% session$dataWrangler$rbind() 

ydata -> y




