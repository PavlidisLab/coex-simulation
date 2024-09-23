pdata <- list()

pdata$figuresDir <- paste0(workspace$workspaceDir, "041_ANALYSIS_02_xsbj_coex_FIGURES/")

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

# ++++++++++++++
# figures
# Now this is where I'm going to look at reproducibility of the xSubject networks

fdata <- list(); gc()

fdata$coexmatFiles <- pdata$files$coexmats %>% filter(level == "sbj")

fdata$celltypes <- fdata$coexmatFiles$cell_type %>% unique() %>% sort()
names(fdata$celltypes) <- fdata$celltypes

fdata$celltypes %>% session$collectionUtils$lapply(function(currCelltype) {
  
  print(paste0("working on cell type: ", currCelltype, "..."))
  
  lnktbls <- fdata$coexmatFiles %>% filter(cell_type == currCelltype)
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
  lnktbls %>% saveRDS(paste0(workspace$outputDir, paste0("sbj_lnktbls_", currCelltype, ".rds")))
  
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

fdata$lnktblFiles <- pdata$files$lnktbls %>% filter(level == "sbj")

fdata$celltypes <- fdata$lnktblFiles$cell_type %>% unique()
names(fdata$celltypes) <- fdata$celltypes

fdata$main <- fdata$celltypes %>% session$collectionUtils$lapply(function(currCelltype) {
  
  print(paste0("working on cell type: ", currCelltype, "..."))
  
  lnkstbls <- (fdata$lnktblFiles %>% filter(cell_type == currCelltype))$outfile
  lnkstbls <- readRDS(paste0(workspace$outputDir, lnkstbls))
  
  datasets <- lnkstbls %>% names()
  names(datasets) <- datasets
  
  print(paste0("preparing lnktbls..."))
  
  lnks_T <- lnkstbls %>% lapply(function(currLnkTbl) { currLnkTbl %>% filter(lnk == 1) })
  lnks_F <- lnkstbls %>% lapply(function(currLnkTbl) { currLnkTbl %>% filter(lnk == 0) })
  
  print(paste0("running fisher's exact tests in parallel..."))
  
  datasets %>% lapply(function(currDat_T) {
    
    datasets_S <- datasets[datasets != currDat_T]
    
    datasets_S %>% mclapply(function(currDat_S) {
      
      session$evaluationUtils$fisher(predicted = lnks_T[[currDat_S]]$pair, 
                                     notPredicted = lnks_F[[currDat_S]]$pair, 
                                     trueSet = lnks_T[[currDat_T]]$pair, 
                                     alternative = "greater")
      
    }, mc.cores = 10)
    
  })
  
})

# COMMIT
fdata$main %>% saveRDS(paste0(workspace$outputDir, paste0("sbj_repro_fisher.rds")))


# save the tbl version of fishers
fdata <- list()

fdata$fishers <- read_rds(paste0(workspace$outputDir, paste0("sbj_repro_fisher.rds"))) 

fdata$main <- fdata$fishers %>% 
  session$collectionUtils$lapplyWithName(function(currCelltype, currCTComps) {
    
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
fdata$output %>% saveRDS(paste0(workspace$outputDir, paste0("sbj_repro_fisher.rds")))
# ++++++++++++++++++++++++


# +++++++++++++++++++++
# now, compute preservation!

fdata <- list(); gc()

fdata$lnktblFiles <- pdata$files$lnktbls

fdata$celltypes <- fdata$lnktblFiles$cell_type %>% unique()
names(fdata$celltypes) <- fdata$celltypes

fdata$main <- fdata$celltypes %>% session$collectionUtils$lapply(function(currCelltype) {
  
  print(paste0("working on cell type: ", currCelltype, "..."))
  
  lnkstbls_sc <- (fdata$lnktblFiles %>% filter(cell_type == currCelltype, level == "sc"))$outfile
  lnkstbls_sbj <- (fdata$lnktblFiles %>% filter(cell_type == currCelltype, level == "sbj"))$outfile
  
  lnkstbls_sc <- readRDS(paste0(workspace$outputDir, lnkstbls_sc))
  lnkstbls_sbj <- readRDS(paste0(workspace$outputDir, lnkstbls_sbj))
  
  datasets <- names(lnkstbls_sc) %>% intersect(names(lnkstbls_sbj)) %>% unique() %>% sort()
  names(datasets) <- datasets
  
  print(paste0("running fisher's exact tests in parallel..."))
  
  datasets %>% mclapply(function(currDat) {
    
    sc <- lnkstbls_sc[[currDat]]
    sc_T <- sc %>% filter(lnk == 1)
    
    sbj <- lnkstbls_sbj[[currDat]]
    sbj_T <- sbj %>% filter(lnk == 1)
    sbj_F <- sbj %>% filter(lnk == 0)
    
    session$evaluationUtils$fisher(predicted = sbj_T$pair, 
                                   notPredicted = sbj_F$pair, 
                                   trueSet = sc_T$pair, 
                                   alternative = "greater")
    
  }, mc.cores = length(datasets))

})

# COMMIT
fdata$main %>% saveRDS(paste0(workspace$outputDir, paste0("sc_sbj_preservation_fisher.rds")))

# save the tbl version of fishers
fdata <- list()

fdata$fishers <- read_rds(paste0(workspace$outputDir, paste0("sc_sbj_preservation_fisher.rds"))) 

fdata$main <- fdata$fishers %>% 
  session$collectionUtils$lapplyWithName(function(currCelltype, currCTComps) {
    
    currCTComps %>% session$collectionUtils$lapplyWithName(function(currDatTruth, currComparison) {
      
      currComparison$stats %>% as.data.frame() %>% t() %>% as_tibble() %>% 
        mutate(pvalue = currComparison$test$p.value) %>% 
        mutate(dat_truth = currDatTruth, cell_type = currCelltype) %>% 
        dplyr::select(cell_type, dat_truth, everything())
      
    }) %>% session$dataWrangler$rbind()
    
  }) %>% session$dataWrangler$rbind()

# ++++++++++++++++++++++++
# COMMIT
fdata$output$fishers <- fdata$fishers
fdata$output$tbl <- fdata$main
fdata$output %>% saveRDS(paste0(workspace$outputDir, paste0("sc_sbj_preservation_fisher.rds")))
# ++++++++++++++++++++++++


# ++++++++++++++
# Now, examine the overall topolog (node degrees, etc)
# and focus on those edges that are reproducible across datasets

fdata <- list(); gc()

fdata$coexmatsFiles <- pdata$files$coexmats %>% filter(level == "sbj")

fdata$celltypes <- fdata$coexmatsFiles$cell_type %>% unique() %>% sort()
names(fdata$celltypes) <- fdata$celltypes

fdata$celltypes %>% session$collectionUtils$lapply(function(currCelltype) {
  
  print(paste0("working on cell type: ", currCelltype, "..."))
  
  lnkmats <- fdata$coexmatsFiles %>% filter(cell_type == currCelltype)
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
  lnkmats %>% saveRDS(paste0(workspace$outputDir, paste0("sbj_lnkmats_", currCelltype, ".rds")))
  
})

pdata$files$lnkmats <- tibble(outfile = list.files(workspace$outputDir)) %>% 
  filter(grepl("lnkmats", outfile)) %>% 
  mutate(outfile_name = gsub(".rds", "", outfile)) %>% 
  mutate(strs = strsplit(outfile_name, "_")) %>% 
  mutate(level = strs %>% sapply(function(currStr) { currStr[1] })) %>% 
  mutate(type = strs %>% sapply(function(currStr) { currStr[2] })) %>% 
  mutate(cell_type = strs %>% sapply(function(currStr) { currStr[3] })) %>% 
  dplyr::select(outfile, level, type, cell_type)


# +++++++++++++++++++++
# compute and examine node degrees 

fdata <- list(); gc()

fdata$lnkmatsFiles <- pdata$files$lnkmats %>% filter(level == "sbj")

fdata$celltypes <- fdata$lnkmatsFiles$cell_type %>% unique() %>% sort()

fdata$files <- fdata$lnkmatsFiles$outfile
names(fdata$files) <- fdata$lnkmatsFiles$cell_type

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
  geom_text(aes(label = n_gene), size = 6, vjust = -0.5) +
  facet_wrap(~cell_type * dataset) + 
  ylim(0, 14000) +
  xlab("Number of genes in network") 

# density distribution of node degrees
fdata$degrees %>% 
  session$graphingUtils$ggplot(aes(x = node_degree_log10)) + 
  geom_density(aes(color = dataset)) + 
  facet_wrap(~cell_type)


# cumulative distribution of node degrees
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
# finally; let's make some figures for ora reproducibility

fdata <- list(); gc()

fdata$sc <- read_rds(paste0(workspace$outputDir, paste0("sc_repro_fisher.rds")))$tbl

fdata$sbj <- read_rds(paste0(workspace$outputDir, paste0("sbj_repro_fisher.rds")))$tbl

fdata$main <- fdata$sc %>% mutate(level = "sc") %>% rbind(fdata$sbj %>% mutate(level = "sbj"))

fdata$main %>% dplyr::select(cell_type, dat_truth, dat_sample, level, or) %>% 
  spread(level, or) %>% 
  na.omit() %>% 
  group_by(cell_type, dat_truth) %>% summarize(sbj = mean(sbj), sc = mean(sc)) %>% mutate(diff = sbj - sc) %>% summarize(diff = mean(diff))

# try showing xCell reproducibility for just one dataset? velmeshev?
fdata$main %>% filter(dat_truth == "velmeshev", level == "sbj") %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = or)) + 
  geom_boxplot(color = "grey50") +
  geom_point(alpha = 0.4, size = 5) +
  scale_y_continuous(trans = "log2") + 
  xlab("Cell type") + 
  ylab("Odds ratio") +
  ggtitle("Reproducibility of xSubject co-expression signals\nin the Velmeshev dataset") +
  session$graphingUtils$tiltX(angle = 90)


fdata$main %>% 
  filter(level == "sbj") %>% 
  session$graphingUtils$ggplot(aes(y = dat_truth, x = or)) + 
  geom_boxplot(color = "grey50") +
  geom_point(alpha = 0.4, size = 5) +
  facet_wrap(~cell_type, ncol = 1) + 
  scale_x_continuous(trans = "log2") + 
  geom_vline(xintercept = 1, linetype = "dashed") + 
  xlab("odds ratio") + 
  ggtitle("Reproducibility of xSubject co-expression signals")



# let's show the means per dataset for MSL presentation
xdata <- list()
xdata$main <- fdata$main %>% mutate(cell_type = workspace$utils$fmtCelltypes(cell_type), 
                                    dat_truth = workspace$utils$fmtDataset(dat_truth))
xdata$means <- xdata$main %>% 
  group_by(cell_type, dat_truth, level) %>% summarize(or = mean(or)) %>% ungroup()
xdata$means %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = or)) + 
  geom_boxplot(aes(color = level), position = position_dodge(width = 0)) +
  geom_point(aes(color = level, group = level), size = 2, position = position_jitter(width = 0.1)) +
  scale_y_continuous(trans = "log2") + 
  labs(x = "Cell type", 
       y = "Recovery of edges by other datasets\n(Odds ratio)", 
       color = "Level") +
  ggtitle("Reproducibility of co-expression edges") +
  session$graphingUtils$tiltX(angle = 90)  + 
  scale_color_brewer(palette = "Dark2")


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_01.eps"

xdata$main <- fdata$sbj %>% mutate(cell_type = workspace$utils$fmtCelltypes(cell_type), 
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
  ggtitle("Reproducibility of xSubject co-expression edges") +
  scale_color_brewer(palette = "Set1") +
  labs(color = "Dataset") +
  scale_y_discrete(limits = rev) +
  theme(axis.text = element_text(size = 13))

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 20, height = 25) 

# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_02.eps"
xdata$figName2 <- "figure_02a.eps"


xdata$main <- fdata$main %>% 
  dplyr::select(cell_type, dat_truth, dat_sample, level, or) %>% 
  spread(level, or) %>% 
  group_by(cell_type, dat_truth) %>% 
  summarize(sbj = mean(sbj), sc = mean(sc)) %>% 
  ungroup()

xdata$main <- xdata$main %>% 
  mutate(cell_type = workspace$utils$fmtCelltypes(cell_type), dataset = workspace$utils$fmtDataset(dat_truth)) %>% 
  dplyr::select(-dat_truth)

# need to compute the paired wilcoxon test
(xdata$pvalues$all <- wilcox.test(xdata$main$sbj, xdata$main$sc, paired = TRUE)$p.value)

xdata$celltypes <- xdata$main$cell_type %>% unique() %>% sort() %>% session$dataWrangler$attachNames()
xdata$pvalues$celltypes <- xdata$celltypes %>% sapply(function(celltype) { 
  dat <- xdata$main %>% filter(cell_type == celltype)
  wilcox.test(dat$sbj, dat$sc, paired = TRUE)$p.value
}) %>% session$dataWrangler$vectorToTibble() %>% 
  dplyr::select(cell_type = variable, pvalue = value) 

xdata$plot2 <- xdata$pvalues$celltypes %>% 
  session$graphingUtils$ggplot(aes(y = cell_type, x = -log10(pvalue))) + 
  geom_bar(stat = "identity") + 
  geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
  xlab("-log10(P-Value)") 
xdata$plot2 

ggsave(filename = paste0(pdata$figuresDir, xdata$figName2), 
       plot = xdata$plot, device = "eps", units = "cm", width = 10, height = 6) 
  
xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = sc, y = sbj)) + 
  geom_point(aes(color = dataset, shape = cell_type), size = 4) + 
  geom_abline(slope = 1, linetype = "dashed") + 
  scale_x_continuous(trans = "log2", limits = c(1, 150)) + 
  scale_y_continuous(trans = "log2", limits = c(1, 150)) +
  xlab("xCell co-expression reproducibility\n(Odds ratio - log2 scale)") + 
  ylab("xSubject co-expression reproducibility\n(Odds ratio - log2 scale)") +   
  ggtitle("Reproducibility of xCell co-expression\nvs. reproducibility of xSubject co-expression") +
  scale_color_brewer(palette = "Set1") 

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 22, height = 15) 

# =================================================


# +++++++++++++++++++++++ 
# compute the 3 numbers, per dataset, cell type: avg sbj_repro, sc_repro, preservation

fdata <- list(); gc()

fdata$sbj <- read_rds(paste0(workspace$outputDir, paste0("sbj_repro_fisher.rds")))$tbl

fdata$sc <- read_rds(paste0(workspace$outputDir, paste0("sc_repro_fisher.rds")))$tbl

fdata$pres <- read_rds(paste0(workspace$outputDir, paste0("sc_sbj_preservation_fisher.rds")))$tbl

fdata$sbjMean <- fdata$sbj %>% 
  dplyr::select(cell_type, dataset = dat_truth, or) %>% 
  group_by(cell_type, dataset) %>% 
  summarize(or = mean(or)) %>% ungroup() %>% 
  mutate(comparison = "sbj_repro")

fdata$scMean <- fdata$sc %>% 
  dplyr::select(cell_type, dataset = dat_truth, or) %>% 
  group_by(cell_type, dataset) %>% 
  summarize(or = mean(or)) %>% ungroup() %>% 
  mutate(comparison = "sc_repro")


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_03.eps"
xdata$figName2 <- "figure_03a.eps"

xdata$main <- fdata$sbjMean %>% dplyr::select(cell_type, dataset, sbj = or) %>% 
  inner_join(fdata$scMean %>% dplyr::select(cell_type, dataset, sc = or), by = c("cell_type", "dataset")) %>% 
  inner_join(fdata$pres %>% dplyr::select(cell_type, dataset = dat_truth, pres = or), by = c("cell_type", "dataset"))

xdata$main <- xdata$main %>% 
  mutate(cell_type = workspace$utils$fmtCelltypes(cell_type), dataset = workspace$utils$fmtDataset(dataset)) 


# need to compute the paired wilcoxon test
xdata$pvalues$all <- wilcox.test(xdata$main$pres, xdata$main$sc, paired = TRUE)$p.value

xdata$celltypes <- xdata$main$cell_type %>% unique() %>% sort() %>% session$dataWrangler$attachNames()
xdata$pvalues$celltypes <- xdata$celltypes %>% sapply(function(celltype) { 
  dat <- xdata$main %>% filter(cell_type == celltype)
  wilcox.test(dat$pres, dat$sc, paired = TRUE)$p.value
}) %>% session$dataWrangler$vectorToTibble() %>% 
  dplyr::select(cell_type = variable, pvalue = value) 

xdata$plot2 <- xdata$pvalues$celltypes %>% 
  session$graphingUtils$ggplot(aes(y = cell_type, x = -log10(pvalue))) + 
  geom_bar(stat = "identity") + 
  geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
  xlab("-log10(P-Value)") 
xdata$plot2 

ggsave(filename = paste0(pdata$figuresDir, xdata$figName2), 
       plot = xdata$plot, device = "eps", units = "cm", width = 10, height = 6) 

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = sc, y = pres)) + 
  geom_point(aes(color = dataset, shape = cell_type), size = 4) + 
  geom_abline(slope = 1, linetype = "dashed") + 
  scale_x_continuous(trans = "log2", limits = c(0.5, 150)) + 
  scale_y_continuous(trans = "log2", limits = c(0.5, 150)) +
  xlab("xCell co-expression reproducibility\n(Odds ratio - log2 scale)") + 
  ylab("Preservation of xCell co-expression\nin xSubject networks\n(Odds ratio - log2 scale)") +
  ggtitle("Cross dataset reproducibility of xCell co-expression\nvs. within dataset preservation ") +
  scale_color_brewer(palette = "Set1") 

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 22, height = 15) 

# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_04.eps"

xdata$main <- fdata$pres %>% 
  dplyr::select(cell_type, dataset = dat_truth, or) %>% 
  mutate(comparison = "sc_sbj_pres") %>% 
  rbind(fdata$sbjMean) %>% 
  rbind(fdata$scMean) 

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = or, color = comparison)) + 
  geom_point(aes(x = or), y = -0.01, shape = 1) +
  geom_density(linewidth = 1) + 
  scale_color_brewer(palette = "Dark2") + 
  ggtitle("Distribution of network similarities between\nthe xCell and xSubject levels") + 
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 14))
  

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 20, height = 15) 

# =================================================


# +++++++++++++++++++++
# also how many edges are discovered in multiple datasets, etc -- AND NOT necessary all between the same datasets?

fdata <- list(); gc()

fdata$lnktblFiles <- pdata$files$lnktbls %>% filter(level == "sbj")

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
fdata$lnkSmrys %>% saveRDS(paste0(workspace$outputDir, "sbj_lnk_smrys.rds"))
# =============

# again single vs. multiple. vs all

fdata <- list(); gc()

fdata$lnkSmrys <- read_rds(paste0(workspace$outputDir, "sbj_lnk_smrys.rds"))

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


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_05.eps"

xdata$main <- fdata$main %>% mutate(cell_type = workspace$utils$fmtCelltypes(cell_type))

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = n_pair)) + 
  geom_bar(aes(fill = factor(class, levels = c("Single", "Multiple", "All"))), stat = "identity") + 
  session$graphingUtils$tiltX(angle = 90) + 
  scale_fill_manual(values = rev(c("grey10", "grey40", "grey80"))) + 
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

xdata$figName <- "figure_06.eps"

xdata$main <- fdata$main %>% mutate(cell_type = workspace$utils$fmtCelltypes(cell_type))

xdata$plot <- xdata$main %>% 
  filter(class == "All") %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = n_pair)) + 
  geom_bar(color = "grey10", stat = "identity") + 
  geom_text(aes(label = n_pair), vjust = -1, size = 6) +
  session$graphingUtils$tiltX(angle = 90) + 
  scale_fill_manual(values = rev(c("grey10", "grey40", "grey80"))) + 
  scale_y_continuous(labels = label_comma(), limits = c(0, 400)) +
  xlab("") + 
  ylab("") 

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 18, height = 12) 

# =================================================
























