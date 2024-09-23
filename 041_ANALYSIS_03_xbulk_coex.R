pdata <- list(); gc()

pdata$figuresDir <- paste0(workspace$workspaceDir, "041_ANALYSIS_03_xbulk_coex_FIGURES/")

pdata$genes <- readRDS(paste0(workspace$outputDir, "genes_metadata.rds")) # genes

pdata$celltypes <- c("excitatory", "inhibitory", "opc", "oligodendrocyte", "astrocyte", "microglia")
names(pdata$celltypes) <- pdata$celltypes

fdata$mkrs <- readRDS(paste0(workspace$outputDir, "umkrs.rds"))

pdata$ctprofiles <- read_rds(paste0(workspace$outputDir, "sc_stats_ctprofiles.rds"))

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

pdata$ctp <- read_rds(paste0(workspace$outputDir, "sc_stats_ctprofiles.rds"))

pdata$mkrs <- read_rds(paste0(workspace$outputDir, "umkrs.rds"))


# +++++++++++++++++++
# quickly examine a few metadata properties 

fdata <- list(); gc()

fdata$exprmats <- list(rosmap = "bk_exprmats_rosmap.rds", velmeshev = "bk_exprmats_velmeshev.rds")

fdata$exprmats <- fdata$exprmats %>% lapply(function(file) { read_rds(paste0(workspace$outputDir, file)) })

fdata$exprmats$rosmap$samples
fdata$exprmats$velmeshev$samples

# +++++++++++++++++++++++++++
# now, compare the genes detected in each cell type and dataset

fdata <- list(); gc()

fdata$datasets <- pdata$files$exprmats %>% filter(level == "bk", dataset %in% c("rosmap", "velmeshev"))
fdata$datasets <- fdata$datasets$outfile %>% session$dataWrangler$attachNames(fdata$datasets$dataset)
fdata$datasets <- fdata$datasets %>% lapply(function(currDataset) { readRDS(paste0(workspace$outputDir, currDataset))$exprmat })

fdata$ctp <- fdata$datasets %>% session$collectionUtils$lapplyWithName(function(datasetName, exprmat) {
  ctp <- exprmat %>%  rowMeans()
  ctp %>% 
    session$dataWrangler$vectorToTibble() %>% 
    dplyr::select(gene_id = variable, expr = value) %>% 
    mutate(level = "bk", dataset = datasetName)
}) %>% session$dataWrangler$rbind()

# =================================  
# COMMIT ==========================
# fdata$ctp %>% saveRDS(paste0(workspace$outputDir, "bk_stats_ctprofiles.rds"))
# =================================

pdata$ctprofiles <- read_rds(paste0(workspace$outputDir, "bk_stats_ctprofiles.rds"))

# +++++++++++++++++++++++++++
# now, compare the gene expressions

# first in terms of detected genes overlap 
# and then in terms of expression levels

fdata <- list()

fdata$ctp$bk <- read_rds(paste0(workspace$outputDir, "bk_stats_ctprofiles.rds"))
fdata$ctp$sc <- read_rds(paste0(workspace$outputDir, "sc_stats_ctprofiles.rds")) %>% filter(dataset %in% fdata$ctp$bk$dataset)

fdata$genesFlat$bk <- fdata$ctp$bk$dataset %>% unique() %>% session$dataWrangler$attachNames() %>%
  session$collectionUtils$lapply(function(datasetName) { fdata$ctp$bk %>% filter(dataset == datasetName) %>% session$dataWrangler$extractColumn("gene_id") %>% unique() })

fdata$genesMat$bk <- fdata$ctp$bk %>% dplyr::select(gene_id, dataset) %>% unique() %>% 
  mutate(in_set = 1) %>% 
  spread(dataset, in_set) %>% 
  session$dataWrangler$fillNa(colNames = unique(fdata$ctp$bk$dataset), value = 0) %>% 
  session$dataWrangler$setColAsRownames("gene_id")


# 1. compare between the 2 bulk tissue datasets

# Velmeshev is effectively a superset of Rosmap; which is reassuring 

fdata$genesMat$bk %>% UpSetR::upset(nsets = 7, text.scale = 2, keep.order = TRUE, nintersects = NA) 

# let's take a look at the correlation in expression

# the two datasets highly agree in terms of co-expression? levels for the intersecting genes

xdata <- list()

xdata$genes <- fdata$genesFlat$bk$rosmap %>% intersect(fdata$genesFlat$bk$velmeshev)

xdata$main <- fdata$ctp$bk %>% filter(gene_id %in% xdata$genes) %>% dplyr::select(gene_id, dataset, expr) %>% spread(dataset, expr)

xdata$pearson <- cor(xdata$main$rosmap, xdata$main$velmeshev)[1, 1]

xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = rosmap, y = velmeshev)) + 
  geom_point(shape = 1) +
  geom_smooth(method = "lm", se = FALSE) + 
  ggtitle("", "Correlation in expression between\nVelmeshev.Bulk and Rosmap.Bulk")

# 2. compare each dataset its single cell counterpart; per cell type

# first let's look at per cell type, how many genes overlap 

xdata <- list()

xdata$dataset <- "rosmap"

xdata$main <- fdata$ctp$sc %>% filter(dataset == xdata$dataset) %>% dplyr::select(gene_id, cell_type, expr_sc = expr) %>% 
  left_join(fdata$ctp$bk %>% filter(dataset == xdata$dataset) %>% dplyr::select(gene_id, expr_bk = expr))

xdata$main <- xdata$main %>% mutate(detected_bk = !is.na(expr_bk))

xdata$smry <- xdata$main %>% group_by(cell_type, detected_bk) %>% summarize(n_gene = n()) %>% ungroup()

(fdata$ctp$sc %>% filter(dataset == xdata$dataset))$gene_id %>% unique() %>% length() # total number of genes at the SC

(fdata$ctp$sc %>% filter(dataset == xdata$dataset))$gene_id %>% unique() %>% 
  intersect((fdata$ctp$bk %>% filter(dataset == xdata$dataset))$gene_id) %>% length() # total number of SC genes also detected at the BK

xdata$main %>% filter(detected_bk) %>% 
  group_by(cell_type) %>% summarize(cor_coef = cor(expr_sc, expr_bk)[1, 1])

xdata$main %>% filter(detected_bk) %>% 
  session$graphingUtils$ggplot(aes(x = expr_bk, y = expr_sc)) +
  geom_point(shape = 1) +
  scale_y_continuous(trans = "log2") +
  geom_smooth(method = "lm", se = FALSE) + 
  facet_wrap(~cell_type) +
  ggtitle("", paste0("Correlation in mean expression between SC and bulk data\n", xdata$dataset))


# 3. out of curiosity, let's see how the genes overlap among the different cell types, within each dataset? 

fdata$ctp$sc %>% filter(dataset == "rosmap") %>% dplyr::select(gene_id, cell_type) %>% 
  mutate(in_set = 1) %>% 
  spread(cell_type, in_set) %>% 
  session$dataWrangler$fillNa(colNames = unique(fdata$ctp$sc$cell_type), value = 0) %>% 
  session$dataWrangler$setColAsRownames("gene_id") %>% 
  UpSetR::upset(nsets = 6, text.scale = 2, keep.order = TRUE, nintersects = NA) 
  
fdata$ctp$sc %>% filter(dataset == "rosmap") %>% dplyr::select(gene_id, cell_type, expr) %>% 
  spread(cell_type, expr) %>% 
  session$dataWrangler$setColAsRownames("gene_id") %>% 
  na.omit() %>% 
  GGally::ggpairs()



# ++++++++++++++
# figures
# Now, examine the overall topolog (node degrees, etc)
# and focus on those edges that are reproducible across datasets

fdata <- list(); gc()

fdata$files <- list(rosmap = "bk_coexmats_rosmap.rds", velmeshev = "bk_coexmats_velmeshev.rds")

fdata$lnkmats <- fdata$files %>% mclapply(function(currFile) {
    
  lnkmat <- read_rds(paste0(workspace$outputDir, currFile))
  lnkmat[which(lnkmat >= 0.99)] <- 1
  lnkmat[which(lnkmat < 0.99)] <- 0
  
  return(lnkmat)
    
}, mc.cores = length(fdata$files))

fdata$lnkmats %>% saveRDS(paste0(workspace$outputDir, paste0("bk_lnkmats.rds")))
  
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

fdata$lnkmats <- read_rds(paste0(workspace$outputDir, "bk_lnkmats.rds"))

fdata$degrees <- fdata$lnkmats %>% session$collectionUtils$lapplyWithName(function(currDataset, currlnkmat) {
  currlnkmat %>% workspace$utils$computeNodeDegrees() %>% mutate(dataset = currDataset)
}) %>% session$dataWrangler$rbind()

# add in log transformed column

fdata$degrees <- fdata$degrees %>% mutate(node_degree_log10 = log10(node_degree + 1))

fdata$degrees %>% group_by(dataset) %>% summarize(node_degree_max = max(node_degree)) %>% arrange(desc(node_degree_max)) %>% ungroup()

# in vs. out of network (having least one edge)
fdata$degrees %>% 
  mutate(in_network = (node_degree > 0)) %>% 
  group_by(dataset, in_network) %>% summarize(n_gene = n()) %>% ungroup() %>% 
  session$graphingUtils$ggplot(aes(x = in_network, y = n_gene)) + 
  geom_bar(stat = "identity") + 
  geom_text(aes(label = n_gene), size = 6, vjust = -0.5) +
  facet_wrap(~dataset) + 
  xlab("Number of genes in network") 


# cumulative distribution of node degrees
fdata$degrees %>% 
  session$graphingUtils$ggplot(aes(x = node_degree_log10, color = dataset)) + 
  stat_ecdf() + 
  geom_vline(xintercept = 2, linetype = "dashed")

# correlation of node degrees
fdata$degrees %>% 
  dplyr::select(gene_id, dataset, node_degree_log10) %>% unique() %>% 
  spread(dataset, node_degree_log10) %>% 
  na.omit() %>% 
  session$dataWrangler$setColAsRownames("gene_id") %>% 
  GGally::ggpairs()


# ++++++++++++++++++++++++++++++=
# compute how well each bulk tissue co-expression network recovers the edges at the subject & sc level, in each of velmeshev and rosmap
# for each dataset: rosmap & velmeshev
# for each cell type: the 6 cell types
# how well does it recover sc links? how well does it recover sbj links? 

fdata <- list(); gc()

# compute lnktbls for the bulk coexmats

fdata$coexmatFiles <- pdata$files$coexmats %>% filter(level == "bk")

fdata$datasets <- fdata$coexmatFiles$outfile
names(fdata$datasets) <- fdata$coexmatFiles$dataset

# fdata$datasets <- fdata$datasets[c("rosmap-ihc", "rosmap-ihcres")] # filter for only target dataset here

1:length(fdata$datasets) %>% mclapply(function(i) {
  
  currFile <- fdata$datasets[i]
  currDataset <- names(fdata$datasets[i])
  
  lnkmat <- read_rds(paste0(workspace$outputDir, currFile))
  lnkmat[which(lnkmat >= 0.99)] <- 1
  lnkmat[which(lnkmat < 0.99)] <- 0
  
  lnktbl <- tibble(pair = lnkmat %>% workspace$utils$getPairIds(), lnk = lnkmat %>% workspace$utils$vectorize())
  
  # ++++++++++++++++++++++++
  # COMMIT
  lnktbl %>% saveRDS(paste0(workspace$outputDir, paste0("bk_lnktbls_", currDataset, ".rds")))
  
}, mc.cores = length(fdata$datasets))

# ++++++++++++++++++++++++++++++=
# compute the level of reproducibility between velmeshev and rosmap at the bulk tissue level

fdata <- list(); gc()

fdata$lnktblFiles <- pdata$files$lnktbls %>% filter(level == "bk") %>% dplyr::select(everything(), dataset = cell_type) %>% 
  filter(dataset %in% c("velmeshev", "rosmap", "velmeshev-mgpres", "rosmap-mgppres"))

fdata$lnktbls <- fdata$lnktblFiles$outfile %>% session$dataWrangler$attachNames(fdata$lnktblFiles$dataset)

fdata$lnktbls <- fdata$lnktbls %>% session$collectionUtils$lapply(function(file) { readRDS(paste0(workspace$outputDir, file)) })

fdata$positives <- fdata$lnktbls %>% session$collectionUtils$lapply(function(lnktbl) { lnktbl %>% filter(lnk == 1) })
fdata$negatives <- fdata$lnktbls %>% session$collectionUtils$lapply(function(lnktbl) { lnktbl %>% filter(lnk == 0) })

fdata$datasets <- fdata$lnktbls %>% names() %>% session$dataWrangler$attachNames()

fdata$fishers <- fdata$datasets %>% session$collectionUtils$lapply(function(dataset_T) {
  
  datasets_P <- fdata$datasets[str_sub(fdata$datasets, 1, 3) != str_sub(dataset_T, 1, 3)]
  
  if (grepl("pres$", dataset_T)) {
    datasets_P <- datasets_P[grepl("pres$", datasets_P)]
  } else {
    datasets_P <- datasets_P[!(grepl("pres$", datasets_P))]
  }
  
  print(paste0("dataset_T = ", dataset_T, " vs. ", datasets_P))
  
  datasets_P %>% lapply(function(datasetP) {
    pred_T <- fdata$positives[[datasetP]]
    pred_F <- fdata$negatives[[datasetP]] 
    
    lnks_T <- fdata$positives[[dataset_T]]
    
    # compute fisher's in parallel for the 6 cell types... 
    session$evaluationUtils$fisher(predicted = pred_T$pair, 
                                   notPredicted = pred_F$pair, 
                                   trueSet = lnks_T$pair,
                                   alternative = "greater")
  })
  
})

# ++++++++++++++++++++++++
# COMMIT
fdata$fishers %>% saveRDS(paste0(workspace$outputDir, paste0("bk_repro_fisher.rds")))
# ++++++++++++++++++++++++++++++=
# save the formatted version

fdata <- list(); gc()

fdata$fishers <- read_rds(paste0(workspace$outputDir, paste0("bk_repro_fisher.rds")))

fdata$main <- fdata$fishers %>% session$collectionUtils$lapplyWithName(function(dataset_T, currfishers1) {
  
  currfishers1 %>% session$collectionUtils$lapplyWithName(function(dataset_P, currFisher) {
    
      currFisher$stats %>% as.data.frame() %>% t() %>% as_tibble() %>% 
        mutate(pvalue = currFisher$test$p.value) %>% 
        mutate(level = "bk", dataset = dataset_P, dat_truth = dataset_T) 

  }) %>% session$dataWrangler$rbind()
  
}) %>% session$dataWrangler$rbind() 

# ++++++++++++++++++++++++
# COMMIT
fdata$output$fishers <- fdata$fishers
fdata$output$tbl <- fdata$main
fdata$output %>% saveRDS(paste0(workspace$outputDir, paste0("bk_repro_fisher.rds")))
# ++++++++++++++++++++++++

# ++++++++++++++++++++++++++++++=
# gather all the lnktbls, compute fisher in parallel

fdata <- list(); gc()

fdata$lnktblFiles <- pdata$files$lnktbls %>% filter(level == "bk") %>% dplyr::select(everything(), dataset = cell_type)

fdata$datasets <- fdata$lnktblFiles$dataset
names(fdata$datasets) <- fdata$datasets

# fdata$datasets <- fdata$datasets[c("rosmap-ihc", "rosmap-ihcres")] # filter for only target dataset here

fdata$main <- fdata$datasets %>% lapply(function(currDataset) {
  
  print(paste0("working on dataset: ", currDataset))
  
  # set up current predictions; which is the bulk tissue networks; 
  # question: how well are the sc & sbj network edges recovered in the bulk tissue networks? 
  pred_lnktbl <- readRDS(paste0(workspace$outputDir, paste0("bk_lnktbls_", currDataset, ".rds")))
  pred_T <- pred_lnktbl %>% filter(lnk == 1)
  pred_F <- pred_lnktbl %>% filter(lnk == 0)
  
  # for each level
  list(sc = "sc", sbj = "sbj") %>% lapply(function(currLevel) {
    
    print(paste0("working on level: ", currLevel))
    
    currFiles <- pdata$files$lnktbls %>% filter(level == currLevel)
    celltypes <- currFiles$outfile
    names(celltypes) <- currFiles$cell_type
    
    print(paste0("running fisher's exact tests in parallel for the different cell types..."))

    celltypes %>% mclapply(function(currFile) {
      
      lnksDataset <- unlist(strsplit(currDataset, "-"))[1]
      lnks_T <- readRDS(paste0(workspace$outputDir, currFile))[[lnksDataset]] %>% filter(lnk == 1)

      # compute fisher's in parallel for the 6 cell types... 
      session$evaluationUtils$fisher(predicted = pred_T$pair, 
                                     notPredicted = pred_F$pair, 
                                     trueSet = lnks_T$pair, 
                                     alternative = "greater")
      
    }, mc.cores = length(celltypes))
    
    
  })
})

# ++++++++++++++++++++++++
# COMMIT
fdata$main %>% saveRDS(paste0(workspace$outputDir, paste0("bk_preservation_fisher.rds")))

# ++++++++++++++++
# save tbl version of fishers 

fdata <- list(); gc()

fdata$main <- readRDS(paste0(workspace$outputDir, paste0("bk_preservation_fisher.rds")))

fdata$main <- fdata$main %>% session$collectionUtils$lapplyWithName(function(currDataset, currfishers1) {
  
  currfishers1 %>% session$collectionUtils$lapplyWithName(function(currLevel, currfishers2) {
    
    currfishers2 %>% session$collectionUtils$lapplyWithName(function(currCelltype, currFisher) {
      currFisher$stats %>% as.data.frame() %>% t() %>% as_tibble() %>% 
        mutate(pvalue = currFisher$test$p.value) %>% 
        mutate(dataset = currDataset, level = currLevel, cell_type = currCelltype) 
    }) %>% session$dataWrangler$rbind()
    
  }) %>% session$dataWrangler$rbind()
  
}) %>% session$dataWrangler$rbind() 

# ++++++++++++++++++++++++
# COMMIT
fdata$output$fishers <- fdata$fishers
fdata$output$tbl <- fdata$main
fdata$output %>% saveRDS(paste0(workspace$outputDir, paste0("bk_preservation_fisher.rds")))
# ++++++++++++++++++++++++


# ++++++++++++++++++++++++++++++=
# look at how well xCell and xSubject are preserved in bulk tissue

fdata <- list(); gc()

fdata$main <- c(xcell_reproducibility =  "sc_repro_fisher.rds",
                xsbj_reproducibility = "sbj_repro_fisher.rds", 
                xsbj_preservation = "sc_sbj_preservation_fisher.rds", 
                xbulk_reproducibility = "bk_repro_fisher.rds", 
                xbulk_preservation = "bk_preservation_fisher.rds")

fdata$main <- fdata$main %>% session$collectionUtils$lapply(function(file) { read_rds(paste0(workspace$outputDir, file)) })

fdata$tbls <- fdata$main %>% session$collectionUtils$lapply(function(fisher) { fisher$tbl })

pdata$fishers <- fdata # SOFT COMMIT
# ++++++++++++++++++++++++++++++++++++++++++++++=

# 3. xBulk co-expression networks are somewhat reproducible (plot the other levels as well for these two datasets )

fdata <- list(); gc()

fdata$datasets <- c("velmeshev", "rosmap") %>% session$dataWrangler$attachNames()

fdata$tbls <- pdata$fishers$tbls

fdata$main <- fdata$tbls$xbulk_reproducibility %>% 
  mutate(cell_type = "bulk") %>% 
  mutate(level = "xbulk") %>% 
  dplyr::select(level, cell_type, dat_ref = dat_truth, dat_test = dataset, or, recovered_n, recovered_frac) %>% 
  mutate(cell_type = dat_ref %>% sapply(function(str) { if (grepl(".*res$", str)) { "bulk_mgpres" } else { "bulk" }})) %>% 
  mutate(dat_ref = str_extract(dat_ref, "^[a-z]+")) %>% 
  rbind(fdata$tbls$xsbj_reproducibility %>%  
          mutate(level = "xsbj") %>% 
          dplyr::select(level, cell_type, dat_ref = dat_truth, dat_test = dat_sample, or, recovered_n, recovered_frac) %>% 
          filter(dat_ref %in% fdata$datasets, dat_test %in% fdata$datasets) %>% 
          mutate(dat_ref = workspace$utils$fmtDataset(dat_ref), 
                 cell_type = workspace$utils$fmtCelltypes(cell_type))) %>% 
  rbind(fdata$tbls$xcell_reproducibility %>%  
          mutate(level = "xcell") %>% 
          dplyr::select(level, cell_type, dat_ref = dat_truth, dat_test = dat_sample, or, recovered_n, recovered_frac) %>% 
          filter(dat_ref %in% fdata$datasets, dat_test %in% fdata$datasets) %>% 
          mutate(dat_ref = workspace$utils$fmtDataset(dat_ref), 
                 cell_type = workspace$utils$fmtCelltypes(cell_type)))

fdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = or)) + 
  geom_bar(stat = "identity") + 
  geom_text(aes(label = round(or, 1)), vjust = -0.5, size = 5) + 
  facet_wrap(~dat_ref*level, scales = "free_x") + 
  session$graphingUtils$tiltX(angle = 90) + 
  ylim(0, 200)


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_01.eps"

xdata$main <- fdata$main %>% mutate(level = level %>% sapply(function(str) { 
  if (str == "xbulk") { "xBulk" } 
  else if (str == "xsbj") { "xSubject" } 
  else { "xCell" }
})) 

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = or)) + 
  geom_bar(stat = "identity") + 
  geom_text(aes(label = round(or, 1)), vjust = -0.5, size = 5) + 
  facet_wrap(~dat_ref*level, scales = "free_x") + 
  session$graphingUtils$tiltX(angle = 90) + 
  ylim(0, 200) +
  ggtitle("Reproducibility of co-expression edges between the Velmeshev and ROSMAP datasets") + 
  theme(plot.title = element_text(size = 20), plot.subtitle = element_text(size = 14)) + 
  ylab("Odds ratio")
  

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 20, height = 30, dpi = 1000) 
# =================================================






# =================================================
# OUTPUT FIGURE 
# TODO move this to the next figure !
xdata <- list()

xdata$figName <- "figure_0XXXX.eps"

xdata$main <- fdata$main %>% 
  mutate(dataset = workspace$utils$fmtDataset(dataset), 
         cell_type = workspace$utils$fmtCelltypes(cell_type)) %>% 
  mutate(level = level %>% sapply(function(str) {
    if (str == "sc") { "xCell" } 
    else { "xSubject" }
  }))
  
xdata$plot <- xdata$main %>% session$graphingUtils$ggplot(aes(x = cell_type, y = or)) + 
  geom_bar(aes(fill = level), position = "dodge", stat = "identity") + 
  facet_wrap(~dataset, ncol = 2) + 
  session$graphingUtils$tiltX(angle = 90) + 
  scale_fill_hue(l = 45) +
  xlab("Cell type") + 
  ylab("Preservation of co-expression edges\n(Odds ratio)") + 
  labs(fill = "Level")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 25, height = 13, dpi = 1000) 
# =================================================


# ++++++++++++================================
# ok, here plot xcell and xsubject reproducibility against preservation at the xbulk level 
# first for rosmap

fdata <- list(); gc()

fdata$dataset <- "rosmap"

fdata$repro <- pdata$fishers$tbls$xcell_reproducibility %>% filter(dat_truth == fdata$dataset) %>% 
  group_by(cell_type) %>% 
  summarize(or = mean(or)) %>% 
  mutate(level = "xCell") %>% 
  rbind(pdata$fishers$tbls$xsbj_reproducibility %>% filter(dat_truth == fdata$dataset) %>% 
          group_by(cell_type) %>% 
          summarize(or = mean(or)) %>% 
          mutate(level = "xSubject")) %>% 
  dplyr::select(cell_type, level, or)

fdata$preserv <- pdata$fishers$tbls$xbulk_preservation %>%
  filter(dataset %in% fdata$dataset) %>% 
  dplyr::select(cell_type, level, or) %>% 
  mutate(level = level %>% sapply(function(str) {
    if (str == "sc") { "xCell" } else { "xSubject" }
  }))


# =================================================
# OUTPUT FIGURE 
xdata <- list()

xdata$figName <- "figure_02.eps"

xdata$main <- fdata$repro %>% 
  dplyr::select(cell_type, level, reproducbility = or) %>% 
  left_join(fdata$preserv %>% dplyr::select(cell_type, level, preservation = or), by = c("cell_type", "level"))

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = preservation, y = reproducbility)) + 
  geom_point(aes(color = level, shape = cell_type), size = 4) + 
  geom_abline(slope = 1, linetype = "dashed") + 
  scale_color_manual(values = c("#bf8502", "#0576b5")) + 
  scale_x_continuous(trans = "log2", limits = c(2, 150)) + 
  scale_y_continuous(trans = "log2", limits = c(2, 150)) + 
  xlab("Preservation of co-expression at the xBulk level\n(Odds ratio - log2 scale)") + 
  ylab("Reproducibility of co-expression at the\nxCell and xSubject levels\n(Odds ratio - log2 scale)") +   
  ggtitle("Reproducibility of xCell or xSubject co-expression\nvs. their preservation at the xBulk level", "ROSMAP") 
  
xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 20, height = 15) 
# =================================================


# first for velmeshev

fdata <- list(); gc()

fdata$dataset <- "velmeshev"

fdata$repro <- pdata$fishers$tbls$xcell_reproducibility %>% filter(dat_truth == fdata$dataset) %>% 
  group_by(cell_type) %>% 
  summarize(or = mean(or)) %>% 
  mutate(level = "xCell") %>% 
  rbind(pdata$fishers$tbls$xsbj_reproducibility %>% filter(dat_truth == fdata$dataset) %>% 
          group_by(cell_type) %>% 
          summarize(or = mean(or)) %>% 
          mutate(level = "xSubject")) %>% 
  dplyr::select(cell_type, level, or)

fdata$preserv <- pdata$fishers$tbls$xbulk_preservation %>%
  filter(dataset %in% fdata$dataset) %>% 
  dplyr::select(cell_type, level, or) %>% 
  mutate(level = level %>% sapply(function(str) {
    if (str == "sc") { "xCell" } else { "xSubject" }
  }))


# =================================================
# OUTPUT FIGURE 
xdata <- list()

xdata$figName <- "figure_03.eps"

xdata$main <- fdata$repro %>% 
  dplyr::select(cell_type, level, reproducbility = or) %>% 
  left_join(fdata$preserv %>% dplyr::select(cell_type, level, preservation = or), by = c("cell_type", "level"))

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = preservation, y = reproducbility)) + 
  geom_point(aes(color = level, shape = cell_type), size = 4) + 
  geom_abline(slope = 1, linetype = "dashed") + 
  scale_color_manual(values = c("#bf8502", "#0576b5")) + 
  scale_x_continuous(trans = "log2", limits = c(2, 150)) + 
  scale_y_continuous(trans = "log2", limits = c(2, 150)) + 
  xlab("Preservation of co-expression at the xBulk level\n(Odds ratio - log2 scale)") + 
  ylab("Reproducibility of co-expression at the\nxCell and xSubject levels\n(Odds ratio - log2 scale)") +   
  ggtitle("Reproducibility of xCell or xSubject co-expression\nvs. their preservation at the xBulk level", "Velmeshev") 

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 20, height = 15) 
# =================================================


# ++++++++++++++++++++++++++++++++++++++++++++++=

fdata <- list(); gc()

fdata$main <- c(rosmap_ihc = "bk_exprmats_rosmap-ihcres.rds",
                rosmap_mgp = "bk_exprmats_rosmap-mgpres.rds", 
                velmeshev_mgp = "bk_exprmats_velmeshev-mgpres.rds")

fdata$main <- fdata$main %>% session$collectionUtils$lapply(function(file) { read_rds(paste0(workspace$outputDir, file)) })

# =======================
pdata$ccvs <- fdata # SOFT COMMIT ###### !!!!!!!!!!!!!!!!!!
# =====================


# ++++++++++++++++======================================
# let's make sure that the markers I used aren't severely co-expressed at the xSubject level
# check for "uniform distribution" of the co-expressions... 

fdata <- list()

fdata$mkrs <- pdata$mkrs$mkrsFlat

fdata$mkrCoexTbl <- fdata$mkrs %>% session$collectionUtils$lapplyWithName(function(celltype, mkrs) {
  
  currCoexmats <- pdata$files$coexmats %>% filter(level == "sbj", dataset %in% c("velmeshev", "rosmap"), cell_type == celltype)
  currCoexmats <- currCoexmats$outfile %>% session$dataWrangler$attachNames()
  
  currCoexmats %>% lapply(function(currCoexmat) {
    coexmat <- read_rds(paste0(workspace$outputDir, currCoexmat))
    genes <- mkrs %>% intersect(rownames(coexmat))
    coexmat <- coexmat[genes, genes]
    coextbl <- coexmat %>% workspace$utils$vectorize() %>% session$dataWrangler$vectorToTibble() %>% 
      dplyr::select(coex = value) %>% 
      mutate(coex_file = currCoexmat, cell_type = celltype) %>% 
      dplyr::select(cell_type, coex_file, coex)
    coextbl <- coexmat %>% workspace$utils$getCoords() %>% cbind(coextbl) %>% as_tibble()
    return(coextbl)
  }) %>% session$dataWrangler$rbind()
  
}) %>% session$dataWrangler$rbind()

# lets get rid of the to the most co-regulated genes in oligodendrocyte and see how things change

fdata$mkrCoexTbl %>%
  filter(cell_type == "oligodendrocyte", coex_file == "sbj_coexmats_velmeshev_oligodendrocyte.rds") %>% 
  filter(coex > 0.9) -> x
x$gene_a %>% 
  c(x$gene_b) %>% as_tibble() %>% group_by(value) %>% summarize(n = n()) %>% arrange(desc(n)) -> y

pdata$ccvs$main$velmeshev_mgp$mgps$rotations$oligodendrocyte %>% 
  session$dataWrangler$setRownameAsColumn("value") %>% dplyr::select(value, PC6) %>% 
  left_join(y, by = "value") %>% 
  session$graphingUtils$ggplot(aes(x = n, y = PC6)) + geom_point()

pdata$ccvs$main$rosmap_mgp$ccvModels$ctpMat

y %>% filter(n < 20) %>% session$dataWrangler$extractColumn("value") -> y

workspace$utils$computeMgps(list(oligo = y), pdata$ccvs$main$velmeshev_mgp$ccvModels$exprmats$orig) -> yy
yy$main %>% 
  session$dataWrangler$setRownameAsColumn("sample") %>% 
  mutate(orig = pdata$ccvs$main$velmeshev_mgp$mgps$estimates$oligodendrocyte[rownames(yy$main)]) -> z

z %>% session$graphingUtils$ggplot(aes(x = oligo, y = orig)) + geom_point()

cor(z$oligo, z$orig)

# =================================================
# OUTPUT FIGURE 
xdata <- list()

xdata$figName <- "sfigure_01.eps"

xdata$plot <- fdata$mkrCoexTbl %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_histogram() + 
  facet_wrap(~cell_type*coex_file) + 
  ggtitle("Within cell type co-expression among marker genes at the xSubject level")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 40, height = 30) 

# ================================================


# ++++++++++++++++++++++=======================================================
# the other thing is to check whether cell type specificity makes a gene's R^2 higher

fdata <- list(); gc()

fdata$minfc <- pdata$mkrs$minfc
fdata$minfc <- fdata$minfc %>% filter(!(gene_id %in% unlist(pdata$mkrs$mkrsFlat)))

# compute the variance
fdata$vars$velmeshev <- fdata$minfc %>% 
  filter(dataset == "velmeshev") %>% group_by(gene_id) %>% summarize(log2_expr_var = log2(var(expr))) %>% ungroup()
fdata$vars$rosmap <- fdata$minfc %>% 
  filter(dataset == "rosmap") %>% group_by(gene_id) %>% summarize(log2_expr_var = log2(var(expr))) %>% ungroup()


# for ihc
fdata$main$ihc <- pdata$ccvs$main$rosmap_ihc$ccvModel$stats$lmStats %>% dplyr::select(gene_id, rsqr)
fdata$main$ihc <- fdata$main$ihc %>% inner_join(fdata$vars$rosmap, by = "gene_id")
fdata$main$ihc <- fdata$main$ihc %>% mutate(method = "ihc")

# for rosmap-mgp
fdata$main$rosmapMgp <- pdata$ccvs$main$rosmap_mgp$ccvModels$stats$lmStats %>% dplyr::select(gene_id, rsqr)
fdata$main$rosmapMgp <- fdata$main$rosmapMgp %>% inner_join(fdata$vars$rosmap, by = "gene_id")
fdata$main$rosmapMgp <- fdata$main$rosmapMgp %>% mutate(method = "rosmap-mgp")

# for velmeshev-mgp
fdata$main$velmeshevMgp <- pdata$ccvs$main$velmeshev_mgp$ccvModels$stats$lmStats %>% dplyr::select(gene_id, rsqr)
fdata$main$velmeshevMgp <- fdata$main$velmeshevMgp %>% inner_join(fdata$vars$velmeshev, by = "gene_id")
fdata$main$velmeshevMgp <- fdata$main$velmeshevMgp %>% mutate(method = "velmeshev-mgp")

# all 
fdata$main$all <- fdata$main$ihc %>% rbind(fdata$main$rosmapMgp) %>% rbind(fdata$main$velmeshevMgp)

fdata$corVals <- fdata$main$all %>% 
  group_by(method) %>% summarize(cor_coef = cor(rsqr, log2_expr_var)[1, 1])

# =================================================
# OUTPUT FIGURE ---- ihc
xdata <- list()

xdata$figName <- "figure_04.png"
xdata$figName2 <- "figure_04_1.png"

xdata$main <- fdata$main$all %>% filter(method == "ihc")

xdata$plot2 <- xdata$main %>% session$graphingUtils$ggplot(aes(x = rsqr)) + geom_density(fill = "grey80") + 
  xlim(0, 1) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

xdata$plot2 

ggsave(filename = paste0(pdata$figuresDir, xdata$figName2), 
       plot = xdata$plot2, device = "png", units = "cm", width = 18, height = 10) 

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = log2_expr_var, y = rsqr)) + 
  geom_point(shape = 1) + 
  geom_smooth(method = "lm", se = FALSE) + 
  ylim(0, 1) + 
  xlab("Variation across cell types\nlog2(variance)") + 
  ylab(bquote("PVE by cell type proportions"~(R^2))) + 
  ggtitle("Association between cross cell type\nvariability and CCV influence across genes", "ROSMAP-IHC")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 18, height = 18) 
# =================================================

# =================================================
# OUTPUT FIGURE ---- rosmap-mgp
xdata <- list()

xdata$figName <- "figure_05.png"
xdata$figName2 <- "figure_05_1.png"

xdata$main <- fdata$main$all %>% filter(method == "rosmap-mgp")

xdata$plot2 <- xdata$main %>% session$graphingUtils$ggplot(aes(x = rsqr)) + geom_density(fill = "grey80") + 
  xlim(0, 1) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

xdata$plot2 

ggsave(filename = paste0(pdata$figuresDir, xdata$figName2), 
       plot = xdata$plot2, device = "png", units = "cm", width = 18, height = 10) 

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = log2_expr_var, y = rsqr)) + 
  geom_point(shape = 1) + 
  geom_smooth(method = "lm", se = FALSE) + 
  ylim(0, 1) + 
  xlab("Variation across cell types\nlog2(variance)") + 
  ylab(bquote("PVE by cell type proportions"~(R^2))) + 
  ggtitle("Association between cross cell type\nvariability and CCV influence across genes", "ROSMAP-MGP")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 18, height = 18) 
# =================================================

# =================================================
# OUTPUT FIGURE ---- ihc
xdata <- list()

xdata$figName <- "figure_06.png"
xdata$figName2 <- "figure_06_1.png"

xdata$main <- fdata$main$all %>% filter(method == "velmeshev-mgp")

xdata$plot2 <- xdata$main %>% session$graphingUtils$ggplot(aes(x = rsqr)) + geom_density(fill = "grey80") + 
  xlim(0, 1) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

xdata$plot2 

ggsave(filename = paste0(pdata$figuresDir, xdata$figName2), 
       plot = xdata$plot2, device = "png", units = "cm", width = 18, height = 10) 

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = log2_expr_var, y = rsqr)) + 
  geom_point(shape = 1) + 
  geom_smooth(method = "lm", se = FALSE) + 
  ylim(0, 1) + 
  xlab("Variation across cell types\nlog2(variance)") + 
  ylab(bquote("PVE by cell type proportions"~(R^2))) + 
  ggtitle("Association between cross cell type\nvariability and CCV influence across genes", "Velmeshev-MGP")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 18, height = 18) 
# =================================================



# ++++++++++++++++++++++++++++++++++++++++++++++=
# 5. CCV correction reduces concordance with xCell network, but improves concordance with xSubject

# let's first start with the IHC dataset 

fdata <- list(); gc()

fdata$main <- pdata$ccvs$main$velmeshev_mgp$ccvModel

fdata$stats <- fdata$main$stats %>% lapply(na.omit)

xdata <- list()

xdata$proportions <- fdata$main$ctpMat %>% 
  session$dataWrangler$setRownameAsColumn("sample") %>% 
  gather(cell_type, proportion, -sample)

xdata$proportions %>% 
  session$graphingUtils$ggplot(aes(x = sample, y = proportion, fill = cell_type)) + 
  geom_bar(stat = "identity") + 
  session$graphingUtils$tiltX(angle = 90) + 
  ylab("Cell type proportion")

xdata <- list()

xdata$mean <- fdata$stats$lmStats$rsqr %>% mean()
  
fdata$stats$lmStats %>% session$graphingUtils$ggplot(aes(x = rsqr)) + geom_histogram() + 
  geom_vline(xintercept = xdata$mean, linetype = "dashed", color = "red", size = 1) + 
  ggtitle("R^2 distribution - IHC")

fdata$stats$lmStats %>% session$graphingUtils$ggplot(aes(x = pvalue)) + geom_histogram() + 
  ggtitle("P-value distribution - IHC")

(fdata$sigGenes <- fdata$stats$lmStats %>% filter(qvalue < 0.1))

fdata$stats$coefStats %>% 
session$graphingUtils$ggplot(aes(x = rsqr_indep)) + geom_histogram() + facet_wrap(~cell_type)

# sanity check that beta values are positive --- genes would be positively correlated with cell type proportions
fdata$stats$coefStats %>% mutate(significant = gene_id %in% fdata$sigGenes$gene_id) %>% 
  session$graphingUtils$ggplot(aes(x = beta, y = significant)) + 
  geom_boxplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  facet_wrap(~cell_type)

fdata$stats$coefStats %>% filter(gene_id %in% fdata$sigGenes$gene_id) %>% arrange(desc(rsqr_indep)) %>% left_join(pdata$genes, by = "gene_id") %>% View()

# TWO GENES TO HIGHLIGHT AS EXAMPLES
# ENSG00000144285 - SCN2A
# ENSG00000198780 - FAM169A 
# pdata$mkrs$minfc %>% filter(gene_id == "ENSG00000198780")

# plot these genes as functions of cell type proportions

xdata <- list()

xdata$ccv <- fdata$main$ctpMat %>% session$dataWrangler$setRownameAsColumn("sample") %>% 
  gather(cell_type, proportion, -sample) 

xdata$gene <- "ENSG00000144285"
xdata$geneSymbol <- pdata$genes %>% filter(gene_id == xdata$gene)

xdata$expr <- fdata$main$exprmats$orig[xdata$gene,] %>% session$dataWrangler$vectorToTibble() %>% dplyr::select(sample = variable, expr = value)

xdata$ccv %>% left_join(xdata$expr, by = "sample") %>% 
  session$graphingUtils$ggplot(aes(x = proportion , y = expr)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  facet_wrap(~cell_type) +
  ggtitle(paste0("Expression of ", xdata$geneSymbol$gene, "\nby neuron fractions (IHC estimates)"))

# plot the co-expression of these genes before and after

xdata <- list()

xdata$genes <- pdata$genes %>% filter(gene_id %in% c("ENSG00000144285", "ENSG00000198780"))

xdata$expr <- fdata$main$exprmats$residual[xdata$genes$gene_id,] %>% t() %>% session$dataWrangler$setRownameAsColumn("sample")
names(xdata$expr) <- c("sample", xdata$genes$gene)

xdata$corCoef <- xdata$expr$FAM169A %>% cor(xdata$expr$SCN1A)

xdata$expr %>%
  session$graphingUtils$ggplot(aes(x = SCN1A, y = FAM169A)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  ggtitle(paste0("Pearson's = ", round(xdata$corCoef, 2)))

pdata$fishers$tbls$xbulk_preservation %>%  
  filter(dataset %in% c("rosmap-ihc", "rosmap-ihcres")) %>% 
  session$graphingUtils$ggplot(aes(x = or, y = cell_type, color = dataset)) + 
  geom_point(size = 5) + 
  facet_wrap(~level, nrow = 2)


# ++++++++++++++++++++++++++++====
# make the rsqr + pvalue distributions here

fdata <- list(); gc()

fdata$main <- pdata$ccvs$main %>% session$collectionUtils$lapplyWithName(function(name, currDat) {
  currDat$ccvModel$stats$lmStats %>% mutate(dataset = name)
}) %>% session$dataWrangler$rbind()


# =================================================
# OUTPUT FIGURE # HERE you can get the rsqr means
xdata <- list()

# xdata$figName <- "figure_02.eps"

xdata$main <- fdata$main %>% na.omit()

xdata$means <- xdata$main %>% group_by(dataset) %>% summarize(rsqr = mean(rsqr))

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = rsqr)) + 
  geom_density(aes(color = dataset), linewidth = 1) +
  geom_vline(data = xdata$means, aes(xintercept = rsqr, color = dataset), linetype = "dashed") + 
  scale_color_brewer(palette = "Dark2")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 20, height = 10) 

xdata$main %>% mutate(sig = qvalue < 0.1) %>% group_by(dataset, sig) %>% summarize(n = n()) # number of significant genes after multiple test correction

# =================================================


# ++++++++++++++++++++++++++++====
# correlations between MGP and IHC estimates in ROSMAP

fdata <- list(); gc()

fdata$ctpmat$ihc <- pdata$ccvs$main$rosmap_ihc$ccvModel$ctpMat
fdata$ctpmat$mgp <- pdata$ccvs$main$rosmap_mgp$ccvModels$ctpMat

fdata$samples <- rownames(fdata$ctpmat$ihc) %>% intersect(rownames(fdata$ctpmat$mg))

fdata$ctpmat <- fdata$ctpmat %>% lapply(function(mat) { mat[fdata$samples, sort(colnames(mat))] })

fdata$ctpmat %>% lapply(dim)

fdata$main <- fdata$ctpmat$mgp %>% cor(fdata$ctpmat$ihc)

fdata$celltypes$ihc <- colnames(fdata$ctpmat$ihc) %>% session$dataWrangler$attachNames()
fdata$celltypes$mgps <- colnames(fdata$ctpmat$mgp) %>% session$dataWrangler$attachNames()

fdata$main <- fdata$celltypes$ihc %>% session$collectionUtils$lapply(function(celltypeIhc) {
  fdata$celltypes$mgps %>% session$collectionUtils$lapply(function(celltypeMpg) {
    corTest <- fdata$ctpmat$ihc[, celltypeIhc] %>% cor.test(fdata$ctpmat$mgp[, celltypeMpg], alternative = "greater")
    
    tibble(ihc_celltype = celltypeIhc, mgp_celltype = celltypeMpg, 
           cor_coef = corTest$estimate, pvalue = corTest$p.value)
  }) %>% session$dataWrangler$rbind()
}) %>% session$dataWrangler$rbind()

fdata$main <- fdata$main %>% mutate(qvalue = p.adjust(pvalue, method = "fdr"))

# =================================================
# OUTPUT FIGURE -------- Paul: don't need to show this
# maybe stick it into the supplementaries so we're upfront about everything
xdata <- list()

xdata$figName <- "sfigure_03.eps"

xdata$main <- fdata$main %>% 
  mutate(mgp_celltype = workspace$utils$fmtCelltypes(mgp_celltype), 
         ihc_celltype = workspace$utils$fmtCelltypes(ihc_celltype)) %>% 
  mutate(mgp_celltype = factor(mgp_celltype), ihc_celltype = factor(ihc_celltype)) 

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = mgp_celltype, y = ihc_celltype)) + 
  geom_tile(aes(fill = cor_coef)) + 
  geom_tile(data = xdata$main %>% filter(pvalue <= 0.05), aes(fill = cor_coef), colour = "black", linewidth = 1) +
  geom_text(aes(label = round(cor_coef, 2)), size = 6) +
  scale_fill_gradient(low = "white", high = "red") + 
  session$graphingUtils$tiltX(angle = 90) + 
  ylim(workspace$utils$fmtCelltypes(c("astrocyte", "neuron", "microglia", "oligodendrocyte", "endothelial"))) + 
  xlab("Cell type compositionestimated by MPG") + 
  ylab("Cell type composition\nestimated by IHC")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 20, height = 13) 

# =================================================



# +++++++++++++++++++++++++++++++++++++++++++++++++++++====================
# see how those change in Pearson's correlation before and after correction

# rosmap.ihc vs ihcres
# rosmap mgp vs mgp res
# velmeshev mgp vs mgp res

fdata <- list(); gc()

fdata$datasets <- list(rosmap = "rosmap", velmeshev = "velmeshev")

fdata$lnktbls$xsbj <- pdata$celltypes %>% mclapply(function(celltype) {
  file <- paste0("sbj_lnktbls_", celltype, ".rds")
  lnktbl <- read_rds(paste0(workspace$outputDir, file))
  lnktbl <- lnktbl[names(fdata$datasets)]
  lnktbl <- lnktbl %>% lapply(function(tbl) { tbl %>% filter(lnk == 1) })
  return(lnktbl)
}, mc.cores = length(pdata$celltypes))

# TODO - run this 
fdata$lnktbls$xcell <- pdata$celltypes %>% session$collectionUtils$lapply(function(celltype) { 
  file <- paste0("sc_lnktbls_", celltype, ".rds")
  lnktbl <- read_rds(paste0(workspace$outputDir, file))
  lnktbl <- lnktbl[names(fdata$datasets)]
  lnktbl <- lnktbl %>% lapply(function(tbl) { tbl %>% filter(lnk == 1) })
  return(lnktbl)
})

fdata$corcoef <- list(rosmap.ihc = "bk_coefmats_rosmap-ihc.rds", 
                      rosmap.ihcres = "bk_coefmats_rosmap-ihcres.rds", 
                      rosmap = "bk_coefmats_rosmap.rds", 
                      rosmap.mgpres = "bk_coefmats_rosmap-mgpres.rds", 
                      velmeshev = "bk_coefmats_velmeshev.rds", 
                      velmeshev.mgpres = "bk_coefmats_velmeshev-mgpres.rds")

fdata$corcoef <- fdata$corcoef %>% mclapply(function(file) {
  read_rds(paste0(workspace$outputDir, file))
}, mc.cores = length(fdata$corcoef))

fdata$corcoefVec <- fdata$corcoef %>% mclapply(function(coexmat) {
  coexmat %>% workspace$utils$vectorize() %>% 
    session$dataWrangler$attachNames(workspace$utils$getPairIds(coexmat))
}, mc.cores = length(fdata$corcoef)) 

# save a copy of fdata to accelerate loading tomorrow
# fdata %>% saveRDS(paste0(workspace$outputDir, "fdata_20231017.rds")); quit(save = "no")
# fdata <- read_rds(paste0(workspace$outputDir, "fdata_20231017.rds"))

# for every xsubject network links, compute the correlations in bulk networks, matching dataset
fdata$main <- fdata$datasets %>% lapply(function(dataset) {
  
  print(paste0("working on ", dataset, "..."))
  
  corcoefVec <- fdata$corcoefVec[grepl(dataset, names(fdata$corcoefVec))]
  
  pdata$celltypes %>% lapply(function(celltype) {
    
    print(paste0("working on ", celltype, "..."))
    
    fdata$lnktbls %>% session$collectionUtils$lapplyWithName(function(level, lnktbls) {
      
      print(paste0("working on ", level, "..."))
      
      lnktbl <- lnktbls[[celltype]][[dataset]]
      
      corcoefVec %>% session$collectionUtils$lapplyWithName(function(network, coexVec) {
        
        print(paste0("working on ", network, "..."))
        
        coexlnks <- coexVec[lnktbl$pair] %>% na.omit() 
        tibble(dataset = dataset, cell_type = celltype, level = level, network = network, 
               n_gene = length(coexlnks), coef_mean = mean(coexlnks))
        
      }, verbose = FALSE) %>% session$dataWrangler$rbind()
      
    }, verbose = FALSE) %>% session$dataWrangler$rbind()

  }) %>% session$dataWrangler$rbind()

}) %>% session$dataWrangler$rbind()

# save the results!!!
fdata$main %>% saveRDS(paste0(workspace$outputDir, "bk_preservation_coef.rds")) # load this for analysis later



# +++++++++++++++++++++++++++++++++++++++++++++++++++++====================
# plot out how the co-expressions change as a result of corrective procedure 

fdata <- list(); gc()

fdata$main <- read_rds(paste0(workspace$outputDir, "bk_preservation_coef.rds")) # load this for analysis later

fdata$main %>% filter(network %in% c("rosmap.ihc", "rosmap.ihcres")) %>% 
  mutate(level = factor(level, levels = c("xsbj", "xcell"))) %>% 
  session$graphingUtils$ggplot(aes(y = cell_type, x = coef_mean)) + 
  geom_point(aes(color = network), size = 5) +
  facet_wrap(~level, ncol = 1) 

fdata$diffs <- fdata$main %>% filter(network %in% c("rosmap.ihc", "rosmap.ihcres")) %>% 
  dplyr::select(cell_type, level, network, coef_mean) %>% 
  spread(network, coef_mean) %>% 
  dplyr::select(everything(), orig = rosmap.ihc, res = rosmap.ihcres) %>% 
  mutate(diff = res - orig, comparison = "ROSMAP-IHC") %>% 
  rbind(fdata$main %>% filter(network %in% c("rosmap", "rosmap.mgpres")) %>% 
          dplyr::select(cell_type, level, network, coef_mean) %>% 
          spread(network, coef_mean) %>% 
          dplyr::select(everything(), orig = rosmap, res = rosmap.mgpres) %>% 
          mutate(diff = res - orig, comparison = "ROSMAP-MPG")) %>% 
  rbind(fdata$main %>% filter(network %in% c("velmeshev", "velmeshev.mgpres")) %>% 
          dplyr::select(cell_type, level, network, coef_mean) %>% 
          spread(network, coef_mean) %>% 
          dplyr::select(everything(), orig = velmeshev, res = velmeshev.mgpres) %>%
          mutate(diff = res - orig, comparison = "Velmeshev-MGP")) 
  

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_07.eps"

xdata$main <- fdata$diffs %>% 
  mutate(cell_type = workspace$utils$fmtCelltypes(cell_type)) %>% 
  mutate(level = level %>% sapply(function(str) {
    if (str == "xcell") { "xCell" }
    else { "xSubject" }
  }))

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = level, y = diff, group = comparison, fill = level)) + 
  geom_bar(position = position_dodge2(preserve = "single", padding = 0.4), stat = "identity") + 
  facet_wrap(~cell_type, nrow = 1) +
  scale_fill_manual(values = c("#bf8502", "#0576b5")) + 
  session$graphingUtils$tiltX(angle = 90) + 
  ylab("Change in average co-expression\nafter CCV correction")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 35, height = 13) 

# =================================================


fdata <- list(); gc()

fdata$main <- pdata$fishers$tbls$xbulk_preservation %>% filter(grepl("-ihc.*$", dataset))

fdata$means <- fdata$main %>% 
  group_by(dataset, level) %>% 
  summarize(or_mean = mean(or)) %>% 
  ungroup()


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_08.eps"

xdata$main <- fdata$main %>% 
  mutate(cell_type = workspace$utils$fmtCelltypes(cell_type)) %>% 
  mutate(level = level %>% sapply(function(str) {
    if (str == "sc") { "xCell" } 
    else { "xSubject" }
  })) %>% 
  mutate(dataset = dataset %>% sapply(function(str) {
    if (str == "rosmap") { "ROSMAP" } 
    else if (str == "rosmap-mgppres") { "ROSMAP MGP-residual" }
    else if (str == "velmeshev") { "Velmeshev" }
    else if (str == "velmeshev-mgpres") { "Velmeshev MGP-residual" }
    else if (str == "rosmap-ihc") { "ROSMAP IHC Samples" }
    else if (str == "rosmap-ihcres") { "ROSMAP IHC-residual" }
    else { "" }
  }))

xdata$main <- xdata$main %>% mutate(dataset = factor(dataset, levels = 
                                         c("ROSMAP IHC Samples", 
                                           "ROSMAP IHC-residual",
                                           "ROSMAP",
                                           "ROSMAP MGP-residual", 
                                           "Velmeshev",
                                           "Velmeshev MGP-residual")))

xdata$plot <- xdata$main %>% session$graphingUtils$ggplot(aes(x = cell_type, y = or)) + 
  geom_bar(aes(fill = level), position = "dodge", stat = "identity") + 
  facet_wrap(~dataset, nrow = 1) + 
  session$graphingUtils$tiltX(angle = 90) + 
  scale_fill_manual(values = c("#bf8502", "#0576b5")) + 
  xlab("Cell type") + 
  ylab("Preservation of co-expression edges\n(Odds ratio)") + 
  labs(fill = "Level")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 45, height = 18, dpi = 1000) 
# =================================================
















# ++++++++++++++++++++++++++++++++++++++++++++++=
# now look at how well MGP corrected the whole network for ROSMAP

# let's first start with the IHC dataset 

fdata <- list(); gc()

fdata$main <- pdata$ccvs$main$rosmap_mgp$ccvModel

fdata$stats <- fdata$main$stats %>% lapply(na.omit)


# compare the MGP estimates to the IHC estimates

xdata <- list()

xdata$ihc <- pdata$ccvs$main$rosmap_ihc$ccvModel$ctpMat
xdata$mgp <- fdata$main$ctpMat

xdata$samples <- rownames(xdata$ihc) %>% intersect(rownames(xdata$mgp))

xdata$ihc <- xdata$ihc[xdata$samples, ]
xdata$mgp <- xdata$mgp[xdata$samples, ]

xdata$main <- xdata$ihc %>% cor(xdata$mgp)

xdata$main <- xdata$main[c("neuron", "oligodendrocyte", "astrocyte", "microglia", "endothelial"), 
           c("excitatory", "inhibitory", "opc",  "oligodendrocyte", "astrocyte", "microglia")]


# Variance explained by MGPs

xdata <- list()

(xdata$mean <- fdata$stats$lmStats$rsqr %>% mean())

fdata$stats$lmStats %>% session$graphingUtils$ggplot(aes(x = rsqr)) + geom_histogram() + 
  geom_vline(xintercept = xdata$mean, linetype = "dashed", color = "red", size = 1) + 
  ggtitle("R^2 distribution - MGP (ROSMAP)")

fdata$stats$lmStats %>% session$graphingUtils$ggplot(aes(x = pvalue)) + geom_histogram() + 
  ggtitle("P-value distribution - MGP (ROSMAP)")

fdata$stats$lmStats %>% summary()

(fdata$sigGenes <- fdata$stats$lmStats %>% filter(qvalue < 0.1))

fdata$stats$coefStats %>% 
  session$graphingUtils$ggplot(aes(x = rsqr_indep)) + geom_histogram() + facet_wrap(~cell_type)

# sanity check that beta values are positive --- genes would be positively correlated with cell type proportions
fdata$stats$coefStats %>% mutate(significant = gene_id %in% fdata$sigGenes$gene_id) %>% 
  session$graphingUtils$ggplot(aes(x = beta, y = significant)) + 
  geom_boxplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  facet_wrap(~cell_type)

fdata$stats$coefStats %>% filter(gene_id %in% fdata$sigGenes$gene_id) %>% arrange(desc(rsqr_indep)) %>% left_join(pdata$genes, by = "gene_id") %>% View()

# TWO GENES TO HIGHLIGHT AS EXAMPLES
# ENSG00000144285 - SCN2A
# ENSG00000198780 - FAM169A 
# pdata$mkrs$minfc %>% filter(gene_id == "ENSG00000198780")

# plot these genes as functions of cell type proportions

xdata <- list()

xdata$ccv <- fdata$main$ctpMat %>% session$dataWrangler$setRownameAsColumn("sample") %>% 
  gather(cell_type, proportion, -sample) 

xdata$gene <- "ENSG00000144285"
xdata$geneSymbol <- pdata$genes %>% filter(gene_id == xdata$gene)

xdata$expr <- fdata$main$exprmats$orig[xdata$gene,] %>% session$dataWrangler$vectorToTibble() %>% dplyr::select(sample = variable, expr = value)

xdata$ccv %>% left_join(xdata$expr, by = "sample") %>% 
  session$graphingUtils$ggplot(aes(x = proportion , y = expr)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  facet_wrap(~cell_type)
  ggtitle(paste0("Expression of ", xdata$geneSymbol$gene, "\nby neuron fractions (IHC estimates)"))

# plot the co-expression of these genes before and after

xdata <- list()

xdata$genes <- pdata$genes %>% filter(gene_id %in% c("ENSG00000144285", "ENSG00000198780"))

xdata$expr <- fdata$main$exprmats$orig[xdata$genes$gene_id,] %>% t() %>% session$dataWrangler$setRownameAsColumn("sample")
names(xdata$expr) <- c("sample", xdata$genes$gene)

xdata$corCoef <- xdata$expr$FAM169A %>% cor(xdata$expr$SCN1A)

xdata$expr %>%
  session$graphingUtils$ggplot(aes(x = SCN1A, y = FAM169A)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  ggtitle(paste0("Pearson's = ", round(xdata$corCoef, 2)))


pdata$fishers$tbls$xbulk_preservation %>%  
  filter(dataset %in% c("rosmap", "rosmap-mgppres")) %>% 
  session$graphingUtils$ggplot(aes(x = or, y = cell_type, color = dataset)) + 
  geom_point(size = 5) + 
  facet_wrap(~level, nrow = 2)



# ++++++++++++++++++++++++++++++++++++++++++++++=
# now look at how well MGP corrected the whole network for VELMESHEV

# let's first start with the IHC dataset 

fdata <- list(); gc()

fdata$main <- pdata$ccvs$main$velmeshev_mgp$ccvModel

fdata$stats <- fdata$main$stats %>% lapply(na.omit)

# Variance explained by MGPs

xdata <- list()

(xdata$mean <- fdata$stats$lmStats$rsqr %>% mean())

fdata$stats$lmStats %>% session$graphingUtils$ggplot(aes(x = rsqr)) + geom_histogram() + 
  geom_vline(xintercept = xdata$mean, linetype = "dashed", color = "red", size = 1) + 
  ggtitle("R^2 distribution - MGP (Velmeshev)")

fdata$stats$lmStats %>% session$graphingUtils$ggplot(aes(x = pvalue)) + geom_histogram() + 
  ggtitle("P-value distribution - MGP (Velmeshev)")

fdata$stats$lmStats %>% summary()

(fdata$sigGenes <- fdata$stats$lmStats %>% filter(qvalue < 0.1))

fdata$stats$coefStats %>% 
  session$graphingUtils$ggplot(aes(x = rsqr_indep)) + geom_histogram() + facet_wrap(~cell_type)

# sanity check that beta values are positive --- genes would be positively correlated with cell type proportions
fdata$stats$coefStats %>% mutate(significant = gene_id %in% fdata$sigGenes$gene_id) %>% 
  session$graphingUtils$ggplot(aes(x = beta, y = significant)) + 
  geom_boxplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  facet_wrap(~cell_type)

# TWO GENES TO HIGHLIGHT AS EXAMPLES
# ENSG00000144285 - SCN2A
# ENSG00000198780 - FAM169A 
# pdata$mkrs$minfc %>% filter(gene_id == "ENSG00000198780")

# plot these genes as functions of cell type proportions

xdata <- list()

xdata$ccv <- fdata$main$ctpMat %>% session$dataWrangler$setRownameAsColumn("sample") %>% 
  gather(cell_type, proportion, -sample) 

xdata$gene <- "ENSG00000198780"
xdata$geneSymbol <- pdata$genes %>% filter(gene_id == xdata$gene)

xdata$expr <- fdata$main$exprmats$orig[xdata$gene,] %>% session$dataWrangler$vectorToTibble() %>% dplyr::select(sample = variable, expr = value)

xdata$ccv %>% left_join(xdata$expr, by = "sample") %>% 
  session$graphingUtils$ggplot(aes(x = proportion , y = expr)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  facet_wrap(~cell_type) +
  ggtitle(paste0("Expression of ", xdata$geneSymbol$gene, "\nby neuron fractions (MGP estimates)"))

# plot the co-expression of these genes before and after

xdata <- list()

xdata$genes <- pdata$genes %>% filter(gene_id %in% c("ENSG00000144285", "ENSG00000198780"))

xdata$expr <- fdata$main$exprmats$residual[xdata$genes$gene_id,] %>% t() %>% session$dataWrangler$setRownameAsColumn("sample")
names(xdata$expr) <- c("sample", xdata$genes$gene)

xdata$corCoef <- xdata$expr$FAM169A %>% cor(xdata$expr$SCN1A)

xdata$expr %>%
  session$graphingUtils$ggplot(aes(x = SCN1A, y = FAM169A)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  ggtitle(paste0("Pearson's = ", round(xdata$corCoef, 2)))

pdata$fishers$tbls$xbulk_preservation %>%  
  filter(dataset %in% c("velmeshev", "velmeshev-mgpres")) %>% 
  session$graphingUtils$ggplot(aes(x = or, y = cell_type, color = dataset)) + 
  geom_point(size = 5) + 
  facet_wrap(~level, nrow = 2)


# ++++++++++++++++++++++++++++++++++++++++++++++=
# Drive home the point by summarizing odds ratio after corrective procedures

fdata <- list(); gc()

fdata$main <- pdata$fishers$tbls$xbulk_preservation %>% filter(dataset %in% c("velmeshev-mgpres", "rosmap-mgppres", "rosmap-ihcres"))

fdata$main %>% 
  session$graphingUtils$ggplot(aes(x = or, y = cell_type, color = level)) + 
  geom_point(size = 4) + 
  facet_wrap(~dataset) + 
  xlab("Odds ratio")

fdata$main %>% group_by(dataset, level) %>% summarize(mean(or)) # maybe add mean values to the figure above

# +++++++++++++++++++++++++++++++++++++++++++===
# include a comparison of... 















  
  


# TODO SANDBOX HERE TO BE DELETED - AD HOC VISUALS FOR PAUL 

fdata <- list()

fdata$dataset <- "lim"
fdata$celltype <- "excitatory"

fdata$main <- pdata$fishers$tbls$xcell_reproducibility %>% filter(dat_truth == fdata$dataset, cell_type == fdata$celltype) %>% 
  dplyr::select(dataset = dat_truth, odds_ratio = or) %>% 
  mutate(test = "xcell_repro") %>% 
  rbind(pdata$fishers$tbls$xsbj_preservation %>% filter(dat_truth == fdata$dataset, cell_type == fdata$celltype) %>% 
          dplyr::select(dataset = dat_truth, odds_ratio = or) %>% 
          mutate(test = "xcell-xsbj_preserv")) %>% 
  rbind(pdata$fishers$tbls$xsbj_reproducibility %>% filter(dat_truth == fdata$dataset, cell_type == fdata$celltype) %>% 
          dplyr::select(dataset = dat_truth, odds_ratio = or) %>% 
          mutate(test = "xsbj_repro")) %>% 
  rbind(pdata$fishers$tbls$xbulk_preservation %>% filter(dataset == fdata$dataset, cell_type == fdata$celltype, level == "sbj") %>% 
          dplyr::select(dataset, odds_ratio = or) %>% 
          mutate(test = "xsbj_xbulk_pres"))
  # rbind(pdata$fishers$tbls$xbulk_reproducibility %>% filter(dat_truth == fdata$dataset) %>% 
  #         dplyr::select(dataset = dat_truth, odds_ratio = or) %>% 
  #         mutate(test = "xbulk_repro"))

fdata$main$test <- factor(fdata$main$test, levels = rev(c("xcell_repro", "xcell-xsbj_preserv", "xsbj_repro", "xsbj_xbulk_pres", "xbulk_repro")))

fdata$main %>% 
  session$graphingUtils$ggplot(aes(x = odds_ratio, y = test)) +
  geom_point(size = 8, alpha = 0.5, color = "black") +
  scale_x_continuous(trans = "log10") 





















