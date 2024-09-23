# in this script perform comparison of co-expression networks with literature curations
# 1. convert the DTRIs into co-regulated pairs (just those with 3 types of experimental evidence) - positive & negative sets
# 2. start with 1 co-expression network, compute AUROC enrichment for those pairs 

pdata <- list()

pdata$figuresDir <- paste0(workspace$workspaceDir, "041_ANALYSIS_07_dtri_comparison_FIGURES/")

pdata$genes <- readRDS(paste0(workspace$outputDir, "genes_metadata.rds")) # genes

pdata$utils$constructCoregmat <- function(coregTf) {
  coregMat <- coregTf %>% dplyr::select(target_a, target_b) %>% mutate(coreg = 1) %>% unique() %>% spread(target_b, coreg) %>% 
    session$dataWrangler$setColAsRownames("target_a") 
  coregMat <- coregMat[sort(rownames(coregMat)), sort(colnames(coregMat))] %>% as.matrix()
  diag(coregMat) <- 0 # the fact that a gene is "co-regulated" with itself is trivial
  coregMat[is.na(coregMat)] <- 0
  # all(rownames(coregMat) == colnames(coregMat))
  # isSymmetric(coregMat) # making sure the matrix is symmetric
  return(coregMat)
}

pdata$ctprofiles <- read_rds(paste0(workspace$outputDir, "sc_stats_ctprofiles.rds"))

# +++++++++++++++++++++++++++++++++ 
# 1. load in the DTRIs
# filter for those with 3 types of experimental evidence
# convert to using ensembl human symbols
# convert to co-regulated pairs 

# combine the two sheets - assign tri new ids
fdata <- list(); gc() 

fdata$curation$plos <- paste0(workspace$miscDir, "dtri_plos.tsv") %>% read.table(sep = "\t", skip = 1, header = TRUE) %>% as_tibble()
fdata$curation$expand <- paste0(workspace$miscDir, "dtri_exp.tsv") %>% read.table(sep = "\t", header = TRUE) %>% as_tibble()

# filter for CNS terms
fdata$cnsTerms <- session$ontologyUtils$getChildNodeRecursive(node = "UBERON_0001017", scope = session$ontologyUtils$ontologies$UBERON) %>%
  c(session$ontologyUtils$getChildNodeRecursive(node = "UBERON_0001017", scope = session$ontologyUtils$ontologies$CL)) %>%
  unique()

fdata$dtris$plos <- fdata$curation$plos %>% filter(Cell_Type %in% fdata$cnsTerms)

# fix formating in ontology terms, replace ":" with "_"
fdata$dtris$expand <- fdata$curation$expand %>% mutate(Cell_Type = Cell_Type %>% sapply(function(str) { sub(":", "_", str) })) 
fdata$dtris$expand <- fdata$dtris$expand %>% filter(Cell_Type %in% fdata$cnsTerms)

fdata$dtris$plos <- fdata$dtris$plos %>% 
  dplyr::select(tf_entrez = TF_Entrez_ID_Human, target_entrez = Target_Entrez_ID_Human, exp_type = Experiment_Type) %>% 
  unique() %>% 
  na.omit()

# convert all mouse entrez to human entrez
fdata$dtris$expand <- fdata$dtris$expand %>% 
  dplyr::select(tf_entrez = TF_Entrez_ID, target_entrez = Target_Entrez_ID, exp_type = Experiment_Type)

fdata$dtris$expand <- fdata$dtris$expand %>% mutate(exp_type = exp_type %>% sapply(function(typeStr) {
  if (grepl("binding", tolower(typeStr))) { "binding" }
  else if (grepl("per", tolower(typeStr))) { "perturbation" }
  else if (grepl("reporter", tolower(typeStr))) { "reporter" }
  else NA
})) %>% na.omit() 

xdata <- list()
xdata$orthologs <- homologene::mouse2human(unique(c(fdata$dtris$expand$tf_entrez, fdata$dtris$expand$target_entrez))) %>% as_tibble()
# retain only one to one mappings
xdata$mouseGeneV <- xdata$orthologs %>% group_by(mouseID) %>% summarize(n = n()) %>% arrange(desc(n)) %>% filter(n == 1)
xdata$humanGeneV <- xdata$orthologs %>% group_by(humanID) %>% summarize(n = n()) %>% arrange(desc(n)) %>% filter(n == 1)
xdata$orthologs <- xdata$orthologs %>% filter(mouseID %in% xdata$mouseGeneV$mouseID, humanID %in% xdata$humanGeneV$humanID)

# convert tf_entrez to human versions
fdata$dtris$expand <- fdata$dtris$expand %>% left_join(xdata$orthologs %>% dplyr::select(tf_entrez = mouseID, tf_entrez_h = humanID), by = "tf_entrez") %>% 
  dplyr::select(tf_entrez.x = tf_entrez_h, tf_entrez.y = tf_entrez, everything()) %>%
  session$dataWrangler$mergeColumnsXy("tf_entrez") 

# convert target_entrez to human versions
fdata$dtris$expand <- fdata$dtris$expand %>% left_join(xdata$orthologs %>% dplyr::select(target_entrez = mouseID, target_entrez_h = humanID), by = "target_entrez") %>% 
  dplyr::select(target_entrez.x = target_entrez_h, target_entrez.y = target_entrez, everything()) %>%
  session$dataWrangler$mergeColumnsXy("target_entrez") 

fdata$dtris$expand <- fdata$dtris$expand %>% dplyr::select(tf_entrez, target_entrez, exp_type) %>% unique() %>% na.omit()

# now combine the two dtri tbls
# ; 1499, original, 4387 after expansion -- total
# ; 330 in total after filtering for CNS annotations
fdata$dtris$all <- fdata$dtris$plos %>% rbind(fdata$dtris$expand) %>% mutate(tf_entrez = as.character(tf_entrez), target_entrez = as.character(target_entrez))

# filter for all the entrez IDs that are present in our platform
fdata$dtris$all <- fdata$dtris$all %>% filter(tf_entrez %in% pdata$genes$entrez, target_entrez %in% pdata$genes$entrez)

# convert entrez ids to ensembl
fdata$dtris$all <- fdata$dtris$all %>% 
  left_join(pdata$genes %>% dplyr::select(tf_entrez = entrez, tf_ensembl = gene_id), by = "tf_entrez") %>% 
  left_join(pdata$genes %>% dplyr::select(target_entrez = entrez, target_ensembl = gene_id), by = "target_entrez") %>% 
  dplyr::select(tf_ensembl, target_ensembl, exp_type) %>% 
  unique() %>% na.omit()

# final network # 4,083 DTRIs in total; 
# 288 DTRIs after filtering + mapping to ensembl IDs
fdata$dtris$all <- fdata$dtris$all %>% group_by(tf_ensembl, target_ensembl) %>% summarize(n_exp_type = n()) %>% ungroup() %>% unique() %>% na.omit()

# COMMIT
fdata$output <- list()
fdata$output$curation <- fdata$curation
fdata$output$dtris <- fdata$dtris
fdata$output %>% saveRDS(paste0(workspace$outputDir, "dtri_curation.rds"))


# ++++++++++++++++++++++++++++++++
# Now, rearrange the DTRIs into co-regulated pairs

fdata <- list(); gc()

fdata$curation <- read_rds(paste0(workspace$outputDir, "dtri_curation.rds"))

fdata$dtris <- fdata$curation$dtris$all

fdata$dtris <- fdata$dtris %>% filter(n_exp_type >= 1) # 288 interactions

fdata$tfs <- fdata$dtris %>% group_by(tf_ensembl) %>% summarize(n_target = n()) %>% arrange(desc(n_target)) %>% filter(n_target >= 2) # 43 TFs 
fdata$tfs$n_target %>% sum() # 247 DTRIs

fdata$coregTbl <- fdata$tfs$tf_ensembl %>% lapply(function(currTf) {
  targets <- fdata$dtris %>% filter(tf_ensembl == currTf)
  targets <- targets$target_ensembl
  coregTbl <- expand.grid(targets, targets) %>% as_tibble() %>% mutate(tf = currTf)
  coregTbl <- coregTbl %>% dplyr::select(tf, target_a = Var1, target_b = Var2)
  return(coregTbl)
}) %>% session$dataWrangler$rbind()

# make it into a matrix
fdata$coregMat <- fdata$coregTbl %>% dplyr::select(target_a, target_b) %>% mutate(co_reg = 1) %>% unique() %>% spread(target_b, co_reg) %>% 
  session$dataWrangler$setColAsRownames("target_a") 

fdata$coregMat <- fdata$coregMat[sort(rownames(fdata$coregMat)), sort(colnames(fdata$coregMat))] %>% as.matrix()
diag(fdata$coregMat) <- 0 # the fact that a gene is "co-regulated" with itself is trivial
rownames(fdata$coregMat) == colnames(fdata$coregMat)
isSymmetric(fdata$coregMat) # making sure the matrix is symmetric
fdata$coregMat[is.na(fdata$coregMat)] <- 0

fdata$coregMat %>% dim() # 184 genes in "co-regulated" pairs

# now turn the matrix back into a tibble, but this time using the sorted ensembl genes + lower.tri that way ids will align
fdata$coregTbl_F <- tibble(pair_id = fdata$coregMat %>% workspace$utils$getPairIds(),
                           coreg = fdata$coregMat %>% workspace$utils$vectorize()) 

fdata$coregTbl_F %>% group_by(coreg) %>% summarize(n()) # 1480 "co-regulated pairs"; 15356 negatives

fdata$coregTbl <- fdata$coregTbl 

# COMMIT
fdata$output <- list()
fdata$output$dtris <- fdata$dtris
fdata$output$coregTf <- fdata$coregTbl
fdata$output$coregmat <- fdata$coregMat
fdata$output$coregTbl <- fdata$coregTbl_F
fdata$output %>% saveRDS(paste0(workspace$outputDir, "dtri_coreg.rds")) # COMMIT TO DISK!!!

# +++++++++++++++++++++++++++++++++++++++++++++++++++++
# next ---- test it out with just a single network & and then write the code to run it on all the co-expression networks...
# might be able to "fast-track" it by taking only the submatrix containing the positive and negative sets

fdata <- list(); gc()

fdata$coreg <- read_rds(paste0(workspace$outputDir, "dtri_coreg.rds"))

fdata$coregtf <- fdata$coreg$coregTf

fdata$tfs <- fdata$coregtf$tf %>% unique() %>% sort()

# fdata$genes <- fdata$coreg$coregmat %>% rownames()

fdata$coexmat <- read_rds(paste0(workspace$outputDir, "sc_coexmats_rosmap_oligodendrocyte.rds"))

fdata$tfsV <- fdata$tfs #%>% intersect(rownames(fdata$coexmat))

fdata$coregtf <- fdata$coregtf %>% filter(tf %in% fdata$tfsV)

fdata$coregmat <- fdata$coregtf %>% pdata$utils$constructCoregmat(); fdata$coregmat %>% isSymmetric() # making sure the matrix is symmetric
  
fdata$coregtbl <- tibble(pair_id = fdata$coregmat %>% workspace$utils$getPairIds(), coreg = fdata$coregmat %>% workspace$utils$vectorize()) 

fdata$genes <- rownames(fdata$coregmat) %>% intersect(rownames(fdata$coexmat)) %>% sort()

fdata$coexmat <- fdata$coexmat[fdata$genes, fdata$genes]

fdata$coextbl <- tibble(pair_id = workspace$utils$getPairIds(fdata$coexmat), 
                        coex = workspace$utils$vectorize(fdata$coexmat))

fdata$coregtbl <- fdata$coregtbl %>% filter(pair_id %in% fdata$coextbl$pair_id)

fdata$coextbl %>% 
  arrange(desc(coex)) %>% 
  left_join(fdata$coregtbl) %>% 
  mutate(coreg = as.logical(coreg)) %>% 
  session$graphingUtils$ggplot(aes(x = coreg, y = coex)) + geom_violin()

session$evaluationUtils$roc(ranking = (fdata$coextbl %>% arrange(desc(coex)))$pair_id, 
                            trueSet = (fdata$coregtbl %>% filter(coreg == 1))$pair_id, 
                            step = 50) -> x

x %>% 
  session$graphingUtils$ggplot(aes(x = false_positive_rate, y = recall)) + 
  geom_line(group = 1) + 
  geom_abline(slope = 1, linetype = "dashed")

session$evaluationUtils$auroc((fdata$coextbl %>% arrange(desc(coex)))$pair_id, 
                              trueSet = (fdata$coregtbl %>% filter(coreg == 1))$pair_id)


# +++++++++++++++++++++++++++++++++++++++++++++++++++++
# here - generate a null using shuffled tf-target pairs, do it 100 times

fdata <- list(); gc()

fdata$coreg <- read_rds(paste0(workspace$outputDir, "dtri_coreg.rds"))

fdata$coregtf <- fdata$coreg$coregTf

fdata$tfs <- fdata$coregtf$tf %>% unique() %>% sort()

# fdata$genes <- fdata$coreg$coregmat %>% rownames()

fdata$datasetName <- "sc_coexmats_ramos_excitatory.rds"

fdata$coexmat <- read_rds(paste0(workspace$outputDir, fdata$datasetName))

fdata$tfsV <- fdata$tfs #%>% intersect(rownames(fdata$coexmat))

fdata$coregtf <- fdata$coregtf %>% filter(tf %in% fdata$tfsV)

fdata$coregmat <- fdata$coregtf %>% pdata$utils$constructCoregmat(); fdata$coregmat %>% isSymmetric() # making sure the matrix is symmetric

fdata$coregtbl <- tibble(pair_id = fdata$coregmat %>% workspace$utils$getPairIds(), coreg = fdata$coregmat %>% workspace$utils$vectorize()) 



fdata$genes <- rownames(fdata$coregmat) %>% intersect(rownames(fdata$coexmat)) %>% sort()

fdata$coexmat <- fdata$coexmat[fdata$genes, fdata$genes]

fdata$coextbl <- tibble(pair_id = workspace$utils$getPairIds(fdata$coexmat), 
                        coex = workspace$utils$vectorize(fdata$coexmat))

fdata$coregtbl <- fdata$coregtbl %>% filter(pair_id %in% fdata$coextbl$pair_id)

fdata$pairsPos <- fdata$coregtbl %>% filter(coreg == 1) %>% mutate(pair_str = strsplit(pair_id, "\\.")) %>% 
  mutate(target_a = pair_str %>% sapply(function(str) { str[1] }), 
         target_b = pair_str %>% sapply(function(str) { str[2] })) %>% 
  dplyr::select(pair_id, target_a, target_b)


set.seed(0); fdata$main <- 1:1000 %>% session$collectionUtils$lapply(function(i) {
  
  currPairsPos <- fdata$pairsPos %>% mutate(target_b = sample(target_b, replace = FALSE))
  # currPairsPos <- fdata$pairsPos
  
  currPairsPos <- currPairsPos %>% 
    mutate(pair_id = paste0(target_a, ".", target_b)) %>% 
    dplyr::select(-pair_id)
  
  currPairsPos <- currPairsPos %>% rbind(tibble(target_a = currPairsPos$target_b, target_b = currPairsPos$target_a)) # flip it 
  
  currCoremat <- currPairsPos %>% pdata$utils$constructCoregmat(); currCoremat %>% isSymmetric()
  currCoregtbl <- tibble(pair_id = currCoremat %>% workspace$utils$getPairIds(), coreg = currCoremat %>% workspace$utils$vectorize()) 
  
  auroc <- session$evaluationUtils$auroc((fdata$coextbl %>% arrange(desc(coex)))$pair_id, 
                                         trueSet = (currCoregtbl %>% filter(coreg == 1))$pair_id)
  
  
  auroc %>% session$dataWrangler$vectorToTibble() %>% spread(variable, value) %>% mutate(iteration = i)
  
}) %>% session$dataWrangler$rbind()

(auroc <- session$evaluationUtils$auroc((fdata$coextbl %>% arrange(desc(coex)))$pair_id, 
                                       trueSet = fdata$pairsPos$pair_id ))

fdata$main %>% 
  session$graphingUtils$ggplot(aes(x = auroc.W)) + 
  geom_histogram() + 
  geom_vline(xintercept = mean(fdata$main$auroc.W)) + 
  geom_vline(xintercept = auroc["auroc.W"]) + 
  ggtitle(fdata$datasetName, 
          paste0("Null distribution generated by shuffling target pairs 1000 times \nMean = ", round(mean(fdata$main$auroc.W), 2),
                     "\nMin = ", round(min(fdata$main$auroc.W), 2), 
                     "\nMax = ", round(max(fdata$main$auroc.W), 2), 
                     "\nObserved = ", round(auroc["auroc.W"], 2), 
                 "; ", nrow(fdata$main %>% filter(auroc.W >= round(auroc["auroc.W"], 2))), "/", nrow(fdata$main), " > Observed"))

# +++++++++++++++++++++++++++++===============================
# load in the file names of co-expression matrices

pdata$files$coexmats <- tibble(outfile = list.files(workspace$outputDir)) %>% 
  filter(grepl("coexmats", outfile)) %>% 
  mutate(outfile_name = gsub(".rds", "", outfile)) %>% 
  mutate(strs = strsplit(outfile_name, "_")) %>% 
  mutate(level = strs %>% sapply(function(currStr) { currStr[1] })) %>% 
  mutate(type = strs %>% sapply(function(currStr) { currStr[2] })) %>% 
  mutate(dataset = strs %>% sapply(function(currStr) { currStr[3] })) %>% 
  mutate(cell_type = strs %>% sapply(function(currStr) { currStr[4] })) %>% 
  dplyr::select(outfile, level, type, dataset, cell_type)


# +++++++++++++++++++++++++++++++++++++==
# now ramp it up to all the coexmats 
# first, put together a list of coexmat files, then go through each and generate the numbers 
# 1. number of genes overlapping, 2. number of pairs (positive and negatve), 3. AUROC, 4. p-value

# here do it for all TFs - next - do it for TFs detected in the particular co-expression netowrk

fdata <- list(); gc()

fdata$coreg <- read_rds(paste0(workspace$outputDir, "dtri_coreg.rds"))
fdata$coregGenes <- fdata$coreg$coregmat %>% rownames()

fdata$coexmatFiles <- pdata$files$coexmats %>% filter(!grepl("psdbk", dataset), dataset != "consensus") # 86 networks

fdata$coexmatFiles <- fdata$coexmatFiles$outfile %>% session$dataWrangler$attachNames()

set.seed(1); fdata$main <- fdata$coexmatFiles %>% session$collectionUtils$lapply(function(coexmatFile) {
  
  coexmat <- read_rds(paste0(workspace$outputDir, coexmatFile))
  genes <- fdata$coregGenes %>% intersect(rownames(coexmat)) %>% unique() %>% sort()
  
  coexmat <- coexmat[genes, genes]
  coextbl <- tibble(pair_id = workspace$utils$getPairIds(coexmat), coex = workspace$utils$vectorize(coexmat))
  coregtbl <- fdata$coreg$coregTbl %>% filter(pair_id %in% coextbl$pair_id)
  
  auroc <- session$evaluationUtils$auroc((coextbl %>% arrange(desc(coex)))$pair_id, trueSet = (coregtbl %>% filter(coreg == 1))$pair_id)
  
  output <- list()
  
  # construct nulls here & compare it to the observed value to obtain a p-value
  output$null <- 1:1000 %>% session$collectionUtils$lapply(function(i) {
    
    pairsPos <- coregtbl %>% filter(coreg == 1) %>%
      mutate(pair_str = strsplit(pair_id, "\\.")) %>% 
      mutate(target_a = pair_str %>% sapply(function(str) { str[1] }), 
             target_b = pair_str %>% sapply(function(str) { str[2] })) %>% 
      dplyr::select(pair_id, target_a, target_b)
    currPairsPos <- pairsPos %>% mutate(target_b = sample(target_b, replace = FALSE))

    currPairsPos <- currPairsPos %>% 
      mutate(pair_id = paste0(target_a, ".", target_b)) %>% 
      dplyr::select(-pair_id)
    
    currPairsPos <- currPairsPos %>% rbind(tibble(target_a = currPairsPos$target_b, target_b = currPairsPos$target_a)) # flip it 
    
    currCoremat <- currPairsPos %>% pdata$utils$constructCoregmat(); currCoremat %>% isSymmetric()
    currCoregtbl <- tibble(pair_id = currCoremat %>% workspace$utils$getPairIds(), coreg = currCoremat %>% workspace$utils$vectorize()) 
    
    auroc <- session$evaluationUtils$auroc((coextbl %>% filter(pair_id %in% currCoregtbl$pair_id) %>% arrange(desc(coex)))$pair_id, 
                                           trueSet = (currCoregtbl %>% filter(coreg == 1))$pair_id)
    
    auroc %>% 
      session$dataWrangler$vectorToTibble() %>% 
      spread(variable, value) %>% mutate(iteration = i) %>%
      mutate(coex_file = coexmatFile)
      
  }) %>% session$dataWrangler$rbind()
  
  output$observed <- tibble(coex_file = coexmatFile,
         n_gene = length(genes), 
         n_pairs_full = auroc[2], 
         n_pairs_positive = auroc[1], 
         auc = auroc[3], 
         pvalue = auroc[4])
  
  return(output)
})

fdata$nulls <- fdata$main %>% lapply(function(item) { item$null }) %>% 
  session$dataWrangler$rbind() %>% 
    left_join(pdata$files$coexmats %>% dplyr::select(coex_file = outfile, level, dataset, cell_type))

fdata$observeds <- fdata$main %>% lapply(function(item) { 
    currObserved <- item$observed 
    pvalueObs <- fdata$nulls %>% filter(coex_file == currObserved$coex_file) %>% filter(auroc.W > currObserved$auc) %>% nrow()
    currObserved %>% mutate(pvalue_em = (pvalueObs / 1000)) # work put the empirical pvalues
  }) %>% 
  session$dataWrangler$rbind() %>% 
  left_join(pdata$files$coexmats %>% dplyr::select(coex_file = outfile, level, dataset, cell_type))

fdata$observeds <- fdata$observeds %>% mutate(qvalue = p.adjust(pvalue, method = "fdr"))

# COMMIT
fdata$observeds %>% saveRDS(paste0(workspace$outputDir, "dtri_cns_auroc.rds"))


# +++++++++++++++++++==========
# visualize what I found with this analysis

fdata <- list(); gc()

fdata$main <- read_rds(paste0(workspace$outputDir, "dtri_cns_auroc.rds"))

fdata$main %>%
  filter(level == "bk") %>% 
  session$graphingUtils$ggplot(aes(x = pvalue_em)) + geom_histogram(bins = 50)

fdata$main %>%
  filter(level == "sc") %>% filter(qvalue < 0.1)

fdata$main %>% 
  filter(level == "sc") %>% group_by(auc > 0.5) %>% summarize(n = n())

fdata$main %>%
  filter(level == "sc") %>% 
  dplyr::select(dataset, cell_type, auc) %>% 
  spread(cell_type, auc) %>% 
  session$dataWrangler$setColAsRownames("dataset") %>% 
  session$graphingUtils$heatmap(cluster_rows = FALSE, cluster_cols = FALSE)


# =================================================
# OUTPUT FIGURE --- xCell level DTRI AUROCs
xdata <- list()

xdata$figName <- "figure_01"

xdata$main <- fdata$main %>% filter(level %in% c("sc")) 
xdata$main <- xdata$main %>% mutate(cell_type = workspace$utils$fmtCelltypes(cell_type), 
                                    dataset = workspace$utils$fmtDataset(dataset))

xdata$sigs_em <- xdata$main %>% filter(pvalue_em <= 0.05)
xdata$sigs <- xdata$main %>% filter(pvalue <= 0.05)

xdata$plot <- xdata$main %>% 
  dplyr::select(dataset, cell_type, auc) %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = auc)) +
  geom_boxplot(color = "grey40") + 
  geom_point(data = xdata$sigs, shape = 1, size = 6) + 
  geom_point(data = xdata$sigs_em, shape = 4, size = 5) + 
  geom_point(aes(color = dataset), size = 3) + 
  scale_color_brewer(palette = "Set1") + 
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(x = "Cell type", y = "Level of enrichment (AUROC)", color = "Dataset") +
  session$graphingUtils$tiltX(angle = 90) + 
  ggtitle("Enrichment for curated co-regulated gene pairs", "xCell level") + 
  ylim(0.45, 0.57)

xdata$plot

# ggsave(filename = paste0(pdata$figuresDir, xdata$figName, ".png"), 
#        plot = xdata$plot, device = "png", units = "cm", width = 24, height = 19, dpi = 400) 

ggsave(filename = paste0(pdata$figuresDir, xdata$figName, ".eps"), 
       plot = xdata$plot, device = "eps", units = "cm", width = 17, height = 20) 

# =================================================


# =================================================
# OUTPUT FIGURE --- xSubject level DTRI AUROCs
xdata <- list()  

xdata$figName <- "figure_02"

xdata$main <- fdata$main %>% filter(level %in% c("sbj")) 
xdata$main <- xdata$main %>% mutate(cell_type = workspace$utils$fmtCelltypes(cell_type), 
                                    dataset = workspace$utils$fmtDataset(dataset))

xdata$sigs_em <- xdata$main %>% filter(pvalue_em <= 0.05)
xdata$sigs <- xdata$main %>% filter(pvalue <= 0.05)

xdata$plot <- xdata$main %>% 
  dplyr::select(dataset, cell_type, auc) %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = auc)) +
  geom_boxplot(color = "grey40") + 
  geom_point(data = xdata$sigs, shape = 1, size = 6) + 
  geom_point(data = xdata$sigs_em, shape = 4, size = 5) + 
  geom_point(aes(color = dataset), size = 3) + 
  scale_color_brewer(palette = "Set1") + 
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(x = "Cell type", y = "Level of enrichment (AUROC)", color = "Dataset") +
  session$graphingUtils$tiltX(angle = 90) + 
  ggtitle("Enrichment for curated co-regulated gene pairs", "xSubject level") + 
  ylim(0.45, 0.57) 

xdata$plot

# ggsave(filename = paste0(pdata$figuresDir, xdata$figName, ".png"), 
#        plot = xdata$plot, device = "png", units = "cm", width = 24, height = 19, dpi = 400) 

ggsave(filename = paste0(pdata$figuresDir, xdata$figName, ".eps"), 
       plot = xdata$plot, device = "eps", units = "cm", width = 17, height = 20) 

# =================================================


# =================================================
# OUTPUT FIGURE 
xdata <- list()

xdata$figName <- "figure_03.eps"
xdata$figName2 <- "figure_03_1.eps"


xdata$main <- fdata$main %>% 
  filter(level %in% c("sc", "sbj")) %>% 
  dplyr::select(dataset, level, cell_type, auc) %>% 
  spread(level, auc)

xdata$main <- xdata$main %>% mutate(cell_type = workspace$utils$fmtCelltypes(cell_type), 
                                    dataset = workspace$utils$fmtDataset(dataset))

# need to compute the paired wilcoxon test
xdata$pvalues$all <- wilcox.test(xdata$main$sbj, xdata$main$sc, paired = TRUE, alternative = "less")$p.value
# 0.0001713364

xdata$celltypes <- xdata$main$cell_type %>% unique() %>% sort() %>% session$dataWrangler$attachNames()
xdata$pvalues$celltypes <- xdata$celltypes %>% sapply(function(celltype) { 
  dat <- xdata$main %>% filter(cell_type == celltype)
  wilcox.test(dat$sbj, dat$sc, paired = TRUE, alternative = "less")$p.value
}) %>% session$dataWrangler$vectorToTibble() %>% 
  dplyr::select(cell_type = variable, pvalue = value) 

xdata$plot2 <- xdata$pvalues$celltypes %>% 
  session$graphingUtils$ggplot(aes(y = cell_type, x = -log10(pvalue))) + 
  geom_bar(stat = "identity") + 
  geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
  xlab("-log10(P-Value)") 
xdata$plot2 

ggsave(filename = paste0(pdata$figuresDir, xdata$figName2), 
       plot = xdata$plot, device = "eps", units = "cm", width = 15, height = 10) 

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = sc, y = sbj)) + 
  geom_point(aes(color = dataset, shape = cell_type), size = 4) + 
  geom_abline(slope = 1, linetype = "dashed") +
  xlim(0.45, 0.57) + 
  ylim(0.45, 0.57) +
  xlab("xCell level enrichment (AUROC)") + 
  ylab("xSubject level enrichment (AUROC)") +   
  ggtitle("Enrichment for co-regulated target pairs in xCell vs. xSubject networks") +
  scale_color_brewer(palette = "Set1") 

xdata$plot

# ggsave(filename = paste0(pdata$figuresDir, xdata$figName, ".png"), 
#        plot = xdata$plot, device = "png", units = "cm", width = 24, height = 19, dpi = 400) 

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 20, height = 16) 

# =================================================

# average auc decreased from 0.514 to 0.501
fdata$main %>% 
  filter(level %in% c("sc", "sbj")) %>% 
  dplyr::select(dataset, level, cell_type, auc) %>% group_by(level) %>% summarize(auc = mean(auc))


# =================================================
# OUTPUT FIGURE --- xBulk level DTRI AUROCs
xdata <- list()

xdata$figName <- "figure_04"

xdata$main <- fdata$main %>% filter(level %in% c("bk")) 
xdata$main <- xdata$main %>% mutate(ccv_crrted = coex_file %>% sapply(function(file) { grepl("res.rds", file) }))

xdata$main <- xdata$main %>% mutate(dataset = dataset %>% sapply(function(str) {
  if (str == "rosmap-ihc") { "ROSMAP IHC" }
  else if (str == "rosmap-ihcres") { "ROSMAP IHC-residual" }
  else if (str == "rosmap") { "ROSMAP" }
  else if (str == "rosmap-mgpres") { "ROSMAP MGP-residual" }
  else if (str == "velmeshev") { "Velmeshev" }
  else if (str == "velmeshev-mgpres") { "Velmeshev MGP-residual" }
  else { NA }
}))

xdata$sigs_em <- xdata$main %>% filter(pvalue_em <= 0.05)
xdata$sigs <- xdata$main %>% filter(pvalue <= 0.05)

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = dataset, y = auc)) + 
  geom_point(data = xdata$sigs, shape = 1, size = 6) + 
  geom_point(data = xdata$sigs_em, shape = 4, size = 5) + 
  geom_point(aes(color = ccv_crrted), size = 3) + 
  scale_color_brewer(palette = "Paired") + 
  session$graphingUtils$tiltX(angle = 90) + 
  xlim(c("ROSMAP IHC", 
         "ROSMAP IHC-residual", 
         "ROSMAP", 
         "ROSMAP MGP-residual", 
         "Velmeshev", 
         "Velmeshev MGP-residual")) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  ylim(0.49, 0.55) + 
  xlab("Dataset") + 
  ylab("Level of enrichment (AUROC)") + 
  ggtitle("Enrichment for validated co-regulated gene pairs", "xBulk level") 

xdata$plot

# ggsave(filename = paste0(pdata$figuresDir, xdata$figName, ".png"), 
#        plot = xdata$plot, device = "png", units = "cm", width = 24, height = 19, dpi = 400) 

ggsave(filename = paste0(pdata$figuresDir, xdata$figName, ".eps"), 
       plot = xdata$plot, device = "eps", units = "cm", width = 17, height = 20) 

# =================================================





  
  


# +++++++++++++++++++++=== 
# let's get the co-expression values for the true pairs to see if there is anything reproducible in particular and what genes are these

fdata <- list(); gc()

fdata$coreg <- read_rds(paste0(workspace$outputDir, "dtri_coreg.rds"))
fdata$coregGenes <- fdata$coreg$coregmat %>% rownames()

fdata$coexmatFiles <- pdata$files$coexmats %>% filter(!grepl("psdbk", dataset), dataset != "consensus") # 86 networks

fdata$coexmatFiles <- fdata$coexmatFiles$outfile %>% session$dataWrangler$attachNames()

fdata$main <- fdata$coexmatFiles %>% session$collectionUtils$lapply(function(coexmatFile) {
  
  coexmat <- read_rds(paste0(workspace$outputDir, coexmatFile))
  genes <- fdata$coregGenes %>% intersect(rownames(coexmat)) %>% unique() %>% sort()
  
  coexmat <- coexmat[genes, genes]
  coextbl <- tibble(pair_id = workspace$utils$getPairIds(coexmat), coex = workspace$utils$vectorize(coexmat))
  coregtbl <- fdata$coreg$coregTbl %>% filter(pair_id %in% coextbl$pair_id)
  
  coregtbl %>% filter(coreg == 1) %>% left_join(coextbl, by = "pair_id") %>% mutate(coex_file = coexmatFile)
  
}) %>% session$dataWrangler$rbind()

fdata$main <- fdata$main %>% left_join(pdata$files$coexmats %>% dplyr::select(coex_file = outfile, level, dataset, cell_type))

# COMMIT
fdata$main %>% saveRDS(paste0(workspace$outputDir, "dtri_coreg_coex_values.rds"))



# +++++++++++++++++++++=== 
# let's get the coef values for the bulk co-expression networks
# this is going to be more informative for what correction does

fdata <- list(); gc()

fdata$coreg <- read_rds(paste0(workspace$outputDir, "dtri_coreg.rds"))
fdata$coregGenes <- fdata$coreg$coregmat %>% rownames()

fdata$coexmatFiles <- tibble(outfile = list.files(workspace$outputDir)) %>% 
  filter(grepl("coefmats_", outfile)) %>% 
  mutate(outfile_name = gsub(".rds", "", outfile)) %>% 
  mutate(strs = strsplit(outfile_name, "_")) %>% 
  mutate(level = strs %>% sapply(function(currStr) { currStr[1] })) %>% 
  mutate(type = strs %>% sapply(function(currStr) { currStr[2] })) %>% 
  mutate(dataset = strs %>% sapply(function(currStr) { currStr[3] })) %>% 
  mutate(cell_type = strs %>% sapply(function(currStr) { currStr[4] })) %>% 
  dplyr::select(outfile, level, type, dataset, cell_type)

fdata$coexmatFiles <- fdata$coexmatFiles$outfile %>% session$dataWrangler$attachNames()

fdata$main <- fdata$coexmatFiles %>% session$collectionUtils$lapply(function(coexmatFile) {
  
  coexmat <- read_rds(paste0(workspace$outputDir, coexmatFile))
  genes <- fdata$coregGenes %>% intersect(rownames(coexmat)) %>% unique() %>% sort()
  
  coexmat <- coexmat[genes, genes]
  coextbl <- tibble(pair_id = workspace$utils$getPairIds(coexmat), coex = workspace$utils$vectorize(coexmat))
  coregtbl <- fdata$coreg$coregTbl %>% filter(pair_id %in% coextbl$pair_id)
  
  coregtbl %>% filter(coreg == 1) %>% left_join(coextbl, by = "pair_id") %>% mutate(coex_file = coexmatFile)
  
}) %>% session$dataWrangler$rbind()


# COMMIT
fdata$main %>% saveRDS(paste0(workspace$outputDir, "dtri_coreg_coef_values.rds"))




# +++++++++++++++++++++++++++++++++++++===
# let's try to identify the co-regulated pairs with the highest co-expressions in each cell type

fdata <- list(); gc()

fdata$coreg <- read_rds(paste0(workspace$outputDir, "dtri_coreg.rds"))
fdata$curation <- read_rds(paste0(workspace$outputDir, "dtri_curation.rds"))

fdata$main <- read_rds(paste0(workspace$outputDir, "dtri_coreg_coex_values.rds"))
fdata$mainCoef <- read_rds(paste0(workspace$outputDir, "dtri_coreg_coef_values.rds"))



# Do this for every cell type? 

fdata$smry$sc <- fdata$main %>% 
  filter(level == "sc") %>% 
  group_by(cell_type, pair_id) %>% 
  summarize(coex = sum(coex)) %>% arrange(desc(coex)) %>% 
  ungroup()

fdata$smry$sbj <- fdata$main %>% 
  filter(level == "sbj") %>% 
  group_by(cell_type, pair_id) %>% 
  summarize(coex = sum(coex)) %>% arrange(desc(coex)) %>% 
  ungroup()

fdata$smry$sc %>% inner_join(fdata$smry$sbj, by = c("cell_type", "pair_id")) %>% 
  arrange(desc(coex.y)) %>% 
  session$graphingUtils$ggplot(aes(x = coex.x, y = coex.y)) + 
  geom_point() + facet_wrap(~cell_type) + 
  xlim(0, 7) + ylim(0, 7) + 
  geom_abline(slope = 1, linetype = "dashed") + 
  geom_vline(xintercept = 6) + 
  geom_hline(yintercept = 6) + 
  xlab("xCell coex (sum)") + 
  ylab("xSbj coex (sum)") + 
  ggtitle("", "co-expression of annotated co-regulated pairs in CNS")

fdata$smry$sc %>% inner_join(fdata$smry$sbj, by = c("cell_type", "pair_id")) %>% 
  group_by(cell_type) %>% summarize(coex.x = mean(coex.x), coex.y = mean(coex.y))

fdata$smry$sc %>% filter(cell_type == "inhibitory") %>% arrange(desc(coex))

# ENSG00000184347.ENSG00000169855 # very interesting case with 2 potential regulators - excitatory neuron???
# ENSG00000148737.ENSG00000110436 # very cool example with "astrocyte" specific biology being picked up ; regulation by PAX6
# ENSG00000123560.ENSG00000105695 # *** this is really nice as a demo for oligodendrocyte specific biology
# ENSG00000136750.ENSG00000128683 # *** this is an example with inhibitory neuron specific biology
# ENSG00000189056.ENSG00000169855 # this example again shows up due to enrichment in inhibitory neurons though may not be inhibitory neuron specific
# ENSG00000154928.ENSG00000118733 # *** THIS IS AN example where things completely fall apart in xSubject and xBulk networks!!!!!


xdata <- list()
xdata$pair <- "ENSG00000118733.ENSG00000105810"

fdata$main %>% 
  filter(level %in% c("sc", "sbj"), pair_id == xdata$pair) %>% 
  dplyr::select(dataset, level, cell_type, coex) %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = coex)) + 
  geom_boxplot() +
  geom_point(aes(color = dataset), size = 2, position = "jitter") + 
  facet_wrap(~level) + 
  session$graphingUtils$tiltX(angle = 90)

fdata$main %>% 
  filter(level %in% c("sc", "sbj"), pair_id == xdata$pair) %>% 
  dplyr::select(dataset, level, cell_type, coex) %>% 
  mutate(dataset = paste0(level, ".", dataset)) %>% 
  dplyr::select(-level) %>% 
  spread(cell_type, coex) %>% 
  session$dataWrangler$setColAsRownames("dataset") %>% 
  session$graphingUtils$heatmap(cluster_rows = FALSE, cluster_cols = FALSE)

fdata$mainCoef %>% 
  filter(pair_id == xdata$pair) %>% filter(!grepl("psdbk", coex_file)) %>% 
  session$graphingUtils$ggplot(aes(x = coex_file, y = coex)) + 
  geom_bar(stat = "identity") + 
  session$graphingUtils$tiltX(angle = 90) + 
  xlim(c("bk_coefmats_rosmap-ihc.rds", 
         "bk_coefmats_rosmap-ihcres.rds", 
         "bk_coefmats_rosmap.rds", 
         "bk_coefmats_rosmap-mgpres.rds", 
         "bk_coefmats_velmeshev.rds", 
         "bk_coefmats_velmeshev-mgpres.rds"))

fdata$main %>% 
  filter(pair_id == xdata$pair) %>% filter(!grepl("psdbk", coex_file), level == "bk")


pdata$ctprofiles %>% filter(gene_id %in% (xdata$pair %>% strsplit("\\.") %>% unlist())) %>% 
  mutate(expr = log2(expr + 1)) %>% 
  group_by(gene_id) %>% mutate(expr = ((expr - mean(expr)) / sd(expr))) %>% ungroup() %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = expr)) + 
  geom_boxplot() + 
  facet_wrap(~gene_id, ncol = 1) +
  session$graphingUtils$tiltX(angle = 90)


pdata$ctprofiles %>% filter(gene_id %in% c("ENSG00000100811")) %>% 
  group_by(gene_id, cell_type) %>% summarize(expr = mean(expr))


fdata$coreg$coregTf <- fdata$coreg$coregTf %>% mutate(pair_id = paste0(target_a, ".", target_b))
fdata$coreg$coregTf <- fdata$coreg$coregTf %>% 
  left_join(pdata$genes %>% dplyr::select(tf = gene_id, tf_gene = gene)) %>% 
  left_join(pdata$genes %>% dplyr::select(target_a = gene_id, target_a_gene = gene)) %>% 
  left_join(pdata$genes %>% dplyr::select(target_b = gene_id, target_b_gene = gene))


x <- readRDS(paste0(workspace$outputDir, "bk_exprmats_rosmap-ihcres.rds"))

(x$ccvModel$ctpMat %>% session$dataWrangler$setRownameAsColumn("sample") %>% 
  left_join(x$ccvModel$exprmats$orig[(xdata$pair %>% strsplit("\\.") %>% unlist()), ] %>% t() %>% session$dataWrangler$setRownameAsColumn("sample")) %>% 
  session$dataWrangler$setColAsRownames("sample")) %>% 
  GGally::ggpairs()


(fdata$coreg$coregTf %>% filter(pair_id == xdata$pair) -> xdata$coreg)
fdata$curation$curation$plos %>% filter(tolower(TF_Symbol_Human) %in% tolower(xdata$coreg$tf_gene), tolower(Target_Symbol_Human) %in% tolower(c(xdata$coreg$target_a_gene, xdata$coreg$target_b_gene))) %>%
  dplyr::select(TF_Symbol_Human, Target_Symbol_Human, PubMed_ID) 
fdata$curation$curation$expand %>% filter(tolower(TF_Gene_Name) %in% tolower(xdata$coreg$tf_gene), tolower(Target_Gene_Name) %in% tolower(c(xdata$coreg$target_a_gene, xdata$coreg$target_b_gene)))


fdata$curation <- read_rds(paste0(workspace$outputDir, "dtri_curation.rds"))




# +++++++++++++++++++++++++++++++++++++===
# let's check the TFs expressions in each dataset

fdata <- list(); gc()

fdata$coreg <- read_rds(paste0(workspace$outputDir, "dtri_coreg.rds"))

fdata$tfs <- fdata$coreg$coregTf %>% group_by(tf) %>% summarize(n_coreg_pairs = n()) %>% arrange(desc(n_coreg_pairs)) %>% 
  left_join(pdata$genes %>% dplyr::select(tf = gene_id, gene)) 
  
fdata$tfs %>%   
  session$graphingUtils$ggplot(aes(x = gene, y = n_coreg_pairs)) + 
  geom_bar(stat = "identity") + 
  session$graphingUtils$tiltX(angle = 90) + 
  xlim(xdata$main$gene)

xdata <- list()

xdata$main <- pdata$ctprofiles %>% filter(gene_id %in% fdata$tfs$tf) %>%
  mutate(expr = log2(expr + 1)) %>% 
  mutate(dataset = paste0(cell_type, ".", dataset)) %>% 
  dplyr::select(-cell_type) %>% 
  spread(dataset, expr) %>%  
  session$dataWrangler$setColAsRownames("gene_id")

xdata$main <- xdata$main[fdata$tfs$tf, ]
rownames(xdata$main) <- fdata$tfs$gene

xdata$main %>% session$graphingUtils$heatmap(cluster_rows = FALSE, cluster_cols = FALSE)



















# +++++++++++++++++++++====
# run it by the consensus networks --- to see if things look even better? --
# *** this is not necessary. 

fdata <- list(); gc()

fdata$coreg <- read_rds(paste0(workspace$outputDir, "dtri_coreg.rds"))
fdata$coregGenes <- fdata$coreg$coregmat %>% rownames()

fdata$coexmatFiles <- pdata$files$coexmats %>% filter(dataset == "consensus")
fdata$coexmatFiles <- fdata$coexmatFiles$outfile %>% session$dataWrangler$attachNames()

fdata$coexmats <- fdata$coexmatFiles %>% lapply(function(file) { read_rds(paste0(workspace$outputDir, file)) })

fdata$main <- fdata$coexmats %>% session$collectionUtils$lapplyWithName(function(coexFile, coexmats) {
  
  coexmats %>% session$collectionUtils$lapplyWithName(function(celltype, coexmat) {
    
    genes <- fdata$coregGenes %>% intersect(rownames(coexmat)) %>% unique() %>% sort()
    
    coexmat <- coexmat[genes, genes]
    coextbl <- tibble(pair_id = workspace$utils$getPairIds(coexmat), coex = workspace$utils$vectorize(coexmat))
    coregtbl <- fdata$coreg$coregTbl %>% filter(pair_id %in% coextbl$pair_id)
    
    auroc <- session$evaluationUtils$auroc((coextbl %>% arrange(desc(coex)))$pair_id, trueSet = (coregtbl %>% filter(coreg == 1))$pair_id)
    
    tibble(level = coexFile,
           cell_type = celltype,
           n_gene = length(genes), 
           n_pairs_full = auroc[2], 
           n_pairs_positive = auroc[1], 
           auc = auroc[3], 
           pvalue = auroc[4])
    
  }, verbose = FALSE) %>% session$dataWrangler$rbind()
  
}) %>% session$dataWrangler$rbind()

fdata$main %>%
  session$graphingUtils$ggplot(aes(x = level, y = auc)) + 
  geom_point() + 
  geom_line(aes(group = cell_type)) + 
  session$graphingUtils$tiltX(angle = 90)













