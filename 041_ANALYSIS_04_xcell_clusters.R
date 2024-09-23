
pdata <- list()

pdata$figuresDir <- paste0(workspace$workspaceDir, "041_ANALYSIS_04_xcell_clusters_FIGURES/")

pdata$celltypes <- c("excitatory", "inhibitory", "opc", "oligodendrocyte", "astrocyte", "microglia")
names(pdata$celltypes) <- pdata$celltypes

pdata$ctprofiles <- read_rds(paste0(workspace$outputDir, "sc_stats_ctprofiles.rds"))

pdata$genes <- readRDS(paste0(workspace$outputDir, "genes_metadata.rds")) # genes

pdata$mkrs <- readRDS(paste0(workspace$outputDir, "umkrs.rds"))

pdata$goTerms <- readRDS(paste0(workspace$outputDir, "gene_ontology.rds"))

pdata$files$coexmats <- tibble(outfile = list.files(workspace$outputDir)) %>% 
  filter(grepl("coexmats", outfile)) %>% 
  mutate(outfile_name = gsub(".rds", "", outfile)) %>% 
  mutate(strs = strsplit(outfile_name, "_")) %>% 
  mutate(level = strs %>% sapply(function(currStr) { currStr[1] })) %>% 
  mutate(type = strs %>% sapply(function(currStr) { currStr[2] })) %>% 
  mutate(dataset = strs %>% sapply(function(currStr) { currStr[3] })) %>% 
  mutate(cell_type = strs %>% sapply(function(currStr) { currStr[4] })) %>% 
  dplyr::select(outfile, level, type, dataset, cell_type)

pdata$files$coefmats <- tibble(outfile = list.files(workspace$outputDir)) %>% 
  filter(grepl("coefmats", outfile)) %>% 
  mutate(outfile_name = gsub(".rds", "", outfile)) %>% 
  mutate(strs = strsplit(outfile_name, "_")) %>% 
  mutate(level = strs %>% sapply(function(currStr) { currStr[1] })) %>% 
  mutate(type = strs %>% sapply(function(currStr) { currStr[2] })) %>% 
  mutate(dataset = strs %>% sapply(function(currStr) { currStr[3] })) %>% 
  mutate(cell_type = strs %>% sapply(function(currStr) { currStr[4] })) %>% 
  dplyr::select(outfile, level, type, dataset, cell_type)

pdata$files$lnkmats <- tibble(outfile = list.files(workspace$outputDir)) %>% 
  filter(grepl("lnkmats", outfile)) %>% 
  mutate(outfile_name = gsub(".rds", "", outfile)) %>% 
  mutate(strs = strsplit(outfile_name, "_")) %>% 
  mutate(level = strs %>% sapply(function(currStr) { currStr[1] })) %>% 
  mutate(type = strs %>% sapply(function(currStr) { currStr[2] })) %>% 
  mutate(cell_type = strs %>% sapply(function(currStr) { currStr[3] })) %>% 
  dplyr::select(outfile, level, type, cell_type)


# +++++++++++++++++++++
# FIRST -- xcell level
# let's build consensus networks; one per cell type, by combining all the datasets

fdata <- list(); gc()

fdata$celltypes <- pdata$celltypes

fdata$genes <- fdata$celltypes %>% session$collectionUtils$lapply(function(currCelltype) {
  pdata$ctprofiles %>% 
    filter(cell_type == currCelltype) %>% 
    dplyr::select(dataset, gene_id) %>%
    mutate(value = 1) %>% 
    spread(dataset, value) %>% 
    na.omit() %>% 
    session$dataWrangler$extractColumn("gene_id") %>% 
    unique() %>% 
    sort()
})

fdata$coexmats <- fdata$genes %>% lapply(function(currGenes) {
  coexmat <- matrix(nrow = length(currGenes), ncol = length(currGenes))
  rownames(coexmat) <- currGenes
  colnames(coexmat) <- currGenes
  return(coexmat)
})

fdata$coexmats <- names(fdata$coexmats) %>% mclapply(function(currCelltype) {
  
  coexmatClone <- fdata$coexmats[[currCelltype]]
  
  lowerCoords <- which(lower.tri(coexmatClone), arr.ind = TRUE) # use the lower tri indices as primary
  upperCoords <- lowerCoords[, c("col", "row")] # figure out the matching coordinates to "flip it back"
  colnames(upperCoords) <- c("row", "col")
  
  print("Computing individual co-expression networks.")
  
  coexRanks <- rep(0, nrow(lowerCoords))
  
  files <- pdata$files$coexmats %>% filter(cell_type == currCelltype, level == "sc")
  
  files$outfile %>% session$collectionUtils$foreach(function(currFile) {
    
    currDCoexmat <- read_rds(paste0(workspace$outputDir, currFile))
    currDCoexmat <- currDCoexmat[rownames(coexmatClone), colnames(coexmatClone)]
    currDCoexRanks <- currDCoexmat[lowerCoords]
    currDCoexRanks[is.na(currDCoexRanks)] <- currDCoexRanks %>% median(na.rm = TRUE)
    
    # normalize into ranks overall again; this would've been off because of data subsetting
    currDCoexRanks <- currDCoexRanks %>% rank(ties.method = "average")
    currDCoexRanks <- currDCoexRanks / max(currDCoexRanks)
    
    coexRanks <<- coexRanks + currDCoexRanks
  })
  
  print("Normalizing co-expressions.")
  
  coexRanks <- coexRanks %>% rank(ties.method = "average")
  coexRanks <- coexRanks / max(coexRanks)
  
  print("Converting to matrix form.")
  
  coexmatClone[lowerCoords] <- coexRanks
  coexmatClone[upperCoords] <- coexRanks
  diag(coexmatClone) <- 1
  
  return(coexmatClone)
}, mc.cores = 6)

names(fdata$coexmats) <- fdata$celltypes

fdata$coexmats %>% saveRDS(paste0(workspace$outputDir, paste0("sc_coexmats_consensus.rds")))
# fdata$coexmats <- read_rds(paste0(workspace$outputDir, paste0("sc_coexmats_consensus.rds")))

# +++++++++++++++++++++
# let's cluster the consensus stuff

fdata <- list(); gc()

fdata$coexmats <- read_rds(paste0(workspace$outputDir, paste0("sc_coexmats_consensus.rds")))

# OK, let's use this space to test out all the different methods

fdata$clusters <- fdata$coexmats %>% mclapply(function(currCoexmat) {
  workspace$utils$computeClusters(currCoexmat, function(dendro, distMat) {
    dynamicTreeCut::cutreeDynamic(dendro, distM = distMat, cutHeight = 0.2, minClusterSize = 30)
  })
}, mc.cores = length(fdata$coexmats))

xdata <- list()

xdata$celltype <- "opc"

plotDendroAndColors(fdata$clusters$clusterObjs[[xdata$celltype]]$dendro, 
                    colors = labels2colors(fdata$clusters$clusterObjs[[xdata$celltype]]$clusters), 
                    rowText = fdata$clusters$clusterObjs[[xdata$celltype]]$clusters, 
                    dendroLabels = FALSE)

# ++++++++++++==== COMMIT
fdata$clusters %>% saveRDS(paste0(workspace$outputDir, paste0("sc_clusters_objects.rds")))
# =================================


# +++++++++++++++++
# characterize the number of clusters and their n_gene distributions across cell types

fdata <- list(); gc()

fdata$clusters <- read_rds(paste0(workspace$outputDir, paste0("sc_clusters_objects.rds")))

# look at the number of genes per cluster, per cell type
fdata$nGenes <- fdata$clusters %>% session$collectionUtils$lapplyWithName(function(currCelltype, currCluster) {
  currCluster$clusterTbl %>% group_by(cluster) %>% summarize(n_gene = n()) %>% mutate(cell_type = currCelltype)
}) %>% session$dataWrangler$rbind()

xdata <- list()

xdata$celltype <- "oligodendrocyte"

xdata$nGenes <- fdata$nGenes %>% filter(cell_type == xdata$celltype)

xdata$nGenes %>% 
  session$graphingUtils$ggplot(aes(x = n_gene, y = cluster)) + 
  geom_text(aes(label = n_gene), hjust = -0.3, size = 4.5) +
  geom_bar(stat = "identity") + 
  xlim(0, max(xdata$nGenes$n_gene) + 10) + 
  ggtitle(paste0("Number of genes per cluster in ", xdata$celltype)) 

# flatten a big list of clusters, attach the cell type names to the clusters (cell type of origin)

fdata$clustersFlat <- fdata$clusters %>% session$collectionUtils$lapplyWithName(function(currCelltype, currCluster) {
  currFlat <- currCluster$clusterFlat[names(currCluster$clusterFlat) != "id_0000"]
  names(currFlat) <- paste0(currCelltype, ".", names(currFlat))
  return(currFlat)
})

fdata$clustersMats <- fdata$clusters %>% session$collectionUtils$lapplyWithName(function(currCelltype, currCluster) {
  currMat <- currCluster$clusterMat[, colnames(currCluster$clusterMat) != "id_0000", drop = FALSE]
  colnames(currMat) <- paste0(currCelltype, ".", colnames(currMat))
  return(currMat)
})

# ++++++++++++==== COMMIT
list(clusterObjs = fdata$clusters, clusterFlats = fdata$clustersFlat, clusterMats = fdata$clustersMats) %>%
  saveRDS(paste0(workspace$outputDir, paste0("sc_clusters_objects.rds")))

# ++++++++++++++++++++++++++==
# compute the cluster quality in consensus SC to ensure quality and cell type specifciity

fdata <- list(); gc()

fdata$clusters <- read_rds(paste0(workspace$outputDir, paste0("sc_clusters_objects.rds")))
fdata$clustersGeneric <- read_rds(paste0(workspace$outputDir, paste0("sc_generic_clusters.rds")))

fdata$clusterMats <- fdata$clusters$clusterMats %>% session$collectionUtils$lapply(function(currClusterMat) {
  
  genericMat <- fdata$clustersGeneric$clusterMat
  genericMat <- genericMat[intersect(rownames(currClusterMat), rownames(genericMat)), ]
  
  genesDiff <- rownames(currClusterMat) %>% setdiff(rownames(genericMat))
  genericMat[genesDiff, ] <- 0
  genericMat <- genericMat[rownames(currClusterMat), ]
  
  return(currClusterMat %>% cbind(genericMat))
})

fdata$clusterFlats <- fdata$clusters$clusterFlats %>% session$collectionUtils$lapply(function(currClusterFlat) {
  genesValid <- currClusterFlat %>% unlist() %>% unname() %>% unique()
  genericFlat <- fdata$clustersGeneric$clusterFlat %>% lapply(function(genes) { genes %>% intersect(genesValid) })
  return(c(currClusterFlat, genericFlat))
})
  
# fdata$coexmatFiles <- pdata$files$coexmats %>% filter(dataset != "consensus")
# fdata$coexmatFiles <- paste0(workspace$outputDir, fdata$coexmatFiles$outfile) %>% session$dataWrangler$attachNames(fdata$coexmatFiles$outfile)

fdata$coexmats <- read_rds(paste0(workspace$outputDir, paste0("sc_coexmats_consensus.rds")))

fdata$egad <- fdata$clusterMats %>% session$collectionUtils$lapply(function(currClusterMat) {
  workspace$utils$computeEgadForEach(fdata$coexmats, currClusterMat) %>% dplyr::select(coex_ct = coexfile, everything())
}) %>% session$dataWrangler$rbind()

fdata$egad %>% session$graphingUtils$ggplot(aes(x = cluster, y = auc)) + geom_point(aes(color = coex_ct)) + session$graphingUtils$tiltX(angle = 90)

fdata$avgcoex <- fdata$clusterFlats %>% session$collectionUtils$lapplyWithName(function(cellType, currClusterFlat) {
  workspace$utils$computeAvgCoexForEach(fdata$coexmats, currClusterFlat) %>% dplyr::select(coex_ct = coexfile, everything())
}) %>% session$dataWrangler$rbind()

fdata$avgcoex %>% 
  session$graphingUtils$ggplot(aes(x = cluster, y = coex)) + geom_point(aes(color = coex_ct)) + session$graphingUtils$tiltX(angle = 90)


# ========== COMMIT
list(egad = fdata$egad, avgcoex = fdata$avgcoex) %>% 
  saveRDS(paste0(workspace$outputDir, paste0("sc_clusters_integrity_across_celltype.rds")))
# =============


# ++++++++++++++++++++++++++==
# compute the quality of the clusters that were discovered at the sc level
# 1. egad, 2. average co-expression
# in sc, sbj, bk consensus and component networks in the same cell type

fdata <- list(); gc()

fdata$clusters <- read_rds(paste0(workspace$outputDir, paste0("sc_clusters_objects.rds")))
# fdata$clustersGeneric <- read_rds(paste0(workspace$outputDir, paste0("sc_generic_clusters.rds")))

fdata$clusterMats <- fdata$clusters$clusterMats %>% session$collectionUtils$lapply(function(currClusterMat) {
  
  genericMat <- fdata$clustersGeneric$clusterMat
  genericMat <- genericMat[intersect(rownames(currClusterMat), rownames(genericMat)), ]
  
  genesDiff <- rownames(currClusterMat) %>% setdiff(rownames(genericMat))
  genericMat[genesDiff, ] <- 0
  genericMat <- genericMat[rownames(currClusterMat), ]
  
  return(currClusterMat %>% cbind(genericMat))
})

fdata$clusterFlats <- fdata$clusters$clusterFlats %>% session$collectionUtils$lapply(function(currClusterFlat) {
  genesValid <- currClusterFlat %>% unlist() %>% unname() %>% unique()
  genericFlat <- fdata$clustersGeneric$clusterFlat %>% lapply(function(genes) { genes %>% intersect(genesValid) })
  return(c(currClusterFlat, genericFlat))
})

fdata$coexmatFiles <- pdata$files$coexmats %>% filter(dataset != "consensus")
# fdata$coexmatFiles <- paste0(workspace$outputDir, fdata$coexmatFiles$outfile) %>% session$dataWrangler$attachNames(fdata$coexmatFiles$outfile)

fdata$coexmatsConsensus <- list(sc = read_rds(paste0(workspace$outputDir, paste0("sc_coexmats_consensus.rds"))), 
                                sbj = read_rds(paste0(workspace$outputDir, paste0("sbj_coexmats_consensus.rds")))) %>% 
  unlist(recursive = FALSE)

# EGAD

fdata$egadComponents <- fdata$clusterMats %>% session$collectionUtils$lapplyWithName(function(cellType, currClusterMat) {
  
  coexmatFiles <- fdata$coexmatFiles %>% filter(cell_type == cellType | level == "bk")
  coexmats <- paste0(workspace$outputDir, coexmatFiles$outfile) %>% session$dataWrangler$attachNames(coexmatFiles$outfile)
  
  workspace$utils$computeEgadForEach(coexmats, currClusterMat) %>% 
    left_join(coexmatFiles %>% dplyr::select(coexfile = outfile, level, dataset, coex_ct = cell_type), by = "coexfile") %>% 
    mutate(cluster_ct = cellType)
  
}) %>% session$dataWrangler$rbind()


fdata$egadConsensus <- fdata$clusterMats %>% session$collectionUtils$lapplyWithName(function(cellType, currClusterMat) {
  
  coexmats <- fdata$coexmatsConsensus[grepl(cellType, names(fdata$coexmatsConsensus))]
  
  workspace$utils$computeEgadForEach(coexmats, currClusterMat) %>% 
    mutate(str = strsplit(coexfile, "\\.")) %>% 
    mutate(level = str %>% sapply(function(currStr) { currStr[1] }), 
           cell_type = str %>% sapply(function(currStr) { currStr[2] })) %>% 
    dplyr::select(-coexfile, -str)
  
}) %>% session$dataWrangler$rbind()

fdata$egad <- fdata$egadComponents %>% dplyr::select(cluster, auc, avg_node_degree, degree_null_auc, n_gene, level, dataset, coex_ct, cluster_ct) %>% 
  rbind(fdata$egadConsensus %>% 
          mutate(dataset = "consensus") %>% 
          dplyr::select(cluster, auc, avg_node_degree, degree_null_auc, n_gene, level, dataset, coex_ct = cell_type) %>% 
          mutate(cluster_ct = coex_ct))


# AVG-COEX

fdata$avgcoexComponents <- fdata$clusterFlats %>% session$collectionUtils$lapplyWithName(function(cellType, currClusterFlat) {
  
  coexmatFiles <- fdata$coexmatFiles %>% filter(cell_type == cellType | level == "bk")
  coexmats <- paste0(workspace$outputDir, coexmatFiles$outfile) %>% session$dataWrangler$attachNames(coexmatFiles$outfile)
  
  workspace$utils$computeAvgCoexForEach(coexmats, currClusterFlat) %>% 
    left_join(coexmatFiles %>% dplyr::select(coexfile = outfile, level, dataset, coex_ct = cell_type), by = "coexfile") %>% 
    mutate(cluster_ct = cellType)
  
}) %>% session$dataWrangler$rbind()


fdata$avgcoexConsensus <- fdata$clusterFlats %>% session$collectionUtils$lapplyWithName(function(cellType, currClusterFlat) {
  
  coexmats <- fdata$coexmatsConsensus[grepl(cellType, names(fdata$coexmatsConsensus))]
  
  workspace$utils$computeAvgCoexForEach(coexmats, currClusterFlat) %>% 
    mutate(str = strsplit(coexfile, "\\.")) %>% 
    mutate(level = str %>% sapply(function(currStr) { currStr[1] }), 
           cell_type = str %>% sapply(function(currStr) { currStr[2] })) %>% 
    dplyr::select(-coexfile, -str)
  
}) %>% session$dataWrangler$rbind()

fdata$avgcoex <- fdata$avgcoexComponents %>% dplyr::select(cluster, coex, n_gene, n_lnks, level, dataset, coex_ct, cluster_ct) %>% 
  rbind(fdata$avgcoexConsensus %>% 
          mutate(dataset = "consensus") %>% 
          dplyr::select(cluster, coex, n_gene, n_lnks, level, dataset, coex_ct = cell_type) %>% 
          mutate(cluster_ct = coex_ct))


# ========== COMMIT
list(egad = fdata$egad, avgCoex = fdata$avgcoex) %>% 
  saveRDS(paste0(workspace$outputDir, paste0("sc_clusters_integrity_within_celltype.rds")))
# =============

# for the bulk comparisons, before and after correction, etc, need to use the raw Pearson's correlations instead
# for all the bulk networks, compute cluster integrity within cell types

fdata <- list(); gc()

fdata$clusters <- read_rds(paste0(workspace$outputDir, paste0("sc_clusters_objects.rds")))
fdata$clustersGeneric <- read_rds(paste0(workspace$outputDir, paste0("sc_generic_clusters.rds")))
fdata$clusterFlats <- c(fdata$clusters$clusterFlats, list(generic = fdata$clustersGeneric$clusterFlat))

fdata$coexmatFiles <- pdata$files$coefmats %>% filter(level == "bk")

fdata$avgcoex <- fdata$clusterFlats %>% session$collectionUtils$lapplyWithName(function(cellType, currClusterFlat) {
  
  coexmatFiles <- fdata$coexmatFiles
  coexmats <- paste0(workspace$outputDir, coexmatFiles$outfile) %>% session$dataWrangler$attachNames(coexmatFiles$outfile)
  
  workspace$utils$computeAvgCoexForEach(coexmats, currClusterFlat) %>% 
    left_join(coexmatFiles %>% dplyr::select(coexfile = outfile, level, dataset, coex_ct = cell_type), by = "coexfile") %>% 
    mutate(cluster_ct = cellType)
  
}) %>% session$dataWrangler$rbind()

# ========== COMMIT
fdata$avgcoex %>% 
  saveRDS(paste0(workspace$outputDir, paste0("sc_clusters_integrity_in_bk_coefmats.rds")))
# =============

# =============================
# ++++++++++++++++++++++++++==
# OK - finally can start plotting some figures & do some writing
# first, characterize the number of clusters / cluster sizes at the single cell level

fdata <- list(); gc()

fdata$clusters <- read_rds(paste0(workspace$outputDir, paste0("sc_clusters_objects.rds")))

fdata$objs <- fdata$clusters$clusterObjs

fdata$clusterFlat <- fdata$clusters$clusterFlats %>% unname() %>% unlist(recursive = FALSE)

fdata$clusterSize <- fdata$clusterFlat %>% session$collectionUtils$lapplyWithName(function(clusterName, genes) {
  tibble(cluster = clusterName, n_gene = length(genes))
}) %>% session$dataWrangler$rbind()

fdata$clusterSize <- fdata$clusterSize %>%
  mutate(cell_type = cluster %>% sapply(function(currStr) { currStr %>% str_extract("^[a-z]+") }))

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "sfigure_01.eps"

xdata$main <- fdata$clusterSize %>% mutate(cell_type = workspace$utils$fmtCelltypes(cell_type))

xdata$means <- fdata$clusterSize$n_gene %>% mean()

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = n_gene)) + 
  geom_histogram(aes(fill = cell_type)) + 
  geom_vline(xintercept = xdata$means, linetype = "dashed") +
  scale_x_continuous(trans = "log10") + 
  scale_fill_brewer(palette = "Set1") +
  ggtitle("Distribution of cluster size") + 
  xlab("Number of genes") + 
  ylab("Number of clusters") + 
  labs(fill = "Cell type")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 25, height = 18, dpi = 1000) 
# =================================================


# =================================================
# OUTPUT FIGURE ------- number of clusters per cell type
xdata <- list()

xdata$figName <- "sfigure_02.eps"

xdata$main <- fdata$clusterSize %>% group_by(cell_type) %>% summarize(n_cluster = n()) %>% 
  mutate(cell_type = workspace$utils$fmtCelltypes(cell_type))

xdata$means <- fdata$clusterSize$n_gene %>% mean()

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = n_cluster)) + 
  geom_bar(stat = "identity", aes(fill = cell_type)) + 
  scale_fill_brewer(palette = "Set1") +
  session$graphingUtils$tiltX(angle = 90) + 
  geom_text(aes(label = n_cluster), vjust = -1, size = 7) + 
  ylim(0, max(xdata$main$n_cluster) + 5) + 
  ggtitle("Number of clusters per cell type") + 
  theme(legend.position = "none") + 
  xlab("Cell type") + 
  ylab("Number of clusters")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 20, height = 18, dpi = 1000) 
# =================================================


# gene, across all networks

xdata <- list()

xdata$main <- fdata$clusters$clusterObjs %>% session$collectionUtils$lapplyWithName(function(currCelltype, currCluster) {
  currCluster$clusterTbl %>% mutate(in_cluster = cluster_int > 0) %>% mutate(cell_type = currCelltype)
}) %>% session$dataWrangler$rbind()

xdata$inGenes <- xdata$main %>% filter(in_cluster) %>% session$dataWrangler$extractColumn("gene_id") %>% unique()
xdata$inGenes %>% length() # number of genes in cluster

xdata$allGenes <- xdata$main$gene_id %>% unique()
xdata$allGenes %>% length()

#=============
# gene coverage
xdata <- list()

xdata$genesOut <- fdata$clusters$clusterObjs %>% sapply(function(currCluster) {
  currCluster$clusterTbl %>% filter(cluster_int == 0) %>% nrow()
})

xdata$genesIn <- fdata$clusterSize %>% group_by(cell_type) %>% summarize(n_gene_in = sum(n_gene))

xdata$genesIn$n_gene_out <- xdata$genesOut[xdata$genesIn$cell_type]

xdata$genesIn <- xdata$genesIn %>% mutate(total = n_gene_in + n_gene_out) %>% mutate(frac_in = n_gene_in / total)

xdata$main <- xdata$genesIn %>% dplyr::select(cell_type, n_gene_in, n_gene_out) %>% 
  gather(status, n_gene, -cell_type) %>% mutate(status = factor(status, levels = c("n_gene_out", "n_gene_in")))

# =================================================
# OUTPUT FIGURE ------- # gene coverage --- computed above
# xdata <- list()

xdata$figName <- "sfigure_03.eps"

xdata$plotMain <- xdata$main %>% mutate(cell_type = workspace$utils$fmtCelltypes(cell_type))

xdata$plotMain <- xdata$plotMain %>% mutate(Status = status %>% sapply(function(str) { 
  if (str == "n_gene_in") { "In" } else { "Out" }
}))

xdata$plot <- xdata$plotMain %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = n_gene)) + 
  geom_bar(aes(fill = Status), stat = "identity") +
  session$graphingUtils$tiltX(angle = 90) + 
  ggtitle("Number of genes assigned to a cluster") + 
  xlab("Cell type") + 
  ylab("Number of genes") + 
  scale_fill_manual(values = c("grey60", "black"))

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 25, height = 20, dpi = 1000) 
# =================================================

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "sfigure_04.eps"

xdata$main <- fdata$clusterSize %>% mutate(cell_type = workspace$utils$fmtCelltypes(cell_type))

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(y = cluster, x = n_gene), size = "SMALL") + 
  geom_bar(aes(fill = cell_type), stat = "identity") + 
  scale_fill_brewer(palette = "Set1") + 
  labs(fill = "Cell type", x = "Number of genes", y = "Cluster")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 20, height = 30, dpi = 1000) 
# =================================================


# =================================================
# OUTPUT FIGURE

# gene intersections across the clusters

xdata <- list()

xdata$main <- fdata$clusters$clusterFlats %>% unname() %>% unlist(recursive = FALSE)

xdata$main <- xdata$main %>% sapply(function(cluster1) {
  xdata$main %>% sapply(function(cluster2) {
    (cluster1 %>% intersect(cluster2) %>% length()) / (cluster1 %>% union(cluster2) %>% length())
  })
})

xdata$figNames <- c("sfigure_05.png", "sfigure_06.png", "sfigure_07.png", "sfigure_08.png", "sfigure_09.png", "sfigure_10.png")

xdata$i <- 6
(xdata$celltype <- sort(names(fdata$clusters$clusterFlats))[xdata$i])
xdata$clusterNames <- fdata$clusters$clusterFlats[[xdata$celltype]] %>% names()
xdata$clusterNamesOthrs <- (colnames(xdata$main))[!(colnames(xdata$main) %in% xdata$clusterNames)]
png(paste0(pdata$figuresDir, xdata$figNames[xdata$i]), units = "cm", 
    width = 55, height = 20, res = 400)
xdata$main[xdata$clusterNames, xdata$clusterNamesOthrs] %>% 
  session$graphingUtils$heatmap(size = session$graphingUtils$MEDIUM, cluster_rows = FALSE, cluster_cols = FALSE, 
                                main = workspace$utils$fmtCelltypes(xdata$celltype))
dev.off()

# =================================================
# =================================================

# =================================================
# OUTPUT FIGURE
# dendrograms

xdata <- list()

xdata$i <- 6

xdata$figNames <- paste0("sfigure_", 11:16)

(xdata$celltype <- sort(pdata$celltypes)[xdata$i])

png(paste0(pdata$figuresDir, paste0(xdata$figNames[xdata$i], "_", xdata$celltype, ".png")), units = "cm", width = 30, height = 20, res = 400)

plotDendroAndColors(fdata$objs[[xdata$celltyp]]$dendro, 
                    colors = labels2colors(fdata$objs[[xdata$celltype]]$clusters), 
                    rowText = fdata$objs[[xdata$celltype]]$clusters, 
                    dendroLabels = FALSE, abHeight = 0.2) %>% print()

dev.off()
# =================================================




# ++++++++++++++++++++++++++==
# here - look at how well single cell networks capture single cell clusters

fdata <- list(); gc()

fdata$main <- readRDS(paste0(workspace$outputDir, paste0("sc_clusters_integrity_within_celltype.rds")))
fdata$main <- fdata$main %>% lapply(function(currTbl) { currTbl %>% filter(level == "sc") })


# avg co-expression
xdata <- list()

xdata$main <- fdata$main$avgCoex

xdata$means <- xdata$main %>% group_by(coex_ct, dataset) %>% summarize(coex = mean(coex))

xdata$main %>% 
  session$graphingUtils$ggplot(aes(y = coex_ct, x = coex, color = coex_ct)) + 
  geom_point(shape = 1) + 
  geom_point(data = xdata$means, size = 5) +
  session$graphingUtils$tiltX(angle = 90) + 
  ggtitle("Average consensus co-expression per cluster by dataset and cell type") +
  facet_wrap(~dataset, ncol = 1) +
  xlim(0, 1)


# egad
xdata <- list()

xdata$main <- fdata$main$egad

xdata$means <- xdata$main %>% group_by(coex_ct, dataset) %>% summarize(auc = mean(auc))

xdata$main %>% 
  session$graphingUtils$ggplot(aes(y = coex_ct, x = auc, color = coex_ct)) + 
  geom_point(shape = 1) + 
  geom_point(data = xdata$means, size = 5) +
  session$graphingUtils$tiltX(angle = 90) + 
  ggtitle("EGAD performance per cluster by dataset and cell type") +
  facet_wrap(~dataset, ncol = 1) +
  xlim(0, 1)


# OK starting with average co-expressions

xdata <- list()

(xdata$celltype <- pdata$celltypes[5])

xdata$main <- fdata$main$avgCoex %>% filter(coex_ct == xdata$celltype)

xdata$consensus <- xdata$main %>% filter(dataset == "consensus")
xdata$component <- xdata$main %>% filter(dataset != "consensus")

xdata$componentMeans <- xdata$component %>% group_by(cluster, coex_ct) %>% summarize(coex = mean(coex))

xdata$component %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_density(aes(color = dataset)) +
  geom_density(data = xdata$consensus, size = 1) +
  xlim(0, 1) + 
  ggtitle(paste0("Average co-expression in component networks")) 

xdata$component %>%
  session$graphingUtils$ggplot(aes(x = cluster, y = coex, color = dataset)) + 
  geom_point(color = "grey70") + 
  geom_line(aes(group = dataset), color = "grey70") + 
  geom_point(data = xdata$componentMeans, color = "black", size = 3) + 
  geom_line(data = xdata$componentMeans, group = 1, color = "black", size = 1) +
  geom_point(data = xdata$consensus, color = "red", size = 3) + 
  geom_line(data = xdata$consensus, group = 1, color = "red", size = 1) +
  session$graphingUtils$tiltX(angle = 90) + 
  ggtitle(paste0("Average co-expression in component networks: ", xdata$celltype)) + 
  ylim(0, 1)

xdata$component %>%
  dplyr::select(cluster, dataset, coex) %>% 
  spread(dataset, coex) %>% 
  session$dataWrangler$setColAsRownames("cluster") %>% 
  session$graphingUtils$heatmap(cluster_row = FALSE, cluster_col = FALSE)

# OK, moving on to using EGAD

xdata <- list()

(xdata$celltype <- pdata$celltypes[6])

xdata$main <- fdata$main$egad #%>% filter(coex_ct == xdata$celltype)

xdata$consensus <- xdata$main %>% filter(dataset == "consensus")
xdata$component <- xdata$main %>% filter(dataset != "consensus")
xdata$componentMeans <- xdata$component %>% group_by(cluster) %>% summarize(auc = mean(auc))

xdata$component %>% 
  session$graphingUtils$ggplot(aes(x = auc)) + 
  geom_density(aes(color = dataset)) +
  geom_density(data = xdata$consensus, size = 1) +
  xlim(0, 1) + 
  ggtitle(paste0("Cluster integrity (EGAD-AUROC) in component networks")) 

xdata$component %>%
  session$graphingUtils$ggplot(aes(x = cluster, y = auc, color = dataset)) + 
  geom_point(color = "grey70") + 
  geom_point(data = xdata$componentMeans, color = "black", size = 2) + 
  geom_line(data = xdata$componentMeans, group = 1, color = "black") +
  session$graphingUtils$tiltX(angle = 90) + 
  ggtitle(paste0("EGAD performance in component networks: ", xdata$celltype)) + 
  geom_hline(yintercept = 0.5, linetype = "dashed")

xdata$component %>%
  dplyr::select(cluster, dataset, auc) %>% 
  spread(dataset, auc) %>% 
  session$dataWrangler$setColAsRownames("cluster") %>% 
  session$graphingUtils$heatmap(cluster_row = FALSE, cluster_col = FALSE)




# ++++++++++++++++++++++++++++++===========
# HERE - plot the EGAD-AUROC performance at the three levels, for all datasets;
# maybe split by cell types - rank by highest min xBulk performance between ROSMAP and Velmeshev
# later add for each cluster 
# -- correlations in cell type profile (rank?)
# -- inter-cell type synchrony
# -- number of genes / fraction of which are ribosomal protein genes

fdata <- list(); gc()

fdata$main <- readRDS(paste0(workspace$outputDir, paste0("sc_clusters_integrity_within_celltype.rds")))$egad

fdata$main <- fdata$main %>% filter(dataset != "consensus", !grepl("generic", cluster)) 

# lets use the bulk level to sort
fdata$xbulk <- fdata$main %>% filter(level == "bk", dataset %in% c("rosmap", "velmeshev")) %>% 
  dplyr::select(cluster, dataset, auc) %>% 
  spread(dataset, auc)

fdata$xbulk <- fdata$xbulk %>% mutate(max_xbulk = pmax(rosmap, velmeshev))

fdata$xbulk <- fdata$xbulk %>% arrange(desc(max_xbulk)) # these are the clusters to focus on

# SOFT COMMIT #===============
fdata$xbulk %>% saveRDS(paste0(workspace$outputDir, "sc_clusters_xbulk_ranking.rds"))
pdata$clusters <- readRDS(paste0(workspace$outputDir, "sc_clusters_xbulk_ranking.rds"))
# ============================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_07.eps" # xcell performance 

xdata$main <- fdata$main %>% filter(level == "sc")

xdata$main <- xdata$main %>% mutate(dataset = dataset %>% workspace$utils$fmtDataset())

xdata$main <- xdata$main %>% # map na values
  dplyr::select(cluster, dataset, auc) %>% spread(dataset, auc) %>% 
  gather(dataset, auroc, -cluster) 

xdata$plot <- xdata$main %>%
  session$graphingUtils$ggplot(aes(x = cluster, y = dataset)) + 
  geom_tile(aes(fill = auroc)) + 
  scale_fill_gradient(low = "white", high = "red", na.value = "grey60", limits = c(0.5, 1)) + 
  session$graphingUtils$tiltX(angle = 90) +
  scale_x_discrete(limits = c(pdata$clusters$cluster[1:20], pdata$clusters$cluster[96:115])) +
  ggtitle("xCell")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 35, height = 13) 

# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_08.eps" # xsbj performance 

xdata$main <- fdata$main %>% filter(level == "sbj")

xdata$main <- xdata$main %>% mutate(dataset = dataset %>% workspace$utils$fmtDataset())

xdata$main <- xdata$main %>% # map na values
  dplyr::select(cluster, dataset, auc) %>% spread(dataset, auc) %>% 
  gather(dataset, auroc, -cluster) 

xdata$plot <- xdata$main %>%
  session$graphingUtils$ggplot(aes(x = cluster, y = dataset)) + 
  geom_tile(aes(fill = auroc)) + 
  scale_fill_gradient(low = "white", high = "red", na.value = "grey60", limits = c(0.5, 1)) + 
  session$graphingUtils$tiltX(angle = 90)  +
  xlim(c(pdata$clusters$cluster[1:20], pdata$clusters$cluster[96:115])) + 
  ggtitle("xSubject")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 35, height = 13) 

# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_09.eps" # xbulk performance 

xdata$main <- fdata$main %>% filter(level == "bk", dataset %in% c("rosmap", "velmeshev"))

xdata$main <- xdata$main %>% mutate(dataset = dataset %>% workspace$utils$fmtDataset())

xdata$main <- xdata$main %>% # map na values
  dplyr::select(cluster, dataset, auc) %>% spread(dataset, auc) %>% 
  gather(dataset, auroc, -cluster) 

xdata$plot <- xdata$main %>%
  session$graphingUtils$ggplot(aes(x = cluster, y = dataset)) + 
  geom_tile(aes(fill = auroc)) + 
  scale_fill_gradient(low = "white", high = "red", na.value = "grey60", limits = c(0.4, 1)) + 
  session$graphingUtils$tiltX(angle = 90)  +
  xlim(c(pdata$clusters$cluster[1:20], pdata$clusters$cluster[96:115])) + 
  ggtitle("xBulk")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 35, height = 10) 

# =================================================

# +++++++++++++++++++++++++++++++++============================
# get the gene sizes + ribosomal protein gene counts & plot it out

fdata <- list(); gc()

fdata$clusters <- read_rds(paste0(workspace$outputDir, paste0("sc_clusters_objects.rds")))
fdata$clusters <- fdata$clusters$clusterFlats %>% unname() %>% unlist(recursive = FALSE)

pdata$rpgs <- pdata$genes %>% filter(grepl("^RPL|^RPS.*$", gene)) %>% filter(!grepl("pseudogene", gene_name)) 

fdata$nGenes <- fdata$clusters %>% session$collectionUtils$lapplyWithName(function(id, clusterGenes) {
  nRpg <- clusterGenes %>% intersect(pdata$rpgs$gene_id) %>% length()
  nGene <- clusterGenes %>% length()
  tibble(cluster = id, n_gene = nGene, n_rpg = nRpg)
}) %>% session$dataWrangler$rbind()

fdata$nGenes %>% arrange(desc(n_rpg))

fdata$nGenes %>% 
  session$graphingUtils$ggplot(aes(x = cluster, y = n_gene)) + 
  geom_bar(stat = "identity") + 
  session$graphingUtils$tiltX(angle = 90) + 
  xlim(c(pdata$clusters$cluster[1:20], pdata$clusters$cluster[96:115]))

fdata$nGenes %>% 
  session$graphingUtils$ggplot(aes(x = cluster, y = n_rpg)) + 
  geom_bar(stat = "identity") + 
  session$graphingUtils$tiltX(angle = 90) + 
  xlim(c(pdata$clusters$cluster[1:20], pdata$clusters$cluster[96:115]))


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_10.eps"

xdata$plot <- fdata$nGenes %>% 
  session$graphingUtils$ggplot(aes(x = cluster)) + 
  geom_bar(aes(y = n_gene), stat = "identity", fill = "grey50") +
  geom_bar(aes(y = n_rpg), stat = "identity", fill = "blue") + 
  session$graphingUtils$tiltX(angle = 90) + 
  scale_x_discrete(limits = c(pdata$clusters$cluster[1:20], pdata$clusters$cluster[96:115]), position = "top") +
  ylab("Number of genes") +
  ggtitle("Number of genes by cluster")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 35, height = 15) 

# =================================================


# =================================================
# ++++++++++++++++++++==============================
# let's see how the clusters change after CCV correction

fdata <- list(); gc()

fdata$main <- read_rds(paste0(workspace$outputDir, "sc_clusters_integrity_in_bk_coefmats.rds"))

fdata$main %>% filter(dataset %in% c("velmeshev", "velmeshev-mgpres")) %>% 
  session$graphingUtils$ggplot(aes(y = cluster, x = coex, color = dataset)) + 
  geom_point(size = 4) +
  ylim(c(pdata$clusters$cluster[1:20], pdata$clusters$cluster[96:115]))

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_11.eps"

xdata$plot <- fdata$main %>% 
  filter(dataset %in% c("rosmap-ihc", "rosmap-ihcres")) %>%
  dplyr::select(cluster, dataset, coex) %>% 
  spread(dataset, coex) %>% 
  mutate(diff = `rosmap-ihcres` - `rosmap-ihc`) %>% 
  session$graphingUtils$ggplot(aes(x = cluster, y = diff)) + 
  geom_bar(stat = "identity") +
  xlim(c(pdata$clusters$cluster[1:20], pdata$clusters$cluster[96:115])) + 
  session$graphingUtils$tiltX(angle = 90) + 
  ggtitle("ROSMAP IHC") 

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 35, height = 11) 

# =================================================

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_12.eps"

xdata$plot <- fdata$main %>% 
  filter(dataset %in% c("rosmap", "rosmap-mgpres")) %>% 
  dplyr::select(cluster, dataset, coex) %>% 
  spread(dataset, coex) %>% 
  mutate(diff = `rosmap-mgpres` - `rosmap`) %>% 
  session$graphingUtils$ggplot(aes(x = cluster, y = diff)) + 
  geom_bar(stat = "identity") +
  xlim(c(pdata$clusters$cluster[1:20], pdata$clusters$cluster[96:115])) + 
  session$graphingUtils$tiltX(angle = 90) +
  ggtitle("ROSMAP MGP")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 35, height = 11) 

# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_13.eps"

xdata$plot <- fdata$main %>% 
  filter(dataset %in% c("velmeshev", "velmeshev-mgpres")) %>%
  dplyr::select(cluster, dataset, coex) %>% 
  spread(dataset, coex) %>% 
  mutate(diff = `velmeshev-mgpres` - `velmeshev`) %>% 
  session$graphingUtils$ggplot(aes(x = cluster, y = diff)) + 
  geom_bar(stat = "identity") +
  xlim(c(pdata$clusters$cluster[1:20], pdata$clusters$cluster[96:115])) + 
  session$graphingUtils$tiltX(angle = 90) +
  ggtitle("Velmeshev MGP")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 35, height = 11) 

# =================================================

# ++++++++++++++++=============================
# display cell type profiles for the clusters

# first Velmeshev (show in the main figures)

fdata <- list()

fdata$ctprofiles <- pdata$mkrs$minfc %>% dplyr::select(gene_id, cell_type, expr, dataset) %>% 
  filter(dataset == "velmeshev") %>% 
  mutate(expr = log2(expr + 1))

fdata$clusters <- read_rds(paste0(workspace$outputDir, paste0("sc_clusters_objects.rds")))
fdata$clusters <- fdata$clusters$clusterFlats %>% unname() %>% unlist(recursive = FALSE)

fdata$clusters <- fdata$clusters %>% session$collectionUtils$lapplyWithName(function(id, genes) {
  tibble(cluster = id, gene_id = genes)
}) %>% session$dataWrangler$rbind()

fdata$main <- fdata$clusters %>% inner_join(fdata$ctprofiles, by = "gene_id")

fdata$main <- fdata$main %>% group_by(cell_type, cluster) %>% summarize(expr = mean(expr)) %>% ungroup()

fdata$main <- fdata$main %>% group_by(cluster) %>% mutate(expr = ((expr - mean(expr)))) %>% ungroup()

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_14.eps"

xdata$main <-  fdata$main %>% mutate(cell_type = workspace$utils$fmtCelltypes(cell_type))

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cluster, y = cell_type)) + 
  geom_tile(aes(fill = expr)) + 
  scale_fill_gradient(low = "black", high = "yellow", na.value = "grey60") + 
  session$graphingUtils$tiltX(angle = 90) + 
  xlim(c(pdata$clusters$cluster[1:20], pdata$clusters$cluster[96:115])) +
  labs(fill = "Relative expression")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 35, height = 11) 

# =================================================

# ++++++++++++++++=============================
# display cell type profiles for the clusters

# first ROSMAP (show in the main figures)

fdata <- list()

fdata$ctprofiles <- pdata$mkrs$minfc %>% dplyr::select(gene_id, cell_type, expr, dataset) %>% 
  filter(dataset == "rosmap") %>% 
  mutate(expr = log2(expr + 1))

fdata$clusters <- read_rds(paste0(workspace$outputDir, paste0("sc_clusters_objects.rds")))
fdata$clusters <- fdata$clusters$clusterFlats %>% unname() %>% unlist(recursive = FALSE)

fdata$clusters <- fdata$clusters %>% session$collectionUtils$lapplyWithName(function(id, genes) {
  tibble(cluster = id, gene_id = genes)
}) %>% session$dataWrangler$rbind()

fdata$main <- fdata$clusters %>% inner_join(fdata$ctprofiles, by = "gene_id")

fdata$main <- fdata$main %>% group_by(cell_type, cluster) %>% summarize(expr = mean(expr)) %>% ungroup()

fdata$main <- fdata$main %>% group_by(cluster) %>% mutate(expr = ((expr - mean(expr)))) %>% ungroup()

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "sfigure_17_rosmap.png"

xdata$main <-  fdata$main %>% mutate(cell_type = workspace$utils$fmtCelltypes(cell_type))

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cluster, y = cell_type)) + 
  geom_tile(aes(fill = expr)) + 
  scale_fill_gradient(low = "black", high = "yellow", na.value = "grey60") + 
  session$graphingUtils$tiltX(angle = 90) + 
  xlim(c(pdata$clusters$cluster[1:20], pdata$clusters$cluster[96:115])) +
  labs(fill = "Relative expression", x = "Cluster", y = "Cell type")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "png", units = "cm", width = 35, height = 17, dpi = 400) 

# =================================================

# load in exprmats

pdata$files$exprmats <- tibble(outfile = list.files(workspace$outputDir)) %>% 
  filter(grepl("exprmats", outfile)) %>% 
  mutate(outfile_name = gsub(".rds", "", outfile)) %>% 
  filter(!grepl("RAW", outfile_name)) %>% 
  mutate(strs = strsplit(outfile_name, "_")) %>% 
  mutate(level = strs %>% sapply(function(currStr) { currStr[1] })) %>% 
  mutate(dataset = strs %>% sapply(function(currStr) { currStr[3] })) %>% 
  dplyr::select(outfile, level, dataset)

# +++++++++++++++++++++++++++++++++============================
# let's work out inter-cell type synchrony

fdata <- list(); gc()

fdata$exprmatFiles <- pdata$files$exprmats %>% filter(level == "sbj") #%>% filter(dataset %in% c("velmeshev", "rosmap"))
fdata$exprmatFiles <- fdata$exprmatFiles$outfile %>% session$dataWrangler$attachNames(fdata$exprmatFiles$dataset)

fdata$ctSyncSmry <- fdata$exprmatFiles %>% session$collectionUtils$lapplyWithName(function(datasetId, datasetFile) {
  currExprmat <- read_rds(paste0(workspace$outputDir, datasetFile))
  
  currGenes <- Reduce(intersect, lapply(currExprmat$exprmats, rownames))
  currSamples <- Reduce(intersect, lapply(currExprmat$exprmats, colnames))
  
  ctSyncTbl <- do.call(rbind, currGenes %>% session$collectionUtils$lapply(function(currGene) {
    coexmat <- currExprmat$exprmats %>% sapply(function(currExprmat) {
      currExprmat[currGene, currSamples]
    }) %>% t() %>% workspace$utils$computeCoex() 
    
    workspace$utils$getCoords(coexmat) %>% 
      mutate(cor_coef = workspace$utils$vectorize(coexmat)) %>% 
      mutate(gene_id = currGene) %>% 
      dplyr::select(gene_id, everything())
  }, verbose = FALSE))
  
  ctSyncSmry <- ctSyncTbl %>%
    group_by(gene_id) %>% 
    summarize(cor_coef = mean(cor_coef)) 
  
  ctSyncSmry %>% mutate(dataset = datasetId)
})

# normalize
fdata$ctSyncSmry <- fdata$ctSyncSmry %>% session$collectionUtils$lapply(function(ctSyncTbl) {
  ctSyncTbl %>% mutate(ct_sync = rank(cor_coef)) %>% mutate(ct_sync = ct_sync / max(ct_sync))
})

fdata$clusters <- read_rds(paste0(workspace$outputDir, paste0("sc_clusters_objects.rds")))
fdata$clusters <- fdata$clusters$clusterFlats %>% unname() %>% unlist(recursive = FALSE)


# ------------

fdata$mains <- fdata$clusters %>% session$collectionUtils$lapplyWithName(function(clusterId, currCluster) {
  
  fdata$ctSyncSmry %>% session$collectionUtils$lapplyWithName(function(currDataset, currCtSync) {
    currCtSync <- currCtSync %>% filter(gene_id %in% currCluster) 
    
    tibble(dataset = currDataset, 
           cluster = clusterId, 
           n_gene_total = length(currCluster), 
           n_gene_valid = nrow(currCtSync), 
           ct_sync = mean(currCtSync$ct_sync))
      
  }) %>% session$dataWrangler$rbind()
  
}) %>% session$dataWrangler$rbind()


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "figure_15.eps" # main figure

xdata$main <- fdata$mains %>% filter(dataset == "velmeshev")

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cluster, y = ct_sync)) + 
  geom_bar(stat = "identity") + 
  xlim(c(pdata$clusters$cluster[1:20], pdata$clusters$cluster[96:115])) + 
  session$graphingUtils$tiltX(angle = 90) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") 

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 35, height = 11) 
# =================================================


# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "sfigure_18_lau.eps" # main figure

xdata$main <- fdata$mains %>% filter(dataset == "lau")

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cluster, y = ct_sync)) + 
  geom_bar(stat = "identity") + 
  xlim(c(pdata$clusters$cluster[1:20], pdata$clusters$cluster[96:115])) + 
  session$graphingUtils$tiltX(angle = 90) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  ggtitle("Lau")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 35, height = 11) 
# =================================================

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "sfigure_19_lim.eps" # main figure

xdata$main <- fdata$mains %>% filter(dataset == "lim")

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cluster, y = ct_sync)) + 
  geom_bar(stat = "identity") + 
  xlim(c(pdata$clusters$cluster[1:20], pdata$clusters$cluster[96:115])) + 
  session$graphingUtils$tiltX(angle = 90) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  ggtitle("Lim")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 35, height = 11) 
# =================================================

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "sfigure_20_nagy.eps" # main figure

xdata$main <- fdata$mains %>% filter(dataset == "nagy")

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cluster, y = ct_sync)) + 
  geom_bar(stat = "identity") + 
  xlim(c(pdata$clusters$cluster[1:20], pdata$clusters$cluster[96:115])) + 
  session$graphingUtils$tiltX(angle = 90) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  ggtitle("Nagy")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 35, height = 11) 
# =================================================

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "sfigure_21_pineda.eps" # main figure

xdata$main <- fdata$mains %>% filter(dataset == "pineda")

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cluster, y = ct_sync)) + 
  geom_bar(stat = "identity") + 
  xlim(c(pdata$clusters$cluster[1:20], pdata$clusters$cluster[96:115])) + 
  session$graphingUtils$tiltX(angle = 90) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  ggtitle("Pineda")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 35, height = 11) 
# =================================================

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "sfigure_22_ramos.eps" # main figure

xdata$main <- fdata$mains %>% filter(dataset == "ramos")

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cluster, y = ct_sync)) + 
  geom_bar(stat = "identity") + 
  xlim(c(pdata$clusters$cluster[1:20], pdata$clusters$cluster[96:115])) + 
  session$graphingUtils$tiltX(angle = 90) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  ggtitle("Ramos")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 35, height = 11) 
# =================================================

# =================================================
# OUTPUT FIGURE
xdata <- list()

xdata$figName <- "sfigure_23_rosmap.eps" # main figure

xdata$main <- fdata$mains %>% filter(dataset == "rosmap")

xdata$plot <- xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cluster, y = ct_sync)) + 
  geom_bar(stat = "identity") + 
  xlim(c(pdata$clusters$cluster[1:20], pdata$clusters$cluster[96:115])) + 
  session$graphingUtils$tiltX(angle = 90) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  ggtitle("ROSMAP")

xdata$plot

ggsave(filename = paste0(pdata$figuresDir, xdata$figName), 
       plot = xdata$plot, device = "eps", units = "cm", width = 35, height = 11) 
# =================================================

 


fdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cluster, y = n_gene_valid)) + 
  geom_bar(stat = "identity") + 
  xlim(c(pdata$clusters$cluster[1:20], pdata$clusters$cluster[96:115])) + 
  session$graphingUtils$tiltX(angle = 90) +
  scale_y_continuous(trans = "log10")














# TODO - actually use the full expression matrix to compute synchrony!!! for both rosmap and Velmeshev; good enough to show these

# +++++++++++++++++++++++++++++++++============================
# let's work out correlations of cell type profiles

fdata <- list(); gc()

fdata$datasets <- pdata$ctprofiles$dataset %>% unique() %>% sort() %>% 
  session$dataWrangler$attachNames()

fdata$ctpCoexmats <- fdata$datasets %>% session$collectionUtils$lapply(function(currDataset) {
  mat <- pdata$ctprofiles %>% filter(dataset == currDataset) %>% 
    dplyr::select(gene_id, cell_type, expr) %>% 
    spread(cell_type, expr) %>% 
    session$dataWrangler$setColAsRownames("gene_id")
  
  mat[is.na(mat)] <- 0
  
  mat <- mat %>% workspace$utils$quantileNormalize()
  
  ctpCoexmat <- mat %>% workspace$utils$computeCoexmat()
  
  ctpCoexmat <- ctpCoexmat %>% workspace$utils$normalizeCoexmat()
  
  return(ctpCoexmat)
})

# +++++++++++++++++++++============
# COMMIT
fdata$ctpCoexmats %>% saveRDS(paste0(workspace$outputDir, "ctp_coexmats.rds"))
# +++++++++++++++++++++============


# +++++++++++++++++++++++++++++++++============================
# now, work out per cluster, average correlation of cell type profiles among the genes

fdata <- list(); gc()

fdata$clusters <- read_rds(paste0(workspace$outputDir, paste0("sc_clusters_objects.rds")))
fdata$clusters <- fdata$clusters$clusterFlats %>% unname() %>% unlist(recursive = FALSE) 

fdata$ctpCoexmats <- read_rds(paste0(workspace$outputDir, "ctp_coexmats.rds"))


fdata$main <- fdata$clusters %>% session$collectionUtils$lapplyWithName(function(clusterId, currCluster) {
  
  currCtpCoexmat <- fdata$ctpCoexmats$velmeshev
  currGenes <- currCluster %>% intersect(rownames(currCtpCoexmat))
  
  ctpCoex <- currCtpCoexmat[currGenes, currGenes] %>% workspace$utils$vectorize() %>% mean()
  
  tibble(cluster = clusterId, n_gene_total = length(currCluster), n_gene_valid = length(currGenes), ctp_coex = ctpCoex)
  
}) %>% session$dataWrangler$rbind()


fdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cluster, y = ctp_coex)) + 
  geom_bar(stat = "identity") + 
  xlim(c(pdata$clusters$cluster[1:20], pdata$clusters$cluster[96:115])) + 
  session$graphingUtils$tiltX(angle = 90) 
  


# !!!!!!!!!!!!!!!!!!!!
# TODO --- this is the code for visualizng cell type expression profile 

xdata <- list()

xdata$ctprofiles <-  pdata$ctprofiles %>% filter(dataset == "velmeshev") %>% 
  mutate(expr = log2(expr + 1)) %>%
  group_by(cell_type) %>% 
  mutate(expr = (expr - mean(expr)) / sd(expr)) %>% 
  ungroup()

xdata$ctprofiles %>% 
  filter(gene_id %in% fdata$clusters$inhibitory.id_0023) %>% 
  dplyr::select(-dataset) %>% 
  spread(cell_type, expr) %>% 
  session$dataWrangler$setColAsRownames("gene_id") %>% 
  session$dataWrangler$fillNa(colNames = pdata$celltypes, 0) %>% 
  # apply(1, function(vals) { (vals - mean(vals)) / sd(vals) }) %>% t() %>% 
  session$graphingUtils$heatmap(cluster_cols = FALSE, size = session$graphingUtils$MEDIUM)

xdata$ctprofiles %>% 
  filter(gene_id %in% pdata$rpgs$ensembl) %>% 
  dplyr::select(-dataset) %>% 
  spread(cell_type, expr) %>% 
  session$dataWrangler$setColAsRownames("gene_id") %>% 
  session$dataWrangler$fillNa(colNames = pdata$celltypes, 0) %>% 
  # apply(1, function(vals) { (vals - mean(vals)) / sd(vals) }) %>% t() %>% 
  session$dataWrangler$setRownameAsColumn("gene_id") %>% 
  gather(cell_type, expr, -gene_id) %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = expr)) + 
  geom_boxplot() + 
  ylim(-2.5, 2.5)

# END --- this is the code for visualizng cell type expression profile 

  
  
  
  # session$dataWrangler$setRownameAsColumn("gene_id") %>% 
  # gather(cell_type, expr_log, -gene_id) %>% 
  # session$graphingUtils$ggplot(aes(x = cell_type, y = gene_id)) + 
  # geom_tile(aes(fill = expr_log)) + scale_fill_gradient(low = "white", high = "red", na.value = "grey60")





















# ++++++++++++++++++++++++++==
# perform GO enrichment analysis on these 9 clusters to figure out the contents; one expectation is ribosomal protein genes; there may be others

fdata <- list()

fdata$clusters <- read_rds(paste0(workspace$outputDir, paste0("sc_clusters_objects.rds")))

fdata$clusterFlat <- fdata$clusters$clusterFlats %>% unname() %>% unlist(recursive = FALSE) 

fdata$clusterTbl <- fdata$clusterFlat %>% session$collectionUtils$lapplyWithName(function(name, genes) {
  tibble(cluster = name, gene_id = genes)
}) %>% session$dataWrangler$rbind()

fdata$genes <- fdata$clusterFlat %>% unlist() %>% unname() %>% unique()

# compute go enrichment for these 10 clusters

fdata$goAnalysis <- fdata$clusterFlat %>% session$collectionUtils$lapplyWithName(function(clusterName, clusterGenes) {
  
  celltype <- str_extract(clusterName, "^[a-z]+")
  bgGenes <- fdata$clusters$clusterObjs[[celltype]]$clusterTbl$gene_id %>% unique() %>% sort() # only use the genes that were "clustered" in the given cell type as the background
  negs <- bgGenes %>% setdiff(clusterGenes)
  
  pdata$go$flat %>% mclapply(function(goGenes) {
    
    trueGenes <- goGenes %>% intersect(bgGenes) 
    
    fisher <- session$evaluationUtils$fisher(predicted = clusterGenes, notPredicted = negs, trueSet = trueGenes, alternative = "greater")
    tibble(cluster = clusterName, odds_ratio = fisher$test$estimate, pvalue = fisher$test$p.value, 
           n_gene_go = length(goGenes), n_gene_bkgrnd = length(bgGenes), n_gene_true = length(trueGenes), n_gene_cluster = length(clusterGenes), 
           n_gene_tp = fisher$tbl[1, 1], n_gene_tn = fisher$tbl[2, 2], n_gene_fp = fisher$tbl[1, 2], n_gene_fn = fisher$tbl[2, 1])
    
  }, mc.cores = 20) %>% 
    session$collectionUtils$lapplyWithName(function(goName, tbl) { tbl %>% mutate(go_id = goName) }, verbose = FALSE) %>% 
    session$dataWrangler$rbind()
  
}) %>% session$dataWrangler$rbind() %>% dplyr::select(cluster, go_id, everything())

fdata$goAnalysis <- fdata$goAnalysis %>% left_join(pdata$go$labels, by = "go_id")
fdata$goAnalysis <- fdata$goAnalysis %>% mutate(qvalue = p.adjust(pvalue, method = "fdr")) # multiple test correction


fdata$goAnalysis %>% group_by(go_id) %>% filter(qvalue == min(qvalue)) %>% ungroup() %>% arrange(qvalue)
fdata$goAnalysis %>% group_by(cluster) %>% filter(qvalue == min(qvalue)) %>% ungroup() %>% arrange(qvalue)

# to summarize this, lets count the number of clusters a term is significant at FDR < 0.1

fdata$goTop <- fdata$goAnalysis %>% mutate(sig = qvalue < 0.1) %>% filter(sig) %>% group_by(go_id) %>% summarize(n_cluster = n()) %>% 
  left_join(pdata$go$labels, by = "go_id") %>% 
  arrange(desc(n_cluster))

fdata$goTop <- fdata$goTop %>% filter(n_cluster >= 5)

fdata$goAnalysis %>% group_by(cluster) %>% filter(odds_ratio == max(odds_ratio))

fdata$goAnalysis %>% group_by(cluster) %>% filter(qvalue == min(qvalue)) %>% ungroup() %>% arrange(qvalue)

# ========== COMMIT
fdata$goAnalysis %>% saveRDS(paste0(workspace$outputDir, paste0("sc_clusters_go_analysis.rds")))
# =============

# ++++++++++++++++++++++++++==
# Let's now look at xSubject level preservation of xCell clusters; xSubject compared to xCell

fdata <- list(); gc()

fdata$go <- readRDS(paste0(workspace$outputDir, paste0("sc_clusters_go_analysis.rds")))

fdata$main <- readRDS(paste0(workspace$outputDir, paste0("sc_clusters_integrity_within_celltype.rds")))

fdata$main <- fdata$main$egad %>% filter(dataset != "consensus", !grepl("generic", cluster), level != "bk") 

fdata$means <- fdata$main %>% group_by(cluster, level) %>% summarize(auc = mean(auc)) %>% ungroup()

fdata$meansMat <- fdata$means %>% spread(level, auc) 
 
fdata$main <- fdata$main %>% dplyr::select(cluster, level, dataset, auc) %>% 
  mutate(network = paste0(level, ".", dataset))

fdata$design <- fdata$main %>% 
  dplyr::select(-auc, -cluster) %>% unique() %>% 
  mutate(level = factor(level, levels = c("sc", "sbj"))) %>% 
  arrange(level, dataset)

fdata$mat <- fdata$main %>% 
  dplyr::select(cluster, network, auc) %>% 
  spread(network, auc) %>% 
  session$dataWrangler$setColAsRownames("cluster")

fdata$clusters <- rownames(fdata$mat) %>% session$dataWrangler$attachNames()

fdata$aucDiff <- fdata$clusters %>% session$collectionUtils$lapply(function(cluster) {
  dat <- fdata$mat[cluster, ] %>% t() %>% session$dataWrangler$setRownameAsColumn("network")
  names(dat) <- c("network", "auc")
  dat <- fdata$design %>% left_join(dat, by = "network") %>% dplyr::select(level, auc) %>% na.omit()
  res <- wilcox.test(auc ~ level, data = dat, paired = TRUE, alternative = "greater")
  stat <- dat %>% group_by(level) %>% summarize(auc = mean(auc)) %>% spread(level, auc) %>% mutate(fc = sbj / sc, diff = sbj - sc)
  stat %>% mutate(cluster = cluster, pvalue = res$p.value) %>% dplyr::select(cluster, everything())
}) %>% session$dataWrangler$rbind()

fdata$aucDiff <- fdata$aucDiff %>% mutate(qvalue = p.adjust(pvalue, method = "fdr"))

fdata$aucDiff %>% filter(qvalue > 0.1) 

fdata$aucDiff %>% session$graphingUtils$ggplot(aes(x = pvalue)) + geom_histogram()


fdata$aucDiff %>% arrange(desc(sbj)) -> x
fdata$aucDiff %>% arrange((diff)) -> x


fdata$main %>% filter(cluster == "oligodendrocyte.id_0007") %>% 
  session$graphingUtils$ggplot(aes(x = level, y = auc)) + 
  geom_line(aes(group = dataset)) + 
  geom_text(data = fdata$main %>% filter(cluster == "oligodendrocyte.id_0007", level == "sc"), aes(label = dataset), vjust = -1, size = 4) +
  geom_point() + 
  ylim(0.5, 1)

fdata$mat[x$cluster, ] %>% session$graphingUtils$heatmap(size = session$graphingUtils$SMALL, cluster_rows = FALSE)









# ++++++++++++++++++++++++++==
# Let's now look at xSubject level preservation of xCell clusters; xBulk compared to xSubject

fdata <- list(); gc()

fdata$go <- readRDS(paste0(workspace$outputDir, paste0("sc_clusters_go_analysis.rds")))

fdata$main <- readRDS(paste0(workspace$outputDir, paste0("sc_clusters_integrity_within_celltype.rds")))

fdata$main <- fdata$main$egad %>% filter(dataset != "consensus", !grepl("generic", cluster), !grepl("-", dataset)) 

fdata$means <- fdata$main %>% group_by(cluster, level) %>% summarize(auc = mean(auc)) %>% ungroup()

fdata$means %>% session$graphingUtils$ggplot(aes(x = auc)) + geom_density(aes(fill = level), alpha = 0.5) # density showing how cluster integrity changes

fdata$meansMat <- fdata$means %>% spread(level, auc) 

fdata$main <- fdata$main %>% dplyr::select(cluster, level, dataset, auc) %>% 
  mutate(network = paste0(level, ".", dataset))

fdata$design <- fdata$main %>% 
  dplyr::select(-auc, -cluster) %>% unique() %>% 
  mutate(level = factor(level, levels = c("sc", "sbj", "bk"))) %>% 
  arrange(level, dataset)

fdata$mat <- fdata$main %>% 
  dplyr::select(cluster, network, auc) %>% 
  spread(network, auc) %>% 
  session$dataWrangler$setColAsRownames("cluster")

fdata$clusters <- rownames(fdata$mat) %>% session$dataWrangler$attachNames()

fdata$smry <- fdata$clusters %>% session$collectionUtils$lapply(function(cluster) {
  dat <- fdata$mat[cluster, ] %>% t() %>% session$dataWrangler$setRownameAsColumn("network")
  names(dat) <- c("network", "auc")
  dat <- fdata$design %>% left_join(dat, by = "network") %>% dplyr::select(level, auc) %>% na.omit()
  stat <- dat %>% group_by(level) %>% summarize(auc = mean(auc)) %>% spread(level, auc) 
  stat %>% mutate(cluster = cluster) %>% dplyr::select(cluster, everything())
}) %>% session$dataWrangler$rbind()

# fdata$aucDiff <- fdata$clusters %>% session$collectionUtils$lapply(function(cluster) {
#   dat <- fdata$mat[cluster, ] %>% t() %>% session$dataWrangler$setRownameAsColumn("network")
#   names(dat) <- c("network", "auc")
#   dat <- fdata$design %>% left_join(dat, by = "network") %>% dplyr::select(level, auc) %>% na.omit()
#   res <- wilcox.test(auc ~ level, data = dat, paired = FALSE, alternative = "two.sided")
#   stat <- dat %>% group_by(level) %>% summarize(auc = mean(auc)) %>% spread(level, auc) %>% mutate(fc = bk / sbj, diff = bk - sbj)
#   stat %>% mutate(cluster = cluster, pvalue = res$p.value) %>% dplyr::select(cluster, everything())
# }) %>% session$dataWrangler$rbind()

fdata$aucDiff <- fdata$aucDiff %>% mutate(qvalue = p.adjust(pvalue, method = "fdr"))

fdata$aucDiff %>% filter(qvalue < 0.1) # nothing is significant at fdr < 0.1; but there is a clear pile up at the low end

fdata$aucDiff %>% session$graphingUtils$ggplot(aes(x = diff)) + geom_histogram()
fdata$aucDiff %>% session$graphingUtils$ggplot(aes(x = pvalue)) + geom_histogram()


fdata$aucDiff %>% arrange(desc(diff)) 
fdata$aucDiff %>% arrange((diff))

fdata$aucDiff %>% arrange(desc(bk)) -> x

fdata$aucDiff %>% filter(bk > 0.8) # there are 24 clusters out of 115 where EGAD-AUROC > 0.8

fdata$mat[x$cluster, rev(colnames(fdata$mat))] %>% 
  session$graphingUtils$heatmap(size = session$graphingUtils$LARGE, cluster_rows = FALSE, cluster_cols = FALSE)

xdata <- list()

xdata$cluster <- "excitatory.id_0020"

fdata$main %>% filter(cluster == xdata$cluster ) %>% 
  session$graphingUtils$ggplot(aes(x = level, y = auc)) + 
  geom_line(aes(group = dataset)) + 
  geom_text(data = fdata$main %>% filter(cluster == xdata$cluster , level == "sbj"), aes(label = dataset), vjust = -1, size = 4) +
  geom_point() + 
  ylim(0.5, 1)


# filter for high performing clusters

fdata$main %>% filter(level == "bk", dataset %in% c("velmeshev", "rosmap"), auc > 0.8) -> x
x %>% group_by(cluster) %>% summarize(auc = mean(auc)) %>% arrange(desc(auc)) -> x

fdata$go %>% filter(cluster %in% x$cluster) %>% group_by(cluster) %>% filter(pvalue == min(pvalue)) %>% arrange(desc(odds_ratio)) %>% ungroup() %>% View()

# let's see what these clusters represent

fdata$aucDiff %>% filter(bk >= 0.8) %>% arrange(desc(fc)) 

fdata$goSmry <- fdata$go %>% group_by(cluster) %>% filter(pvalue == min(pvalue)) %>% ungroup()

fdata$aucDiff %>% dplyr::select(cluster, bk) %>% left_join(fdata$goSmry %>% dplyr::select(cluster, odds_ratio, pvalue, go_label), by = "cluster") %>% 
  arrange(desc(bk)) %>% View()

fdata$mat %>% session$dataWrangler$setRownameAsColumn("cluster") %>% dplyr::select(cluster, bk.velmeshev, bk.rosmap) %>% 
  left_join(fdata$goSmry %>% dplyr::select(cluster, odds_ratio, pvalue, go_label), by = "cluster") %>% 
  arrange(desc(bk.velmeshev)) -> x




fdata$avgcoef <- readRDS(paste0(workspace$outputDir, paste0("sc_clusters_integrity_in_bk_coefmats.rds")))

fdata$dataset <- c("velmeshev", "rosmap") %>% session$dataWrangler$attachNames()

fdata$mgpDiff <- fdata$dataset %>% session$collectionUtils$lapply(function(currDataset) {
  main <- fdata$avgcoef %>% filter(dataset %in% c(paste0(currDataset), paste0(currDataset, "-mgpres"))) %>% 
    dplyr::select(cluster, dataset, coex) %>% 
    spread(dataset, coex)
  colnames(main) <- c("cluster", "original", "residual")
  main <- main %>% mutate(diff = residual - original) %>% mutate(dataset = currDataset, procedure = "mgp") 
  return(main)
}) %>% session$dataWrangler$rbind()

# fdata$mgpDiffSmry <- fdata$mgpDiff %>% dplyr::select(cluster, dataset, diff) %>% spread(dataset, diff)

fdata$psdbkDiff <- fdata$dataset %>% session$collectionUtils$lapply(function(currDataset) {
  main <- fdata$avgcoef %>% filter(dataset %in% c(paste0(currDataset), paste0(currDataset, "-psdbkres"))) %>% 
    dplyr::select(cluster, dataset, coex) %>% 
    spread(dataset, coex)
  colnames(main) <- c("cluster", "original", "residual")
  main <- main %>% mutate(diff = residual - original) %>% mutate(dataset = currDataset, procedure = "psdbk")
  return(main)
}) %>% session$dataWrangler$rbind()

# fdata$psdbkDiffSmry <- fdata$mgpDiff %>% dplyr::select(cluster, dataset, diff) %>% spread(dataset, diff)

fdata$diffs <- fdata$mgpDiff %>% rbind(fdata$psdbkDiff)


xdata <- list()

xdata$ylim <- (x %>% arrange(desc(auc)))$cluster %>% unique()

fdata$diffs %>% filter(cluster %in% xdata$ylim) %>% 
  filter(procedure == "mgp") %>% 
  session$graphingUtils$ggplot(aes(x = diff, y = cluster)) +
  geom_bar(stat = "identity") + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  facet_wrap(~dataset) + 
  ylim(rev(xdata$ylim)) + 
  xlim(-0.4, 0.2) + 
  ggtitle("Average cluster co-expression change after CCV correction")

fdata$go %>% filter(cluster == "astrocyte.id_0007") %>% arrange(pvalue) %>% dplyr::select(cluster, go_id, go_label, odds_ratio, pvalue)


x %>% left_join(fdata$mgpDiffSmry, by = "cluster") %>% View()


fdata$ccv %>% session$graphingUtils$ggplot(aes(x = original, y = residual)) + geom_text(aes(label = cluster)) 

# ================= END OF ONGOING ANALYSIS = THE REMAINDER REQUIRE EXTENSIVE REFACTORING






































































# ++++++++++++++++++++++++++==
# characterize how clusters are cell type specific

fdata <- list(); gc()

fdata$clusters <- read_rds(paste0(workspace$outputDir, paste0("sc_clusters_objects.rds")))

fdata$main <- readRDS(paste0(workspace$outputDir, paste0("sc_clusters_integrity_across_celltype.rds")))

# let's first characterize cell type specificity for each cell type

# EGAD

fdata$egad <- fdata$main$egad %>% filter(!grepl("^generic", cluster))

fdata$egad <- fdata$egad %>% 
  dplyr::select(cluster, coex_ct, auc) %>% spread(coex_ct, auc) %>% na.omit() %>% 
  session$dataWrangler$setColAsRownames("cluster") 

fdata$fcEgad <- pdata$celltypes %>% sapply(function(currCelltype) {
  currCt <- fdata$egad[, currCelltype]
  othrCtMax <- fdata$egad[, colnames(fdata$egad) != currCelltype] %>% apply(1, max) 
  minfc <- currCt / othrCtMax
  return(minfc)
}) %>% session$dataWrangler$setRownameAsColumn("cluster")

xdata <- list()

xdata$main <- fdata$fcEgad %>% reshape2::melt(id = "cluster") %>% as_tibble() %>% 
  dplyr::select(cluster, cell_type = variable, min_fc = value) %>% 
  mutate(ct_cluster = cluster %>% sapply(function(currStr) { str_extract(currStr, "^[a-z]+") }))

xdata$main %>% 
  session$graphingUtils$ggplot(aes(y = cluster, x = min_fc)) + 
  geom_bar(stat = "identity", fill = "grey70") + 
  geom_bar(data = xdata$main %>% filter(cell_type == ct_cluster), stat = "identity") +
  facet_wrap(~cell_type, nrow = 1) + 
  geom_vline(xintercept = 1, linetype = "dashed") + 
  theme(axis.text = element_text(size = 12)) +
  ggtitle("EGAD")


# avg coex

fdata$avgcoex <- fdata$main$avgcoex %>% filter(!grepl("^generic", cluster))

fdata$avgcoex <- fdata$avgcoex %>% 
  dplyr::select(cluster, coex_ct, coex) %>% spread(coex_ct, coex) %>% na.omit() %>% 
  session$dataWrangler$setColAsRownames("cluster") 

fdata$fcAvgcoex <- pdata$celltypes %>% sapply(function(currCelltype) {
  currCt <- fdata$avgcoex[, currCelltype]
  othrCtMax <- fdata$avgcoex[, colnames(fdata$avgcoex) != currCelltype] %>% apply(1, max) 
  minfc <- currCt / othrCtMax
  return(minfc)
}) %>% session$dataWrangler$setRownameAsColumn("cluster")

xdata <- list()

xdata$main <- fdata$fcAvgcoex %>% reshape2::melt(id = "cluster") %>% as_tibble() %>% 
  dplyr::select(cluster, cell_type = variable, min_fc = value) %>% 
  mutate(ct_cluster = cluster %>% sapply(function(currStr) { str_extract(currStr, "^[a-z]+") }))

xdata$main %>% 
  session$graphingUtils$ggplot(aes(y = cluster, x = min_fc)) + 
  geom_bar(stat = "identity", fill = "grey70") + 
  geom_bar(data = xdata$main %>% filter(cell_type == ct_cluster), stat = "identity") +
  facet_wrap(~cell_type, nrow = 1) + 
  geom_vline(xintercept = 1, linetype = "dashed") + 
  theme(axis.text = element_text(size = 12)) + 
  ggtitle("Intra-cluster co-expression")


# plot out particular cluster performance ****

xdata <- list()

xdata$cluster <- "inhibitory.id_0018"

xdata$avgcoex <- fdata$main$avgcoex %>% filter(cluster == xdata$cluster)
xdata$egad <- fdata$main$egad %>% filter(cluster == xdata$cluster)

xdata$main <- xdata$avgcoex %>% dplyr::select(coex_ct, coex) %>% 
  left_join(xdata$egad %>% dplyr::select(coex_ct, auc), by = c("coex_ct"))

xdata$main <- xdata$main %>% reshape2::melt(id = "coex_ct") %>% as_tibble() %>% 
  dplyr::select(cell_type = coex_ct, metric = variable, value) %>% 
  left_join(xdata$avgcoex %>% dplyr::select(cell_type = coex_ct, n_lnks), by = "cell_type")

xdata$main %>% session$graphingUtils$ggplot(aes(x = cell_type, y = value, color = metric)) + 
  geom_point(aes(size = n_lnks), alpha = 0.8) + 
  session$graphingUtils$tiltX(angle = 90) + 
  ylim(0, 1) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  ggtitle(paste0("Cluster integrity of ", xdata$cluster))

# check expression biases
pdata$ctprofiles %>% 
  filter(gene_id %in% fdata$clusters$clusterFlats$astrocyte[[xdata$cluster]]) %>% 
  session$graphingUtils$ggplot(aes(y = cell_type, x = expr)) + geom_boxplot() + facet_wrap(~dataset, ncol = 1) + 
  scale_x_continuous(trans = "log10") +
  ggtitle(paste0("expression levels of ", xdata$cluster))


# ++++++++++++++++++++++++++==
# next, look at clusters that are generic

fdata <- list(); gc()

fdata$clusters <- read_rds(paste0(workspace$outputDir, paste0("sc_clusters_objects.rds")))

fdata$main <- readRDS(paste0(workspace$outputDir, paste0("sc_clusters_integrity_across_celltype.rds")))

# for the generic clusters; pick the one with the highest number of genes

fdata$egad <- fdata$main$egad %>% group_by(cluster, coex_ct) %>% filter(n_gene == max(n_gene)) %>% summarize(auc = mean(auc)) %>% ungroup()

fdata$egad <- fdata$egad %>% 
  dplyr::select(cluster, coex_ct, auc) %>% spread(coex_ct, auc) %>% na.omit() %>% 
  session$dataWrangler$setColAsRownames("cluster") 

fdata$egadGeneric <- fdata$egad %>% 
  apply(1, function(x) { min(x, na.rm = TRUE) }) %>% session$dataWrangler$vectorToTibble() %>% 
  dplyr::select(cluster = variable, auc_generic = value) %>% 
  arrange(desc(auc_generic)) 

fdata$genericClusters <- fdata$egadGeneric %>% filter(auc_generic > 0.9) %>% arrange(cluster) # these are the 9 clusters

fdata$egad[fdata$egadGeneric$cluster, ] %>% 
  session$graphingUtils$heatmap(size = session$graphingUtils$MEDIUM, cluster_row = FALSE, cluster_col = FALSE)


fdata$avgcoex <- fdata$main$avgCoex

fdata$avgcoex <- fdata$main$avgcoex %>% group_by(cluster, coex_ct) %>% filter(n_gene == max(n_gene)) %>% summarize(coex = mean(coex)) %>% ungroup()

fdata$avgcoex <- fdata$avgcoex %>% 
  dplyr::select(cluster, coex_ct, coex) %>% spread(coex_ct, coex) %>% na.omit() %>% 
  session$dataWrangler$setColAsRownames("cluster") 

fdata$avgcoex[fdata$egadGeneric$cluster, ] %>% 
  session$graphingUtils$heatmap(size = session$graphingUtils$MEDIUM, cluster_row = FALSE, cluster_col = FALSE)

fdata$avgcoexGeneric <- fdata$avgcoex %>% 
  apply(1, function(x) { min(x, na.rm = TRUE) }) %>% session$dataWrangler$vectorToTibble() %>% 
  dplyr::select(cluster = variable, coex_generic = value) %>% 
  arrange(desc(coex_generic)) 

xdata <- list()

xdata$celltype <- "astrocyte"

plotDendroAndColors(fdata$clusters$clusterObjs[[xdata$celltype]]$dendro, 
                    colors = labels2colors(fdata$clusters$clusterObjs[[xdata$celltype]]$clusters), 
                    rowText = fdata$clusters$clusterObjs[[xdata$celltype]]$clusters, 
                    dendroLabels = FALSE)


# ++++++++++++++++++++++++++==
# now, let's put together the generic clusters
# cluster flat + cluster mat so that they could be saved and used for later analyses

# cluster                 auc_generic
# <chr>                         <dbl>
# 1 astrocyte.id_0002             0.947
# 2 astrocyte.id_0004             0.943
# 3 excitatory.id_0001            0.959
# 4 inhibitory.id_0001            0.948
# 5 microglia.id_0001             0.964
# 6 oligodendrocyte.id_0001       0.942
# 7 oligodendrocyte.id_0007       0.937
# 8 opc.id_0001                   0.942
# 9 opc.id_0002                   0.926



# grab all the genes that belong to any cluster for each GO term

fdata$goTop <- fdata$goTop$go_id
names(fdata$goTop) <- fdata$goTop

fdata$goSets <- fdata$goTop %>% session$collectionUtils$lapply(function(goId) {
  genes <- pdata$go$flat[[goId]] %>% intersect(fdata$genes)
  tibble(go_id = goId, gene_id = genes)
}) %>% session$dataWrangler$rbind() %>% 
  left_join(pdata$go$labels, by = "go_id")

fdata$goSets %>% dplyr::select(gene_id, go_id) %>% spread(go_id)

fdata$goSetsMat <- fdata$goSets %>% 
  dplyr::select(gene_id, go_id) %>% 
  mutate(in_set = 1) %>% 
  spread(go_id, in_set)
fdata$goSetsMat <- fdata$goSetsMat %>% session$dataWrangler$setColAsRownames("gene_id")
fdata$goSetsMat[is.na(fdata$goSetsMat)] <- 0

# show upset plot
fdata$goSetsMat %>% UpSetR::upset(nsets = 10, text.scale = 2, keep.order = TRUE, nintersects = NA)


# OK, I have narrowed it down to 3 go terms
# GO:0006412; translation
# GO:1902600; proton transmembrane transport

fdata$goSetsMat %>% UpSetR::upset(sets = c("GO:0006412", "GO:1902600"), 
                                  text.scale = 2, keep.order = TRUE, nintersects = NA)

pdata$go$tbl %>% filter(go_id %in% c("GO:0006412", "GO:1902600")) %>% group_by(go_id) %>% summarize(n = n())

# ok pull together the genes in these two clusters now. 


fdata$genericGoIds <- c(generic.translation = "GO:0006412", generic.proton_transport = "GO:1902600")

fdata$genericClustersTbl <- fdata$goSets %>% filter(go_id %in% fdata$genericGoIds) 

fdata$genericClustersFlat <- fdata$genericClustersFlat %>% lapply(function(currGoId) { 
  fdata$goSets %>% filter(go_id == currGoId) %>% session$dataWrangler$extractColumn("gene_id")
    
})

fdata$genericClustersMat <- fdata$genericClustersTbl %>% dplyr::select(gene_id, go_id) %>% mutate(in_set = 1) %>% spread(go_id, in_set)
fdata$genericClustersMat <- fdata$genericClustersMat %>% session$dataWrangler$setColAsRownames("gene_id")
fdata$genericClustersMat[is.na(fdata$genericClustersMat)] <- 0
fdata$genericClustersMat <- fdata$genericClustersMat[, fdata$genericGoIds]
colnames(fdata$genericClustersMat) <- names(fdata$genericGoIds)

# add in negative set over here - which is all the genes? generic genes?

# ========== COMMIT
list(clusterTbl = fdata$genericClustersTbl, clusterFlat = fdata$genericClustersFlat, clusterMat = fdata$genericClustersMat) %>% 
  saveRDS(paste0(workspace$outputDir, paste0("sc_generic_clusters.rds")))
# =============


































# TODO SANDBOX - TO BE DELETED

fdata <- list(); gc()

fdata$main <- readRDS(paste0(workspace$outputDir, paste0("sc_clusters_integrity_in_bk_coefmats.rds")))

fdata$dataset <- "velmeshev"
fdata$celltype <- "inhibitory"

fdata$ccv <- fdata$main %>% filter(dataset %in% c(paste0(fdata$dataset), paste0(fdata$dataset, "-mgpres"))) %>% 
  dplyr::select(cluster, dataset, coex) %>% 
  spread(dataset, coex)

colnames(fdata$ccv) <- c("cluster", "original", "residual")

fdata$ccv %>% session$graphingUtils$ggplot(aes(x = original, y = residual)) + geom_text(aes(label = cluster)) + geom_abline(slope = 1, linetype = "dashed")

fdata$ccv <- fdata$ccv %>% mutate(diff = residual - original) %>% mutate(procedure = "mgp")

fdata$psdbk <- fdata$main %>% filter(dataset %in% c(paste0(fdata$dataset, "-psdbk"), paste0(fdata$dataset, "-psdbkres"))) %>% 
  dplyr::select(cluster, dataset, coex) %>% 
  spread(dataset, coex)

colnames(fdata$psdbk) <- c("cluster", "original", "residual")

fdata$psdbk %>% session$graphingUtils$ggplot(aes(x = original, y = residual)) + geom_text(aes(label = cluster)) + geom_abline(slope = 1, linetype = "dashed")

fdata$psdbk <- fdata$psdbk %>% mutate(diff = residual - original) %>% mutate(procedure = "psdbk")

fdata$main2 <- fdata$ccv %>% rbind(fdata$psdbk) %>% 
  dplyr::select(cluster, procedure, diff) 

ordered <- fdata$psdbk %>% arrange(diff)

fdata$ccv %>% rbind(fdata$psdbk) %>% 
  dplyr::select(cluster, procedure, diff) %>% 
  session$graphingUtils$ggplot(aes(x = diff, y = cluster, color = procedure)) + geom_point(size = 4) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  ylim(unique(rev(xdata$order$cluster[1:50]))) + 
  ggtitle(paste0(fdata$dataset))

fdata$psdbk %>% 
  session$graphingUtils$ggplot(aes(x = diff, y = cluster)) + geom_point(size = 4) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  ylim(unique(rev(ordered$cluster[1:50]))) + 
  ggtitle(paste0(fdata$dataset))


# inhibitory.id_0012 # represents one that 
# 
fdata$clusters <- read_rds(paste0(workspace$outputDir, paste0("sc_clusters_objects.rds")))
fdata$clustersGeneric <- read_rds(paste0(workspace$outputDir, paste0("sc_generic_clusters.rds")))
fdata$clusterFlats <- c(fdata$clusters$clusterFlats, list(generic = fdata$clustersGeneric$clusterFlat))

pdata$genes %>% filter(gene_id %in% fdata$clusters$clusterFlats$inhibitory$inhibitory.id_0012) %>% View()

fdata$xSbjSd <- read_rds(paste0(workspace$outputDir, paste0("bk_exprmats_velmeshev-psdbkres.rds")))$xSbjSd
fdata$xCtSd <- read_rds(paste0(workspace$outputDir, paste0("umkrs.rds")))$minfc

fdata$xSbjSd$minfc %>% filter(gene_id %in% fdata$clusters$clusterFlats$inhibitory$inhibitory.id_0003)

fdata$psdbkres <- read_rds(paste0(workspace$outputDir, paste0("bk_exprmats_velmeshev-psdbkres.rds")))

fdata$mgpres <- read_rds(paste0(workspace$outputDir, paste0("bk_exprmats_velmeshev-mgpres.rds")))

fdata$mgpres$mgps$rotations %>% session$collectionUtils$lapplyWithName(function(celltype, rotations) { 
  rotations[, "PC1"] %>% session$dataWrangler$setRownameAsColumn("gene_id") %>% mutate(cell_type = celltype)
}) %>% session$dataWrangler$rbind() -> x

fdata$psdbkres$xSbjSd$minfc %>% right_join(x, by = c("cell_type", "gene_id")) -> y

y %>% mutate(sd = expr) %>% 
  filter(cell_type == "excitatory") %>% 
  session$graphingUtils$ggplot(aes(x = expr, y = .)) + geom_point() + facet_wrap(~cell_type)


fdata$xCtSd %>%
  filter(cell_type == "excitatory") %>% 
  mutate(in_cluster = gene_id %in% fdata$clusterFlats$excitatory$excitatory.id_0006) %>% 
  session$graphingUtils$ggplot(aes(x = in_cluster, y = min_fc_rank)) + geom_boxplot()

# TODO ---- THIS PART IS VERY IMPORTANT FOR VISUALIZING THE CONTRIBUTION OF DIFFERENT CELL TYPES

xdata$main <- fdata$psdbkres$lmModel$coefStats %>% dplyr::select(gene_id, cell_type, rsqr_inde) %>% 
  spread(cell_type, rsqr_inde) %>% 
  session$dataWrangler$setColAsRownames("gene_id") %>% apply(1, function(vals) { (vals - mean(vals)) / sd(vals) }) %>% t()
xdata$genes <- fdata$clusterFlats$inhibitory$inhibitory.id_0012 %>% intersect(rownames(xdata$main))
xdata$full <- fdata$psdbkres$lmModel$lmStats %>% filter(gene_id %in% xdata$genes) %>% dplyr::select(gene_id, rsqr) %>% 
  session$dataWrangler$setColAsRownames("gene_id")
xdata$main <- xdata$main[xdata$genes, ] %>% cbind(xdata$full[xdata$genes, "rsqr", drop = FALSE])
xdata$main[xdata$genes, pdata$celltypes] %>% session$graphingUtils$heatmap(size = session$graphingUtils$MEDIUM, show_rownames = FALSE)
# xdata$main[xdata$genes, ] %>% session$graphingUtils$heatmap(size = session$graphingUtils$MEDIUM)


fdata$psdbkres$lmModel$coefStats %>% dplyr::select(gene_id, cell_type, rsqr_boost) %>% 
  filter(gene_id %in% fdata$clusterFlats$generic$generic.translation) %>% 
  session$graphingUtils$ggplot(aes(x = rsqr_boost)) + geom_density(aes(color = cell_type))


# TODO ----------- THIS CODE IS REUSABLE NOW LET'S TAKE A LOOK AT CELL TYPE SYNCHRONY

fdata$psdbkres$exprmat$sbj %>% lapply(dim)


# compute cell type synchrony -- separately and for all genes 
# VELMESHEV

fdata <- list()

fdata$psdbkres <- read_rds(paste0(workspace$outputDir, paste0("bk_exprmats_rosmap-psdbkres.rds")))
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
  summarize(cor_coef = mean(cor_coef)) 

fdata$ctSyncSmry <- fdata$ctSyncSmry %>% 
  left_join(pdata$genes %>% dplyr::select(gene_id = ensembl, gene), by = "gene_id")

fdata$ctSyncSmry <- fdata$ctSyncSmry %>% 
  mutate(cor_perc = rank(cor_coef) / length(cor_coef))

fdata$ctSyncSmry <- fdata$ctSyncSmry %>% dplyr::select(gene_id, gene, cor_coef, cor_perc)

fdata$ctSyncTbl %>% filter(gene_id %in% fdata$clusterFlats$microglia$microglia.id_0006) %>% 
  group_by(gene_a, gene_b) %>% summarize(mean(cor_coef))

fdata$ctSyncSmry %>% filter(gene_id %in% fdata$clusterFlats$excitatory$excitatory.id_0023) %>% 
  session$graphingUtils$ggplot(aes(x = cor_perc)) + geom_density()

fdata$x <- fdata$clusterFlats %>% unname() %>% unlist(recursive = FALSE) %>% 
  session$collectionUtils$lapplyWithName(function(id, currCluster) {
    fdata$ctSyncSmry %>% filter(gene_id %in% currCluster) %>% 
      session$dataWrangler$extractColumn("cor_perc") %>% mean() -> x
    tibble(cluster = id, sync = x)
  }) %>% session$dataWrangler$rbind()

fdata$x %>% left_join(fdata$psdbk, by = "cluster") %>% session$graphingUtils$ggplot(aes(x = sync, y = diff)) + 
  geom_point() + 
  geom_label(data = fdata$x %>% left_join(fdata$psdbk, by = "cluster") %>% 
               filter(cluster %in% c("generic.translation", "excitatory.id_0034", 
                                     "inhibitory.id_0012", "excitatory.id_0010", 
                                     "microglia.id_0001", "astrocyte.id_0006",
                                     "oligodendrocyte.id_0014", "excitatory.id_0023", "astrocyte.id_0011")), aes(label = cluster))

fdata$x %>% session$graphingUtils$ggplot(aes(x = sync)) + geom_histogram() + 
  geom_histogram(data = fdata$x %>% filter(cluster %in% c("excitatory.id_0010", "excitatory.id_0020", "inhibitory.id_0012")) , fill = "red")


fdata$x %>% session$graphingUtils$ggplot(aes(x = sync)) + geom_histogram() +
  geom_histogram(data = fdata$x %>% filter(cluster %in% c("microglia.id_0001", "opc.id_0001", "astrocyte.id_0002")), fill = "red") 

fdata$x %>% session$graphingUtils$ggplot(aes(x = sync)) + geom_histogram() + 
  geom_histogram(data = fdata$x %>% filter(cluster %in% c("excitatory.id_0023", "microglia.id_0006", "astrocyte.id_0007")), fill = "red")
  
  











# look at the correlation between synchrony and cluster integrity

fdata$clusterSync <- fdata$clusterFlats %>% unname() %>% unlist(recursive = FALSE) %>% session$collectionUtils$lapplyWithName(function(id, genes) {
  fdata$ctSyncSmry %>% filter(gene_id %in% genes) %>% summarize(cor_coef = mean(cor_coef), cor_perc = mean(cor_perc)) %>% 
    mutate(cluster = id)
}) %>% session$dataWrangler$rbind()

fdata$clusterSync <- fdata$clusterSync %>% dplyr::select(cluster, everything())

fdata$clusterSync %>% left_join(ungroup(xdata$egad), by = "cluster") %>% 
  filter(level == "bk") %>% 
  session$graphingUtils$ggplot(aes(x = cor_perc, y = auc)) + geom_point() +
  facet_wrap(~level)

fdata$clusterSync %>% left_join(ungroup(xdata$egad), by = "cluster") %>% filter(level == "bk") -> x
x$cor_perc %>% cor.test(x$auc)


fdata$clusterSync %>% left_join(ungroup(xdata$avgcoex), by = "cluster") %>% 
  filter(level == "bk") %>% 
  session$graphingUtils$ggplot(aes(x = cor_perc, y = coex)) + geom_point() +
  facet_wrap(~level)

fdata$clusterSync %>% left_join(ungroup(xdata$avgcoex), by = "cluster") %>% 
  filter(level == "bk") -> x
x$cor_perc %>% cor.test(x$coex)

fdata$clusterSync %>% left_join(z, by = "cluster") %>% 
  session$graphingUtils$ggplot(aes(x = cor_perc, y = abs(diff))) + 
  geom_point()


# HERE IS THE STORY: 
# a lot of subject level stuff is getting passed through to the bulk tissue level (so xSc > xSubject is the real constraint)
# for the xCell clusters anyway
# ----> I know this because the cluster integrity is highly correlated so those that are highly co-expressed are also highly co-expressed at the xBulk level
# ----> and getting rid of the xSubject signal leads to reduction in xBulk co-expression
# ----> So how are those high co-expressions being propagated through? By synchrony! 


xdata <- list(); gc()

xdata$main <- readRDS(paste0(workspace$outputDir, paste0("sc_clusters_integrity_within_celltype.rds")))

xdata$egad <- xdata$main$egad %>% 
  filter(!grepl("-", cluster)) %>% 
  group_by(cluster, level) %>% summarize(auc = mean(auc)) %>% ungroup()

xdata$order <- xdata$egad %>% filter(level == "sbj") %>% arrange(desc(auc))

xdata$egad %>% 
  session$graphingUtils$ggplot(aes(y = cluster, x = auc)) + 
  geom_point(aes(color = level), size = 4, alpha = 0.5) + 
  ylim(unique(rev(xdata$order$cluster[1:30]))) 

xdata$egad %>% session$graphingUtils$ggplot(aes(x = level, y = auc)) + geom_violin()

xdata$avgcoex <- xdata$main$avgCoex %>% 
  filter(!grepl("-", dataset)) %>% 
  group_by(cluster, level) %>% summarize(coex = mean(coex)) 

xdata$order <- xdata$egad %>% filter(level == "bk") %>% arrange(desc(auc)) 

xdata$order <- xdata$avgcoex %>% filter(level == "bk") %>% arrange(desc(coex))

xdata$avgcoex %>% 
  group_by(cluster, level) %>% summarize(coex = mean(coex)) %>% 
  session$graphingUtils$ggplot(aes(y = cluster, x = coex)) + 
  geom_point(aes(color = level), size = 4, alpha = 0.5) + 
  ylim(unique(rev(xdata$order$cluster[1:50]))) + 
  ggtitle("velmeshev")

xdata$dat <- xdata$main$avgCoex %>% 
  filter(!grepl("-", dataset)) %>% 
  group_by(level, cluster) %>% 
  summarize(coex = mean(coex)) %>% 
  ungroup()

xdata$dat %>% group_by(cluster) %>% summarize(coex = mean(coex)) %>% arrange(desc(coex)) -> x

xdata$dat %>% filter(level == "bk") %>% arrange(desc(coex)) -> x


xdata$dat %>% session$graphingUtils$ggplot(aes(y = cluster, x = coex, color = level)) + geom_point(size = 4) + ylim(rev(x$cluster[1:30]))


xdata$main$avgCoex %>% filter(dataset %in% c("rosmap-ihc", "rosmap-ihcres")) %>% 
  session$graphingUtils$ggplot(aes(y = cluster, x = coex, color = dataset)) + geom_point(size = 4) + ylim(rev(x$cluster[1:30]))

ydata <- list()

ydata$main <- read_rds(paste0(workspace$outputDir, "sc_clusters_integrity_in_bk_coefmats.rds"))

ydata$main %>% filter(dataset %in% c("velmeshev-psdbk", "velmeshev-psdbkres")) %>% 
  session$graphingUtils$ggplot(aes(y = cluster, x = coex, color = dataset)) + geom_point(size = 4) + ylim(rev(xdata$order$cluster))

ydata$main %>% filter(dataset %in% c("velmeshev", "velmeshev-mgpres")) %>% dplyr::select(cluster, dataset, coex) %>% 
  spread(dataset, coex) %>% 
  mutate(diff = `velmeshev-mgpres` - `velmeshev`) %>% 
  session$graphingUtils$ggplot(aes(x = diff, y = cluster)) + 
  geom_bar(stat = "identity") +
  ylim(rev(xdata$order$cluster)) 

ydata$main <- read_rds(paste0(workspace$outputDir, "sc_clusters_integrity_in_bk_coefmats.rds"))

xdata$main$egad %>% filter(dataset %in% c("rosmap-psdbk", "rosmap-psdbkres")) %>% 
  session$graphingUtils$ggplot(aes(y = cluster, x = auc, color = dataset)) + geom_point(size = 4) + ylim(rev(x$cluster[1:30]))

# XXXXXX

































