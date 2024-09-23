
pdata <- list()

pdata$celltypes <- c("excitatory", "inhibitory", "opc", "oligodendrocyte", "astrocyte", "microglia")
names(pdata$celltypes) <- pdata$celltypes

pdata$ctprofiles <- read_rds(paste0(workspace$outputDir, "sc_stats_ctprofiles.rds"))

pdata$genes <- readRDS(paste0(workspace$outputDir, "genes_metadata.rds")) # genes

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

pdata$files$lnkmats <- tibble(outfile = list.files(workspace$outputDir)) %>% 
  filter(grepl("lnkmats", outfile)) %>% 
  mutate(outfile_name = gsub(".rds", "", outfile)) %>% 
  mutate(strs = strsplit(outfile_name, "_")) %>% 
  mutate(level = strs %>% sapply(function(currStr) { currStr[1] })) %>% 
  mutate(type = strs %>% sapply(function(currStr) { currStr[2] })) %>% 
  mutate(cell_type = strs %>% sapply(function(currStr) { currStr[3] })) %>% 
  dplyr::select(outfile, level, type, cell_type)


# +++++++++++++++++++++++=
# now, let's cluster the xSubject level networks, in the same way

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
  
  files <- pdata$files$coexmats %>% filter(cell_type == currCelltype, level == "sbj")
  
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

fdata$coexmats %>% saveRDS(paste0(workspace$outputDir, paste0("sbj_coexmats_consensus.rds")))
# fdata$coexmats <- read_rds(paste0(workspace$outputDir, paste0("sbj_coexmats_consensus.rds")))

# +++++++++++++++++++++
# let's cluster the consensus stuff

fdata <- list(); gc()

fdata$coexmats <- read_rds(paste0(workspace$outputDir, paste0("sbj_coexmats_consensus.rds")))

fdata$clusters <- fdata$coexmats %>% mclapply(function(currCoexmat) {
  workspace$utils$computeClusters(currCoexmat, function(dendro, distMat) {
    dynamicTreeCut::cutreeDynamic(dendro, distM = distMat, cutHeight = 0.2, minClusterSize = 30)
  })
}, mc.cores = length(fdata$coexmats))

xdata <- list()

xdata$celltype <- "inhibitory"

plotDendroAndColors(fdata$clusters[[xdata$celltype]]$dendro, 
                    colors = labels2colors(fdata$clusters[[xdata$celltype]]$clusters), 
                    rowText = fdata$clusters[[xdata$celltype]]$clusters, 
                    dendroLabels = FALSE)

# ++++++++++++==== COMMIT
fdata$clusters %>% saveRDS(paste0(workspace$outputDir, paste0("sbj_clusters_objects.rds")))


# +++++++++++++++++
# characterize the number of clusters and their n_gene distributions across cell types

fdata <- list(); gc()

fdata$clusters <- read_rds(paste0(workspace$outputDir, paste0("sbj_clusters_objects.rds")))

# look at the number of genes per cluster, per cell type
fdata$nGenes <- fdata$clusters %>% session$collectionUtils$lapplyWithName(function(currCelltype, currCluster) {
  currCluster$clusterTbl %>% group_by(cluster) %>% summarize(n_gene = n()) %>% mutate(cell_type = currCelltype)
}) %>% session$dataWrangler$rbind()

xdata <- list()

xdata$celltype <- "inhibitory"

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
  saveRDS(paste0(workspace$outputDir, paste0("sbj_clusters_objects.rds")))


# ++++++++++++++++++++++++++==
# let's look at just the cluster quality in consensus networks to ensure quality and cell type specifciity

fdata <- list(); gc()

fdata$clusters <- read_rds(paste0(workspace$outputDir, paste0("sbj_clusters_objects.rds")))
# fdata$clustersGeneric <- read_rds(paste0(workspace$outputDir, paste0("sc_generic_clusters.rds")))

# fdata$clusterMats <- fdata$clusters$clusterMats %>% session$collectionUtils$lapply(function(currClusterMat) {
#   
#   genericMat <- fdata$clustersGeneric$clusterMat
#   genericMat <- genericMat[intersect(rownames(currClusterMat), rownames(genericMat)), ]
#   
#   genesDiff <- rownames(currClusterMat) %>% setdiff(rownames(genericMat))
#   genericMat[genesDiff, ] <- 0
#   genericMat <- genericMat[rownames(currClusterMat), ]
#   
#   return(currClusterMat %>% cbind(genericMat))
# })
# 
# fdata$clusterFlats <- fdata$clusters$clusterFlats %>% session$collectionUtils$lapply(function(currClusterFlat) {
#   genesValid <- currClusterFlat %>% unlist() %>% unname() %>% unique()
#   genericFlat <- fdata$clustersGeneric$clusterFlat %>% lapply(function(genes) { genes %>% intersect(genesValid) })
#   return(c(currClusterFlat, genericFlat))
# })

# fdata$coexmatFiles <- pdata$files$coexmats %>% filter(dataset != "consensus")
# fdata$coexmatFiles <- paste0(workspace$outputDir, fdata$coexmatFiles$outfile) %>% session$dataWrangler$attachNames(fdata$coexmatFiles$outfile)

fdata$coexmats <- read_rds(paste0(workspace$outputDir, paste0("sbj_coexmats_consensus.rds")))

fdata$egad <- fdata$clusters$clusterMats %>% session$collectionUtils$lapply(function(currClusterMat) {
  workspace$utils$computeEgadForEach(fdata$coexmats, currClusterMat) %>% dplyr::select(coex_ct = coexfile, everything())
}) %>% session$dataWrangler$rbind()

fdata$egad %>% session$graphingUtils$ggplot(aes(x = cluster, y = auc)) + geom_point(aes(color = coex_ct)) + session$graphingUtils$tiltX(angle = 90)

fdata$avgcoex <- fdata$clusters$clusterFlats %>% session$collectionUtils$lapplyWithName(function(cellType, currClusterFlat) {
  workspace$utils$computeAvgCoexForEach(fdata$coexmats, currClusterFlat) %>% dplyr::select(coex_ct = coexfile, everything())
}) %>% session$dataWrangler$rbind()

fdata$avgcoex %>% 
  session$graphingUtils$ggplot(aes(x = cluster, y = coex)) + geom_point(aes(color = coex_ct)) + session$graphingUtils$tiltX(angle = 90)


# ========== COMMIT
list(egad = fdata$egad, avgcoex = fdata$avgcoex) %>% 
  saveRDS(paste0(workspace$outputDir, paste0("sbj_clusters_integrity_across_celltype.rds")))
# =============

# ++++++++++++++++++++++++++==
# let's look at the quality of the clusters that were discovered at the sc level
# 1. egad, 2. average co-expression
# in sc, sbj consensus and component networks in the same cell type

fdata <- list(); gc()

fdata$clusters <- read_rds(paste0(workspace$outputDir, paste0("sbj_clusters_objects.rds")))

# fdata$clustersGeneric <- read_rds(paste0(workspace$outputDir, paste0("sc_generic_clusters.rds")))

# fdata$clusterMats <- fdata$clusters$clusterMats %>% session$collectionUtils$lapply(function(currClusterMat) {
#   
#   genericMat <- fdata$clustersGeneric$clusterMat
#   genericMat <- genericMat[intersect(rownames(currClusterMat), rownames(genericMat)), ]
#   
#   genesDiff <- rownames(currClusterMat) %>% setdiff(rownames(genericMat))
#   genericMat[genesDiff, ] <- 0
#   genericMat <- genericMat[rownames(currClusterMat), ]
#   
#   return(currClusterMat %>% cbind(genericMat))
# })
# 
# fdata$clusterFlats <- fdata$clusters$clusterFlats %>% session$collectionUtils$lapply(function(currClusterFlat) {
#   genesValid <- currClusterFlat %>% unlist() %>% unname() %>% unique()
#   genericFlat <- fdata$clustersGeneric$clusterFlat %>% lapply(function(genes) { genes %>% intersect(genesValid) })
#   return(c(currClusterFlat, genericFlat))
# })

fdata$coexmatFiles <- pdata$files$coexmats %>% filter(dataset != "consensus")
# fdata$coexmatFiles <- paste0(workspace$outputDir, fdata$coexmatFiles$outfile) %>% session$dataWrangler$attachNames(fdata$coexmatFiles$outfile)

fdata$coexmatsConsensus <- list(sc = read_rds(paste0(workspace$outputDir, paste0("sc_coexmats_consensus.rds"))), 
                                sbj = read_rds(paste0(workspace$outputDir, paste0("sbj_coexmats_consensus.rds")))) %>% 
  unlist(recursive = FALSE)

# EGAD

fdata$egadComponents <- fdata$clusters$clusterMats %>% session$collectionUtils$lapplyWithName(function(cellType, currClusterMat) {
  
  coexmatFiles <- fdata$coexmatFiles %>% filter(cell_type == cellType | level == "bk")
  coexmats <- paste0(workspace$outputDir, coexmatFiles$outfile) %>% session$dataWrangler$attachNames(coexmatFiles$outfile)
  
  workspace$utils$computeEgadForEach(coexmats, currClusterMat) %>% 
    left_join(coexmatFiles %>% dplyr::select(coexfile = outfile, level, dataset, coex_ct = cell_type), by = "coexfile") %>% 
    mutate(cluster_ct = cellType)
  
}) %>% session$dataWrangler$rbind()


fdata$egadConsensus <- fdata$clusters$clusterMats %>% session$collectionUtils$lapplyWithName(function(cellType, currClusterMat) {
  
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

fdata$avgcoexComponents <- fdata$clusters$clusterFlats %>% session$collectionUtils$lapplyWithName(function(cellType, currClusterFlat) {
  
  coexmatFiles <- fdata$coexmatFiles %>% filter(cell_type == cellType | level == "bk")
  coexmats <- paste0(workspace$outputDir, coexmatFiles$outfile) %>% session$dataWrangler$attachNames(coexmatFiles$outfile)
  
  workspace$utils$computeAvgCoexForEach(coexmats, currClusterFlat) %>% 
    left_join(coexmatFiles %>% dplyr::select(coexfile = outfile, level, dataset, coex_ct = cell_type), by = "coexfile") %>% 
    mutate(cluster_ct = cellType)
  
}) %>% session$dataWrangler$rbind()


fdata$avgcoexConsensus <- fdata$clusters$clusterFlats %>% session$collectionUtils$lapplyWithName(function(cellType, currClusterFlat) {
  
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
  saveRDS(paste0(workspace$outputDir, paste0("sbj_clusters_integrity_within_celltype.rds")))
# =============


# ++++++++++++++++++++++++++==
# OK - finally can start plotting some figures & do some writing
# first, characterize the number of clusters / cluster sizes at the single cell level

fdata <- list(); gc()

fdata$clusters <- read_rds(paste0(workspace$outputDir, paste0("sbj_clusters_objects.rds")))

fdata$objs <- fdata$clusters$clusterObjs

fdata$clusterFlat <- fdata$clusters$clusterFlats %>% unname() %>% unlist(recursive = FALSE)

fdata$clusterSize <- fdata$clusterFlat %>% session$collectionUtils$lapplyWithName(function(clusterName, genes) {
  tibble(cluster = clusterName, n_gene = length(genes))
}) %>% session$dataWrangler$rbind()

fdata$clusterSize <- fdata$clusterSize %>%
  mutate(cell_type = cluster %>% sapply(function(currStr) { currStr %>% str_extract("^[a-z]+") }))

# cluster sizes

fdata$clusterSize %>% 
  session$graphingUtils$ggplot(aes(y = cluster, x = n_gene), size = "SMALL") + 
  geom_bar(aes(fill = cell_type), stat = "identity")

fdata$clusterSize %>% session$graphingUtils$ggplot(aes(x = n_gene)) + geom_histogram(aes(fill = cell_type))

fdata$clusterSize$n_gene %>% mean()

# number of clusters

xdata <- list()

xdata$main <- fdata$clusterSize %>% group_by(cell_type) %>% summarize(n_cluster = n())

xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = cell_type, y = n_cluster)) + 
  geom_bar(stat = "identity") + 
  session$graphingUtils$tiltX(angle = 90) + 
  geom_text(aes(label = n_cluster), vjust = -1, size = 5) + 
  ylim(0, max(xdata$main$n_cluster) + 5)


# gene, across all networks

xdata <- list()

xdata$main <- fdata$clusters$clusterObjs %>% session$collectionUtils$lapplyWithName(function(currCelltype, currCluster) {
  currCluster$clusterTbl %>% mutate(in_cluster = cluster_int > 0) %>% mutate(cell_type = currCelltype)
}) %>% session$dataWrangler$rbind()

xdata$inGenes <- xdata$main %>% filter(in_cluster) %>% session$dataWrangler$extractColumn("gene_id") %>% unique()
xdata$inGenes %>% length() # number of genes in cluster

xdata$allGenes <- xdata$main$gene_id %>% unique()
xdata$allGenes %>% length()

# gene coverage

xdata <- list()

xdata$genesOut <- fdata$clusters$clusterObjs %>% sapply(function(currCluster) {
  currCluster$clusterTbl %>% filter(cluster_int == 0) %>% nrow()
})

xdata$genesIn <- fdata$clusterSize %>% group_by(cell_type) %>% summarize(n_gene_in = sum(n_gene))

xdata$genesIn$n_gene_out <- xdata$genesOut[xdata$genesIn$cell_type]

xdata$genesIn <- xdata$genesIn %>% mutate(total = n_gene_in + n_gene_out) %>% mutate(frac_in = n_gene_in / total)

xdata$main <- xdata$genesIn %>% dplyr::select(cell_type, n_gene_in, n_gene_out) %>% 
  reshape2::melt() %>% as_tibble() %>% mutate(variable = factor(variable, levels = c("n_gene_out", "n_gene_in")))

xdata$main %>% session$graphingUtils$ggplot(aes(x = cell_type, y = value)) + 
  geom_bar(aes(fill = variable), stat = "identity") +
  session$graphingUtils$tiltX(angle = 90)


# look at the number of genes per cluster, per cell type
fdata$nGenes <- fdata$objs %>% session$collectionUtils$lapplyWithName(function(currCelltype, currCluster) {
  currCluster$clusterTbl %>% group_by(cluster) %>% summarize(n_gene = n()) %>% mutate(cell_type = currCelltype)
}) %>% session$dataWrangler$rbind()

fdata$celltype <- pdata$celltypes[6]

plotDendroAndColors(fdata$objs[[fdata$celltype]]$dendro, 
                    colors = labels2colors(fdata$objs[[fdata$celltype]]$clusters), 
                    rowText = fdata$objs[[fdata$celltype]]$clusters, 
                    dendroLabels = FALSE)



# ++++++++++++++++++++++++++==
# here - look at how well single cell networks capture single cell clusters

fdata <- list(); gc()

fdata$main <- readRDS(paste0(workspace$outputDir, paste0("sbj_clusters_integrity_within_celltype.rds")))
fdata$main <- fdata$main %>% lapply(function(currTbl) { currTbl %>% filter(level == "sbj") })


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

xdata$celltype <- pdata$celltypes[2]

xdata$main <- fdata$main$avgCoex %>% filter(coex_ct == xdata$celltype)

xdata$consensus <- xdata$main %>% filter(dataset == "consensus")
xdata$component <- xdata$main %>% filter(dataset != "consensus")

xdata$componentMeans <- xdata$component %>% group_by(cluster) %>% summarize(coex = mean(coex))

xdata$component %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_density(aes(color = dataset, fill = dataset)) +
  facet_wrap(~dataset, ncol = 1) +
  xlim(0, 1) + 
  ggtitle(paste0("Average co-expression in component networks:\n", xdata$celltype))

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

xdata$celltype <- pdata$celltypes[1]

xdata$main <- fdata$main$egad %>% filter(cell_type == xdata$celltype)

xdata$consensus <- xdata$main %>% filter(dataset == "consensus")
xdata$component <- xdata$main %>% filter(dataset != "consensus")
xdata$componentMeans <- xdata$component %>% group_by(cluster) %>% summarize(auc = mean(auc))

xdata$component %>% 
  session$graphingUtils$ggplot(aes(x = auc)) + 
  geom_density(aes(color = dataset, fill = dataset)) +
  facet_wrap(~dataset, ncol = 1) +
  xlim(0, 1) + 
  ggtitle(paste0("Average EGAD performance in component networks:\n", xdata$celltype))

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


# ++++++++++++++++++++++++++==
# characterize how clusters are cell type specific

fdata <- list(); gc()

fdata$clusters <- read_rds(paste0(workspace$outputDir, paste0("sbj_clusters_objects.rds")))

fdata$main <- readRDS(paste0(workspace$outputDir, paste0("sbj_clusters_integrity_across_celltype.rds")))

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

fdata$clusters <- read_rds(paste0(workspace$outputDir, paste0("sbj_clusters_objects.rds")))

fdata$main <- readRDS(paste0(workspace$outputDir, paste0("sbj_clusters_integrity_across_celltype.rds")))

# for the generic clusters; pick the one with the highest number of genes

fdata$egad <- fdata$main$egad %>% group_by(cluster, coex_ct) %>% filter(n_gene == max(n_gene)) %>% summarize(auc = mean(auc)) %>% ungroup()

fdata$egad <- fdata$egad %>% 
  dplyr::select(cluster, coex_ct, auc) %>% spread(coex_ct, auc) %>% na.omit() %>% 
  session$dataWrangler$setColAsRownames("cluster") 

fdata$egadGeneric <- fdata$egad %>% 
  apply(1, function(x) { min(x, na.rm = TRUE) }) %>% session$dataWrangler$vectorToTibble() %>% 
  dplyr::select(cluster = variable, auc_generic = value) %>% 
  arrange(desc(auc_generic)) 

fdata$genericClusters <- fdata$egadGeneric %>% filter(auc_generic > 0.9) %>% arrange(cluster) # these are 27 clusters where the minimum auroc > 0.9

fdata$egad[fdata$egadGeneric$cluster[1:50], ] %>% 
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

# perform GO enrichment analysis on these 9 clusters to figure out the contents; one expectation is ribosomal protein genes; there may be others

fdata <- list()

fdata$clusters <- read_rds(paste0(workspace$outputDir, paste0("sbj_clusters_objects.rds")))

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































