# 1. prepare the expression matrices; QCs, etc
# 2. prepare cell level co-expression networks
# 3. prepare the cross-sample level co-expression networks
# 4. prepare the bulk tissue co-expression networks

library(qlcMatrix)
library(EGAD)
library(GO.db)
library(WGCNA)
library(data.table)

source("GlobalSession.R"); session <- initGlobalSession()

workspace <- list()

workspace$workspaceDir <- paste0(session$WORKSPACES, "041_CH4_FINAL/")
workspace$miscDir <- paste0(workspace$workspaceDir, "misc/")
workspace$spaceDir <- "/space/scratch/ericchu/r_cache/041_CH4_FINAL/"
workspace$dataDir <- paste0(workspace$spaceDir, "data/")
workspace$outputDir <- paste0(workspace$spaceDir, "output/") # for storing data outputs from data setup
workspace$sdataDir <- paste0(workspace$outputDir, "sdatafiles/")

workspace$utils <- list()

# CLEAN
rm(list = setdiff(ls(), c("session", "workspace", "pdata", "fdata")))
gc()

workspace$utils$size <- function(obj) {
  format(object.size(obj), units = "auto")
}

workspace$utils$stdVec <- function(vals) {
  (vals - min(vals)) / sd(vals)
}

workspace$utils$readMatSparse <- function(file) {
  dat <- Matrix::readMM(file)
  return(dat)
}

workspace$utils$sparseToDense <- function(dat) {
  dat <- dat %>% as.matrix()
  return(dat)
}

workspace$utils$ctmatToExprmat <- function(ctmat) {
  
  tbl <- tibble(j = ctmat@j, x = ctmat@x)
  tbl <- tbl %>% group_by(j) %>% mutate(x_norm = (x / sum(x)) * 1e6) %>% ungroup()
  ctmat@x <- tbl$x_norm
  
  return(ctmat)
}

workspace$utils$ctmatToExprmat2 <- function(ctmat) {
  
  # for dense matrices
  ctmatClone <- ctmat
  libSizes <- ctmatClone %>% apply(2, sum)
  result <- t( t(ctmatClone) / libSizes) * 1e6 # convert the unit into CPM

  return(result)
}

workspace$utils$cleanCtmat <- function(ctmat, geneThr = 0.05, sampleThr = 0.05) {
  # geneThr = % of samples must express a gene for it to be valid
  # sampleThr = % of n_genes beyond which samples are filtered out
  # filter out genes that are detected in less than X% of all cells of this type
  
  # filter for genes detected in at least N cells
  
  nCellsTot <- ncol(ctmat)
  nCellsThr <- geneThr * nCellsTot
  
  genesNN <- tibble(i = ctmat@i) %>% arrange(i) %>% group_by(i) %>% summarize(nn = n())
  genesNN <- genesNN %>% mutate(i = i + 1) # index starts at 0
  genesNN <- genesNN %>% filter(nn > nCellsThr)
  
  genes <- rownames(ctmat)
  genesV <- genes[genesNN$i]
  
  genesNN$gene <- genesV
  
  ctmat <- ctmat[genesNN$gene, ] # update count matrix to enable filtering cells
  
  # filter for cells with N genes
  
  cellsNN <- tibble(j = ctmat@j) %>% arrange(j) %>% group_by(j) %>% summarize(nn = n())
  cellsNN <- cellsNN %>% mutate(j = j + 1) # index starts at 0
  
  cellsNN <- cellsNN %>% 
    mutate(perc = percent_rank(nn)) %>% 
    filter(perc < (1 - sampleThr), perc > sampleThr)

  cells <- colnames(ctmat)
  cellsV <- cells[cellsNN$j]
  
  cellsNN$cell <- cellsV
  
  ctmat <- ctmat[, cellsNN$cell] # update count matrix to enable filtering cells
  
  return(ctmat)
  
}

workspace$utils$cleanCtmat2 <- function(ctmat, geneThr = 0.05, sampleThr = 0.05) {
  
  # for dense matrices
  
  # geneThr = % of samples must express a gene for it to be valid
  # sampleThr = % of samples to filter out, based on number of genes detected

  # filter for genes
  
  nSamplesTot <- ncol(ctmat)
  nSamplesThr <- geneThr * nSamplesTot
  
  genesNN <- ctmat %>% 
    apply(1, function(vals) { sum(vals > 0) }) %>% 
    session$dataWrangler$vectorToTibble() %>% 
    dplyr::select(gene = variable, nn = value)
    
  genesNN <- genesNN %>% filter(nn > nSamplesThr) %>% arrange(gene)
  
  ctmat <- ctmat[genesNN$gene, ] # update count matrix to enable filtering cells
  
  # filter for samples
  
  samplesNN <- ctmat %>% 
    apply(2, function(vals) { sum(vals > 0) }) %>% 
    session$dataWrangler$vectorToTibble() %>% 
    dplyr::select(sample = variable, nn = value) %>% 
    mutate(perc = percent_rank(nn))
  
  samplesNN <- samplesNN %>% 
    filter(perc <= (1 - sampleThr), perc >= sampleThr) %>% 
    arrange(sample)
  
  ctmat <- ctmat[, samplesNN$sample] # update count matrix to enable filtering cells
  
  return(ctmat)
  
}

workspace$utils$aggregSbj <- function(exprmat, samples) {
  # samples is a tibble with 2 columns: [sample, sample_bk]
  # exprmat is the expression matrix where columns are samples
  exprmat <- exprmat[, samples$sample]
  
  sampleGrps <- samples %>% 
    dplyr::select(sample, sample_bk) %>% 
    mutate(value = 1) %>% 
    spread(sample_bk, value) %>% 
    as.data.frame()
  
  rownames(sampleGrps) <- sampleGrps[, 1]
  sampleGrps <- sampleGrps[, -1]
  sampleGrps <- sampleGrps %>% as.matrix()
  sampleGrps[is.na(sampleGrps)] <- 0
  
  exprmat <- exprmat[, rownames(sampleGrps)]
  exprmatBk <- exprmat %*% sampleGrps
  exprmatBk <- exprmatBk %>% as.matrix()
  
  # compute the average
  nCels <- sampleGrps %>% apply(2, sum)
  exprmatBk <- t( t(exprmatBk) / nCels)
  
  return(exprmatBk)
}

workspace$utils$quantileNormalize <- function(exprmat) {
  # script obtained from https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/#:~:text=In%20statistics%2C%20quantile%20normalization%20is,arithmetical%20mean)%20of%20the%20distributions.
  df_rank <- exprmat %>% apply(2, rank,ties.method = "min")
  df_sorted <- exprmat %>% apply(2, sort)
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean = df_mean)
  rownames(df_final) <- rownames(exprmat)
  
  return(as.matrix(df_final))
}

workspace$utils$logTrans <- function(exprmat) {
  log2(exprmat + 1)
}

workspace$utils$readCounts <- function(file) {
  dat <- fread(file, data.table = FALSE) %>% as.matrix()
  rownames(dat) <- dat[, 1]
  dat <- dat[, -1]
  return(dat)
}

workspace$utils$asNumeric <- function(mat) {
  nMat <- mat %>% as.numeric() %>% matrix(ncol = ncol(mat))
  dimnames(nMat) <- dimnames(mat)
  return(nMat)
}

workspace$utils$computeCoexmat <- function(exprmat) {
  exprmat <- exprmat[sort(rownames(exprmat)), sort(colnames(exprmat))]
  coexMat <- exprmat %>% t() %>% cor(method = "pearson")
  return(coexMat)
}

workspace$utils$normalizeCoexmat <- function(coexmat) {
  
  rankCoexmat <- matrix(nrow = nrow(coexmat), ncol = ncol(coexmat))
  rownames(rankCoexmat) <- rownames(coexmat)
  colnames(rankCoexmat) <- colnames(coexmat)
  
  lowerCoords <- which(lower.tri(rankCoexmat), arr.ind = TRUE) # use the lower tri indices as primary
  upperCoords <- lowerCoords[, c("col", "row")] # figure out the matching coordinates to "flip it back"
  colnames(upperCoords) <- c("row", "col")
  
  currCoexVec <- coexmat[lowerCoords]
  currCoexVec[is.na(currCoexVec)] <- currCoexVec %>% median(na.rm = TRUE)
  
  coexRanks <- currCoexVec %>% rank(ties.method = "average")
  coexRanks <- coexRanks / max(coexRanks)
  
  rankCoexmat[lowerCoords] <- coexRanks
  rankCoexmat[upperCoords] <- coexRanks
  diag(rankCoexmat) <- 1
  
  return(rankCoexmat)
}


workspace$utils$coexMatToTbl <- function(coexMat) {
  tibble(gene_a = rownames(coexMat)[row(coexMat)[lower.tri(coexMat)]], 
         gene_b = colnames(coexMat)[col(coexMat)[lower.tri(coexMat)]], 
         cor_coef = coexMat[lower.tri(coexMat)]) %>% 
    mutate(pair_id = paste0(gene_a, ".", gene_b)) %>%
    dplyr::select(pair_id, everything())
}

workspace$utils$computePercMat <- function(coexMat) {
  coexMatCopy <- coexMat
  diag(coexMatCopy) <- NA
  coexVec <- coexMatCopy %>% as.vector()
  percMat <- coexVec %>% percent_rank() %>% matrix(ncol = ncol(coexMat))
  rownames(percMat) <- rownames(coexMat)
  colnames(percMat) <- colnames(coexMat)
  diag(percMat) <- 1
  return(percMat)
}

workspace$utils$computeLnksMat <- function(percMat, thresh = 0.99) {
  lnksMat <- percMat
  lnksMat[lnksMat >= thresh] <- 1
  lnksMat[lnksMat < thresh] <- 0
  return(lnksMat)
}

workspace$utils$vectorize <- function(coexMat) {
  coexMat[lower.tri(coexMat, diag = FALSE)]
}

workspace$utils$getCoords <- function(coexMat) {
  geneA <- rownames(coexMat)[row(coexMat)[lower.tri(coexMat, diag = FALSE)]]
  geneB <- colnames(coexMat)[col(coexMat)[lower.tri(coexMat, diag = FALSE)]]
  tbl <- tibble(gene_a = geneA, gene_b = geneB)
  return(tbl)
}

workspace$utils$getPairIds <- function(coexMat) {
  geneA <- rownames(coexMat)[row(coexMat)[lower.tri(coexMat)]]
  geneB <- colnames(coexMat)[col(coexMat)[lower.tri(coexMat)]]
  pairId <- paste0(geneA, ".", geneB)
  return(pairId)
}

workspace$utils$computeAggScCoex <- function(samples, exprmat) {
  
  # samples must have columns: sample, sample_bk
  
  samplesBk <- samples$sample_bk %>% unique() %>% sort()
  names(samplesBk) <- samplesBk
  
  genes <- rownames(exprmat) %>% sort()
  
  rankCoexmat <- matrix(nrow = length(genes), ncol = length(genes))
  rownames(rankCoexmat) <- genes
  colnames(rankCoexmat) <- genes
  
  lowerCoords <- which(lower.tri(rankCoexmat), arr.ind = TRUE) # use the lower tri indices as primary
  upperCoords <- lowerCoords[, c("col", "row")] # figure out the matching coordinates to "flip it back"
  colnames(upperCoords) <- c("row", "col")
  
  print("Computing individual co-expression networks.")
  
  coexRanks <- rep(0, nrow(lowerCoords))
  
  count <- 1
  
  samplesBk %>% session$collectionUtils$foreach(function(currSamplesBk) {
    
    print(paste0("working on sample ", count, " of ", length(samplesBk)))
    
    currSamples <- samples %>% filter(sample_bk == currSamplesBk)
    currExprmat <- exprmat[genes, currSamples$sample] %>% t()
    currCoexmat <- currExprmat %>% corSparse()
    
    currCoexVec <- currCoexmat[lowerCoords]
    currCoexVec[is.na(currCoexVec)] <- currCoexVec %>% median(na.rm = TRUE)
    
    currCoexRanks <- currCoexVec %>% rank(ties.method = "min")
    currCoexRanks <- currCoexRanks / max(currCoexRanks)
    
    coexRanks <<- coexRanks + currCoexRanks
    
    count <<- count + 1
    
  }) 
  
  print("Computing average co-expression.")
  
  coexRanksAgg <- coexRanks / length(samplesBk)
  
  print("Normalizing co-expressions.")
  
  coexRanksAgg <- coexRanksAgg %>% rank(ties.method = "average")
  coexRanksAgg <- coexRanksAgg / max(coexRanksAgg)
  
  print("Converting to matrix form.")
  
  rankCoexmat[lowerCoords] <- coexRanksAgg
  rankCoexmat[upperCoords] <- coexRanksAgg
  diag(rankCoexmat) <- 1
  
  return(rankCoexmat)
}

workspace$utils$computeClusters <- function(coexmat, cutreeFunc) { 
  # cutreeFunc is a function that takes in a dendro and outputs cluster labels in numeric
  
  coexmatDist <- 1 - coexmat
  dendro <- hclust(as.dist(coexmatDist), method = "average") 
  
  clusters <- cutreeFunc(dendro, coexmatDist)

  clusterTbl <- tibble(gene_id = dendro$labels, cluster = paste0("id_", formatC(clusters, digits = 3, flag = "0")), cluster_int = clusters) 
  
  clusterFlat <- clusterTbl$cluster %>% unique() %>% sort()
  names(clusterFlat) <- clusterFlat
  clusterFlat <- clusterFlat %>% lapply(function(currCluster) { (clusterTbl %>% filter(cluster == currCluster))$gene_id %>% unique() %>% sort() })
  
  clusterMat <- clusterTbl %>% 
    dplyr::select(gene_id, cluster) %>% 
    mutate(in_set = 1) %>% 
    spread(cluster, in_set)
    
  clusterMat <- clusterMat %>% session$dataWrangler$setColAsRownames("gene_id")
  clusterMat <- clusterMat %>% as.matrix()
  clusterMat[is.na(clusterMat)] <- 0
  
  return(list(dendro = dendro, clusters = clusters, clusterTbl = clusterTbl, clusterFlat = clusterFlat, clusterMat = clusterMat))
}

workspace$utils$computeNodeDegrees <- function(lnkmat) {
  lnkmatClone <- lnkmat[sort(rownames(lnkmat)), sort(colnames(lnkmat))]
  diag(lnkmatClone) <- 0
  
  degrees <- rowSums(lnkmatClone) + colSums(lnkmatClone)
  degrees <- degrees %>% session$dataWrangler$vectorToTibble() %>% dplyr::select(gene_id = variable, node_degree = value)
  
  return(degrees)
}

workspace$utils$computeAvgCoex <- function(coexmat, clusterFlats) {
  # clusters is a list<cluster> where cluster is a vector of genes
  names(clusterFlats) %>% lapply(function(currClusterId) {
    
    currCluster <- clusterFlats[[currClusterId]]
    
    genesV <- currCluster %>% intersect(rownames(coexmat)) %>% sort()
    currmat <- coexmat[genesV, genesV]
    vec <- currmat[lower.tri(currmat, diag = FALSE)]
    
    tibble(cluster = currClusterId, n_gene = length(genesV), n_lnks = length(vec), coex = mean(vec))
  }) %>% session$dataWrangler$rbind()
}

workspace$utils$computeAvgCoexForEach <- function(coexmats, clusterFlats) {
  
  # coexmats can be list of files paths pointing matrix rds files or an actual list of matrices
  
  if (is.character(coexmats)) {
    
    coexmats %>% mclapply(function(coexmat) {
      workspace$utils$computeAvgCoex(read_rds(coexmat), clusterFlats)
    }, mc.cores = min(5, length(coexmats))) %>% 
      session$collectionUtils$lapplyWithName(function(file, tbl) { tbl %>% mutate(coexfile = file) }, verbose = FALSE) %>% 
      session$dataWrangler$rbind() 
    
  } else {
    
    coexmats %>% mclapply(function(coexmat) {
      workspace$utils$computeAvgCoex(coexmat, clusterFlats)
    }, mc.cores = min(5, length(coexmats))) %>% 
      session$collectionUtils$lapplyWithName(function(file, tbl) { tbl %>% mutate(coexfile = file) }, verbose = FALSE) %>% 
      session$dataWrangler$rbind() 
    
  }
  
}



workspace$utils$computeEgad <- function(coexmat, clusterMat) {
  
  Vgenes <- intersect(rownames(clusterMat), rownames(coexmat))
  clusterMatClone <- clusterMat[Vgenes, ]
  # 
  # diffGenes <- rownames(coexmat) %>% setdiff(Vgenes) 
  # 
  # if (length(diffGenes) > length(Vgenes)) { 
  #   diffGenes <- diffGenes %>% sample(length(Vgenes))
  # }
  # clusterMatClone[diffGenes, ] <- 0
  # 
  nGenes <- clusterMatClone %>% colSums() %>%
    session$dataWrangler$vectorToTibble() %>% dplyr::select(cluster = variable, n_gene = value)
  
  EGAD::run_GBA(coexmat, clusterMatClone, max = 10000, min = 5)[[1]] %>% 
    session$dataWrangler$setRownameAsColumn("cluster") %>% 
    left_join(nGenes, by = "cluster")
    
}

workspace$utils$computeEgadForEach <- function(coexmats, clusterMat) {
  
  # coexmats can be list of files paths pointing matrix rds files or an actual list of matrices
  
  if (is.character(coexmats)) {
    
    coexmats %>% mclapply(function(coexmat) {
      workspace$utils$computeEgad(read_rds(coexmat), clusterMat)
    }, mc.cores = min(5, length(coexmats))) %>% 
      session$collectionUtils$lapplyWithName(function(file, tbl) { tbl %>% mutate(coexfile = file) }, verbose = FALSE) %>% 
      session$dataWrangler$rbind() 

  } else {
    
    coexmats %>% mclapply(function(coexmat) {
      workspace$utils$computeEgad(coexmat, clusterMat)
    }, mc.cores = min(5, length(coexmats))) %>% 
      session$collectionUtils$lapplyWithName(function(file, tbl) { tbl %>% mutate(coexfile = file) }, verbose = FALSE) %>% 
      session$dataWrangler$rbind() 
    
  }

}










workspace$utils$computeDensity <- function(coexmat, clusterFlats) {
  names(clusterFlats) %>% lapply(function(currClusterId) {
    
    currCluster <- clusterFlats[[currClusterId]]
    
    genesV <- currCluster %>% intersect(rownames(coexmat)) %>% sort()
    lnkmat <- coexmat[genesV, genesV]
    
    lnkmat[which(lnkmat >= 0.99)] <- 1
    lnkmat[which(lnkmat < 0.99)] <- 0
    
    
    subnetwork <- lnkmat[lower.tri(lnkmat, diag = FALSE)]
    subnetwork[is.na(subnetwork)] <- 0
    
    nLinks = sum(subnetwork)
    nPairs = length(subnetwork)
    
    tibble(cluster = currClusterId, n_links = nLinks, n_pairs = nPairs, density = nLinks / nPairs) 
    
  }) %>% session$dataWrangler$rbind()
}


workspace$utils$plotDendro <- function(clusters) {
  plotDendroAndColors(clusters$dendro, 
                      colors = labels2colors(clusters$clusters), 
                      rowText = clusters$clusters, 
                      dendroLabels = FALSE)
}

workspace$utils$fitCCVModels <- function(ctpMat, exprmat) {
  
  genes <- rownames(exprmat) 
  names(genes) <- genes
  
  print(paste0("fitting linear models for ", length(genes), " genes..."))
  
  lms <- genes %>% session$collectionUtils$lapply(function(currGene) {
    exprBk <- exprmat[currGene, ]
    covariates <- ctpMat[names(exprBk), ] %>% as.data.frame()
    lm(exprBk ~ ., covariates)
  }, verbose = FALSE)
  
  print(paste0("extracting lm stats..."))
  
  lmStats <- do.call(rbind, lms %>% session$collectionUtils$lapplyWithName(function(currGene, currLm) {
    lmSummary <- summary(currLm)
    fstat_pval <- pf(lmSummary$fstatistic[1], lmSummary$fstatistic[2], lmSummary$fstatistic[3], lower.tail = FALSE)
    tibble(gene_id = currGene, rsqr = lmSummary$r.squared, pvalue = fstat_pval) # extract the p-value as well?
  }, verbose = FALSE)) %>% 
    mutate(qvalue = p.adjust(pvalue, method = "fdr")) %>% 
    mutate(logpvalue = -log10(pvalue), 
           logqvalue = -log10(qvalue))
  
  print(paste0("extracting coef stats..."))
  
  coefStats <- do.call(rbind, lms %>% session$collectionUtils$lapplyWithName(function(currGene, currLm) {
    lmSummary <- summary(currLm)
    coefTbl <- lmSummary$coefficients %>% as_tibble(rownames = "cell_type")
    names(coefTbl) <- c("cell_type", "beta", "se", "t_stat", "pvalue")
    coefTbl <- coefTbl %>% mutate(gene_id = currGene) %>% dplyr::select(gene_id, everything())
  }, verbose = FALSE)) %>% 
    filter(cell_type != "(Intercept)") %>% 
    mutate(qvalue = p.adjust(pvalue, method = "fdr")) %>% 
    mutate(logpvalue = -log10(pvalue), 
           logqvalue = -log10(qvalue))
  
  celltypes <- coefStats$cell_type %>% unique() %>% sort()
  names(celltypes) <- celltypes 
  
  print("computing rsqr_indep...")
  
  rsqrsIndep <- do.call(rbind, genes %>% session$collectionUtils$lapply(function(currGene) {
    celltypes %>% sapply(function(currCelltype) {
      ctp <- ctpMat[, currCelltype]
      expr <- exprmat[currGene, rownames(ctpMat)]
      cor(ctp, expr) ^ 2
    }) %>% t() %>% as_tibble()
  }, verbose = FALSE)) %>% 
    mutate(gene_id = genes) %>% 
    gather(cell_type, rsqr_indep, -gene_id) %>% 
    dplyr::select(gene_id, everything())
  
  coefStats <- coefStats %>% left_join(rsqrsIndep, by = c("gene_id", "cell_type"))
  
  print("extracting residual exprmat...")
  
  exprmatRes <- lms %>% sapply(function(currLm) { currLm$residuals } ) %>% t()
  
  output <- list()
  
  output$ctpMat <- ctpMat
  output$stats$lmStats <- lmStats
  output$stats$coefStats <- coefStats
  output$exprmats$orig <- exprmat
  output$exprmats$residual <- exprmatRes
  
  return(output)  
}


workspace$utils$computeMgps <- function(mkrsFlat, exprmat) { 
  
  exprmatCopy <- exprmat %>% session$dataWrangler$setRownameAsColumn("gene_id")
  mgp <- markerGeneProfile::mgpEstimate(exprmatCopy, mkrsFlat, geneColName = "gene_id", geneTransform = NULL) 
  
  mgp$main <- mgp$estimates %>% sapply(function(vals) { vals })
  
  mgp$mkrs$flat <- mgp$usedMarkerExpression %>% lapply(function(exprmat) { rownames(exprmat) %>% unique() %>% sort() })
  mgp$mkrs$tbl <- mgp$mkrs$flat %>% session$collectionUtils$lapplyWithName(function(celltype, mkrs) { tibble(cell_type = celltype, gene_id = mkrs) }) %>% session$dataWrangler$rbind()
  mgp$mkrs$mat <- mgp$mkrs$tbl %>% mutate(inset = 1) %>% spread(cell_type, inset) %>% session$dataWrangler$fillNa(colNames = names(mgp$estimates), value = 0) %>% session$dataWrangler$setColAsRownames("gene_id")
  
  return(mgp)
  
  # saving the code below to do it via PCA just in case
  # actually Ogan's function does exactly the same thing so I shouldn't need to reimplement any of it. 
  # fdata$pcas <- fdata$mkrs$mkrsFlat %>% lapply(function(currMkrs) {
  #   currMkrs <- currMkrs %>% intersect(fdata$mkrsV)
  #   fdata$exprmatMkrs[currMkrs, ] %>% t() %>% prcomp(center = TRUE, scale = TRUE)
  # })
  
  # fdata$samples <- fdata$exprdat$samples %>% arrange(sample)
  
  # fdata$mgps <- do.call(rbind, fdata$pcas %>% session$collectionUtils$lapplyWithName(function(currName, currPca) {
  #   currPca$x[, "PC1"] %>% as_tibble(rownames = "sample") %>% 
  #     dplyr::select(sample, mgp = value) %>% 
  #     mutate(cell_type = currName)
  # }))
  
  # fdata$mgps <- fdata$mgps %>% 
  #   spread(cell_type, mgp) %>% 
  #   session$dataWrangler$setColAsRownames("sample")
    
}

workspace$utils$computeRSqr <- function(exprmatObs, exprmatCovs) {
  
  genes <- exprmatObs %>% rownames() %>% sort()
  names(genes) <- genes 
  
  samples <- exprmatObs %>% colnames() %>% sort() # important to sort so correlations are computed correctly
  names(samples) <- samples
  
  rsqr <- genes %>% sapply(function(currGene) {
    
    if (class(exprmatCovs) == "list") {
      currRsqrs <- exprmatCovs %>% sapply(function(currExprmatCov) {
        cor(exprmatObs[currGene, samples], currExprmatCov[currGene, samples]) ^ 2
      })
      names(currRsqrs) <- paste0("rsqr_", names(currRsqrs))
    } else {
      currRsqrs <- cor(exprmatObs[currGene, samples], exprmatCovs[currGene, samples]) ^ 2
      names(currRsqrs) <- "rsqr"
    }
    
    return(currRsqrs)
    
  }) %>% t() %>% 
    as_tibble(rownames = "gene_id")
  
  return(rsqr)
}


workspace$utils$computeLms <- function(exprmatBk, exprmatsPsdbk, verbose = TRUE) {
  
  samples <- Reduce(intersect, c(exprmatsPsdbk, list(bk = exprmatBk)) %>% lapply(colnames)) %>% sort()
  genes <- Reduce(intersect, c(exprmatsPsdbk, list(bk = exprmatBk)) %>% lapply(rownames)) %>% sort()
  names(genes) <- genes
  
  # filter for the right samples & genes and sorted
  
  currExprmatBk <- exprmatBk[genes, samples]
  currExprmatsPsbk <- exprmatsPsdbk %>% lapply(function(exprmat) { exprmat[genes, samples] })
  
  # assume that the dimensions match
  genes <- currExprmatBk %>% rownames()
  names(genes) <- genes
  
  lms <- genes %>% session$collectionUtils$lapply(function(currGene) {
    exprBk <- currExprmatBk[currGene, ]
    covariates <- currExprmatsPsbk %>% sapply(function(exprmat) { exprmat[currGene, ] }) %>% as.data.frame()
    lm(exprBk ~ ., covariates)
  }, verbose = verbose)
  
  lmStats <- do.call(rbind, lms %>% session$collectionUtils$lapplyWithName(function(currGene, currLm) {
    lmSummary <- summary(currLm)
    fstat_pval <- pf(lmSummary$fstatistic[1], lmSummary$fstatistic[2], lmSummary$fstatistic[3], lower.tail = FALSE)
    tibble(gene_id = currGene, rsqr = lmSummary$r.squared, pvalue = fstat_pval, aic = extractAIC(currLm)[2]) # compute AIC here as well
  }, verbose = verbose)) %>% 
    mutate(qvalue = p.adjust(pvalue, method = "fdr")) %>% 
    mutate(logpvalue = -log10(pvalue), 
           logqvalue = -log10(qvalue))
  
  return(list(models = lms, stats = lmStats))
  
}

workspace$utils$fitModelsSimple <- function(exprmatBk, exprmatsPsdbk) {
  
  samples <- Reduce(intersect, c(exprmatsPsdbk, list(bk = exprmatBk)) %>% lapply(colnames)) %>% sort()
  genes <- Reduce(intersect, c(exprmatsPsdbk, list(bk = exprmatBk)) %>% lapply(rownames)) %>% sort()
  names(genes) <- genes
  
  # filter for the right samples & genes and sorted
  
  currExprmatBk <- exprmatBk[genes, samples]
  currExprmatsPsbk <- exprmatsPsdbk %>% lapply(function(exprmat) { exprmat[genes, samples] })
  
  # now, per gene, fit linear model adding up all cell type signals
  
  lms <- workspace$utils$computeLms(currExprmatBk, currExprmatsPsbk, verbose = FALSE)
  
  return(lms$stats)
  
}



workspace$utils$fitModels <- function(exprmatBk, exprmatsPsdbk) {
  
  # filter for genes with at least some expression
  exprmatBk <- exprmatBk[ (apply(exprmatBk, 1, sum) > 1),  ]
  exprmatsPsdbk <- exprmatsPsdbk %>% lapply(function(exprmat) { exprmat[ (apply(exprmat, 1, sum) > 1), ] }) 
  
  samples <- Reduce(intersect, c(exprmatsPsdbk, list(bk = exprmatBk)) %>% lapply(colnames)) %>% sort()
  genes <- Reduce(intersect, c(exprmatsPsdbk, list(bk = exprmatBk)) %>% lapply(rownames)) %>% sort()
  names(genes) <- genes
  
  # filter for the right samples & genes and sorted
  
  currExprmatBk <- exprmatBk[genes, samples]
  currExprmatsPsbk <- exprmatsPsdbk %>% lapply(function(exprmat) { exprmat[genes, samples] })
  
  # now, per gene, fit linear model adding up all cell type signals
  
  lms <- workspace$utils$computeLms(currExprmatBk, currExprmatsPsbk)
  
  exprmatRes <- lms$models %>% sapply(function(currLm) { currLm$residuals } ) %>% t()
  
  # extract the coefficient stats
  
  coefStats <- do.call(rbind, lms$models %>% session$collectionUtils$lapplyWithName(function(currGene, currLm) {
    lmSummary <- summary(currLm)
    coefTbl <- lmSummary$coefficients %>% as_tibble(rownames = "cell_type")
    names(coefTbl) <- c("cell_type", "beta", "se", "t_stat", "pvalue")
    coefTbl <- coefTbl %>% mutate(gene_id = currGene) %>% dplyr::select(gene_id, everything())
  })) %>% 
    filter(cell_type %in% names(currExprmatsPsbk)) %>% 
    mutate(qvalue = p.adjust(pvalue, method = "fdr")) %>% 
    mutate(logpvalue = -log10(pvalue), 
           logqvalue = -log10(qvalue))
  
  
  lmsIndep <- names(currExprmatsPsbk) %>% session$collectionUtils$lapply(function(currCelltype) {
    print(paste0("working on independent models for: ", currCelltype))
    currLms <- workspace$utils$computeLms(currExprmatBk, currExprmatsPsbk[currCelltype])
    currLms$stats %>% 
      mutate(cell_type = currCelltype) %>% 
      dplyr::select(gene_id, cell_type, rsqr_indep = rsqr, pvalue_indep = pvalue, qvalue_indep = qvalue, aic_indep = aic) 
  }) %>% session$dataWrangler$rbind()
  
  coefStats <- coefStats %>% left_join(lmsIndep, by = c("gene_id", "cell_type"))
  
  return(list(lmStats = lms$stats, coefStats = coefStats, dims = list(samples = samples, genes = genes), 
              exprmats = list(orig = currExprmatBk, residual = exprmatRes)))
}


workspace$utils$fmtCelltypes <- function(strs) { 
  strs %>% sapply(function(str) {
    if (str == "opc") { "OPC" }
    else { str %>% str_to_sentence() }
  })
}

workspace$utils$fmtDataset <- function(strs) { 
  strs %>% sapply(function(str) {
    if (str == "rosmap") { "ROSMAP" }
    else { str %>% str_to_sentence() }
  })
}







