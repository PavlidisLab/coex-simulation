source("GlobalSession.R"); session <- initGlobalSession()

library(mvtnorm)

workspace <- list()

workspace$workspaceDir <- paste0(session$WORKSPACES, "042_CH3_FINAL/")
workspace$dataDir <- paste0(workspace$workspaceDir, "data/")
workspace$cacheLocDir <- paste0(workspace$workspaceDir, "cache/") # for writing out stuff
workspace$cacheSpaceDir <- "/space/scratch/ericchu/r_cache/042_CH3_FINAL/" # for large & generated files. 
workspace$outputDir <- paste0(workspace$cacheSpaceDir, "output/") # for storing data outputs from data setup

workspace$utils <- list()
workspace$data <- list()

source(paste0(workspace$workspaceDir, "the_simulator.R"))
workspace$initSimulator <- initSimulator
workspace$api <- workspace$initSimulator()$api

# CLEAN
rm(list = setdiff(ls(), c("session", "workspace", "pdata", "fdata")))
gc()

# utils implementation 

workspace$utils$size <- function(obj) {
  format(object.size(obj), units = "auto")
}

workspace$utils$readCounts <- function(file) {
  dat <- fread(file, data.table = FALSE)
  rownames(dat) <- dat[, 1]
  dat <- dat[, -1]
  dat <- dat %>% as.matrix()
  return(dat)
}

workspace$utils$normalizeCtmat <- function(ctmat) {
  
  # normFactors <- ctmat %>% calcNormFactors() # TODO GET RID OF THIS. 
  libSizes <- ctmat %>% apply(2, sum)
  
  # enforce the same ordering
  # libSizes <- libSizes[names(normFactors)]
  
  effLibSizes <- libSizes #* normFactors
  
  result <- t( t(ctmat) / effLibSizes) * 1e6 # convert the unit into CPM
  
  return(result)
}

workspace$utils$cleanCtmat <- function(ctmat, geneThr = 0.1, sampleThr = 2) {
  # geneThr = % of samples must express a gene for it to be valid
  # sampleThr = SD of n_genes beyond which samples are filtered out
  # filter out genes that are detected in less than X% of all cells of this type
  
  genesValid <- ctmat %>% apply(1, function(vals) { length(which(vals > 0)) }) %>% 
    as_tibble(rownames = "gene_id") %>% 
    dplyr::select(gene_id, n_cel = value)
  
  genesValid <- genesValid %>% mutate(frac_cel = n_cel / ncol(ctmat)) %>% filter(frac_cel >= geneThr) # 10%
  
  ctmat <- ctmat[genesValid$gene_id, ] # update count matrix to enable filtering cells
  
  samplesValid <- ctmat %>% apply(2, function(vals) { length(which(vals > 0)) }) %>% 
    as_tibble(rownames = "sample") %>% 
    dplyr::select(sample, n_gene = value)
  
  nGenes <- list(m = mean(samplesValid$n_gene), sd = sd(samplesValid$n_gene)) 
  
  samplesValid <- samplesValid %>% 
    filter(n_gene > (nGenes$m - (nGenes$sd * sampleThr)),
           n_gene < (nGenes$m + (nGenes$sd * sampleThr))) 
  
  return(list(vGenes = sort(genesValid$gene_id), vSamples = sort(samplesValid$sample)))
  
}

workspace$utils$computeCoexmat <- function(exprmat) {
  exprmat <- exprmat[sort(rownames(exprmat)), sort(colnames(exprmat))]
  coexMat <- exprmat %>% t() %>% cor(method = "pearson")
  return(coexMat)
}

workspace$utils$coexToTbl <- function(coexMat) {
  tibble(gene_a = rownames(coexMat)[row(coexMat)[upper.tri(coexMat)]], 
         gene_b = colnames(coexMat)[col(coexMat)[upper.tri(coexMat)]], 
         cor_coef = coexMat[upper.tri(coexMat)]) %>% 
    mutate(pair_id = paste0(gene_a, ".", gene_b)) %>%
    dplyr::select(pair_id, everything())
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


workspace$utils$aggregSbj <- function(exprmat, samples) {
  # samples is a tibble with 2 columns: [sample, subject]
  # cpm is the expression matrix where columns are samples
  
  samples <- samples %>% mutate(weight = 1)
  
  sampleGrps <- samples %>% 
    dplyr::select(sample, subject, weight) %>% 
    group_by(subject, sample) %>% 
    summarize(value = sum(weight)) %>% 
    ungroup()
  
  sampleGrps <- sampleGrps %>% 
    spread(subject, value) %>% 
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

workspace$utils$computeRSqr <- function(exprmatObs, exprmatCovs) {
  
  genes <- exprmatObs %>% rownames() %>% sort()
  names(genes) <- genes 
  
  samples <- exprmatObs %>% colnames() %>% sort() # important to sort so correlations are computed correctly
  names(samples) <- samples
  
  if (any(class(exprmatCovs) == "list")) {
    genes %>% sapply(function(currGene) {
      currRsqrs <- exprmatCovs %>% sapply(function(currExprmatCov) {
        cor(exprmatObs[currGene, samples], currExprmatCov[currGene, samples]) ^ 2
      })
      names(currRsqrs) <- paste0("rsqr_", names(currRsqrs))
      return(currRsqrs)
    }) %>% t() %>% session$dataWrangler$setRownameAsColumn("gene_id")
  } else {
    currRsqrs <- genes %>% sapply(function(currGene) {
      cor(exprmatObs[currGene, samples], exprmatCovs[currGene, samples]) ^ 2
    })
    currRsqrs <- currRsqrs %>% session$dataWrangler$vectorToTibble() %>% dplyr::select(gene_id = variable, rsqr = value)
    return(currRsqrs)
  }
  
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

workspace$utils$getExprGenes <- function(exprmat, gene_a, gene_b) {
  exprtbl <- exprmat[c(gene_a, gene_b), ] %>% t()
  colnames(exprtbl) <- c("gene_a", "gene_b")
  exprtbl <- exprtbl %>% session$dataWrangler$setRownameAsColumn("cell")
}


workspace$utils$computeGamMM <- function(paramMat) {
  
  paramMat %>% apply(1, function(vals) {
    a <- vals["shape"]
    s <- vals["scale"]
    m <- a * s
    v <- a * (s^2)
    m <- unname(m)
    v <- unname(v)
    result <- c(mean = m, variance = v)
    return(result)
  }) %>% t()
  
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













