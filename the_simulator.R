initSimulator <- function() {
  
  private <- list()
  
  private$utils <- list()
  
  private$utils$initThis <- function() {
    
    this <- list()
    
    this$genes <- NULL
    this$celltypes <- NULL
    
    # container for reference data, each entry is an expression count matrix from which the four components of the model are learned: 
    # 1. subject and 2. cell level marginal distributions and 3. subject and 4. cell level correlations
    # 5. additionally ctMat is the container for the CT expression profiles (define the cell type means per gene)
    
    this$refdat <- list()
    
    this$refdat <- list(exprmat = list(cel = NULL, sbj = NULL), 
                        coexPrograms = list(sbj = NULL, cel = NULL), 
                        cteprf = NULL) 
    
    # container for parameters for marginal distributions; list of N (per cell type) data matrices, each matrix has X genes and 2 columns with sbj/cel level marginal params
    this$mdParams <- list(cel = NULL, sbj = NULL)
    
    return(this)
  }
  
  
  private$utils$estNbParams <- function(exprmat) {
    
    # output a matrix with 2 columns: mu and phi, and G rows
    
    print(paste0("Estimating gene wise mu and phi (nb) for expression matrix of G = ", nrow(exprmat), " genes and S = ", ncol(exprmat), " samples."))
    
    if (nrow(exprmat) < 10) {
      phi <- private$utils$estPhi(exprmat)
    } else {
      phi <- estimateDisp(exprmat)$tagwise.dispersion
    }
    
    mu <- exprmat %>% apply(1, mean)
    
    # enforce the same order
    # phi <- phi[names(mu)]
    
    result <- cbind(mu, phi)
    
    return(result)
  }
  
  private$utils$muphiToPr <- function(muphi) {
    
    mu <- muphi[, "mu"]
    phi <- muphi[, "phi"]
    
    size <- 1 / phi
    p <- size / (size + mu) # use R parameterization for probability of success (instead of failure)
    r <- size
    
    result <- cbind(p = p, r = r)
    
    return(result)
  }
  
  
  private$utils$nbToGamma <- function(nbParams) {
    
    params <- private$utils$muphiToPr(nbParams)
    
    p <- params[, "p"]
    r <- params[, "r"]
    
    scale <- (1 - p) / p
    shape <- r
    
    result <- cbind(scale = scale, shape = shape)
    
    # hack: add some negligible baseline to scale; exactly 0 leads to NaN values in the simulation. 
    # result[, "scale"] <- result[, "scale"] + 1e-10
    
    return(result)
  }
  
  private$utils$computeGamMM <- function(paramMat) {
    
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
  
  private$utils$gamMMToParams <- function(mmMat) {
    
    mmMat %>% apply(1, function(vals) {
      m <- vals["mean"]
      v <- vals["variance"]
      a <- (m^2) / v
      s <- v / m 
      a <- unname(a)
      s <- unname(s)
      result <- c(shape = a, scale = s)
      return(result)
    }) %>% t()
    
  }
  
  private$utils$generateVals <- function(n, mdParams, corMat, idNames) {
    
    # TODO - move the non-coexpressed genes to the top because they should always be present
    
    # nCell -> n
    # currParams -> mdParams
    # currCormat -> corMat
    # idNames <- paste0(currCelltype, ".", currSbj, ".cel")
    
    # nBkSamples -> n
    # currParams -> mdParams
    # currSbjCor -> corMat
    # subjects -> idNames
    
    if (length(corMat) <= 0) {
    
      copula <- data.frame(matrix(nrow = n, ncol = 0))
      
    } else {
      
      # output - matrix where rows are genes and columns are samples/subjects, each cell is CPM expression
      copula <- rmvnorm(n = n, sigma = corMat, checkSymmetry = FALSE)  
      
    }
      
    if (length(idNames) == 1) {
      rownames(copula) <- paste0(idNames, "_", formatC(1:n, digits = ceiling(log10(n)), flag = "0"))
    } else {
      rownames(copula) <- idNames
    }
    
    colnames(copula) <- colnames(corMat)
    
    # standardize the data so that everything has an SD == 1
    # copula <- copula %>% apply(2, function(vals) { vals / sd(vals) }) # this step is redundant
    
    # now, convert to uniform distributions
    copula <- copula %>% apply(2, function(vals) { pnorm(vals) })
    
    # for now, cols are genes, rows are samples
    
    # here, go through each gene and transform the sample into the appropriate marginal distribution
    genesCopula <- colnames(copula) %>% session$dataWrangler$attachNames()
    result <- genesCopula %>% sapply(function(geneId) {
      currParam <- mdParams[geneId,]
      scale <- currParam[["scale"]]
      shape <- currParam[["shape"]]
      currSample <- copula[, geneId] %>% qgamma(shape = shape, scale = scale)
      return(currSample)
    }) %>% t() 
    
    # now, for the genes not included in the copula, sample directly using rgamma. 
    
    nonCorGenes <- rownames(mdParams) %>% setdiff(rownames(result))
    
    # in case there aren't any genes uncoexpressed, return right away. 
    if (length(nonCorGenes) <= 0) { return(result[rownames(mdParams),]) }
    
    result2 <- nonCorGenes %>% sapply(function(geneId) {
      currParam <- mdParams[geneId, ]
      scale <- currParam[["scale"]]
      shape <- currParam[["shape"]]
      currSample <- rgamma(n, shape = shape, scale = scale)
      return(currSample)
    }) %>% t()
    
    if (length(result > 0)) {
      colnames(result2) <- colnames(result)
      exprmat <- result %>% rbind(result2)
    } else {
      if (length(idNames) == 1) {
        colnames(result2) <- paste0(idNames, "_", formatC(1:n, digits = ceiling(log10(n)), flag = "0"))
      } else {
        colnames(result2) <- idNames
      }
      exprmat <- result2
    }
    
    exprmat <- exprmat[rownames(mdParams),] # reorder genes to what was specified in the md parameters
    
    return(exprmat)
  }
  
  
  # container for all public / main functions to be called by the client
  # all public functions should take in and return the a simulator object
  api <- list() 
  
  
  api$setRefDatMd <- function(this, exprmat, level = c("sbj", "cel")) {
    
    # samples tibble should provide cell type annotations
    
    # confirm that we have annotations for all the subjects
    # print out subjects & cell type tallies
    
    print(paste0("setting reference data for level: ", level, "..."))
    
    print(paste0("reference data matrix has g = ", nrow(exprmat), " genes and s = ", ncol(exprmat), " samples. "))
    
    this$refdat$exprmat[[level]] <- exprmat
    
    return(this)
  }
  
  api$setRefDatCteprf <- function(this, cteprf) {
    
    print(paste0("setting reference cteprf..."))
    
    print(paste0("reference data matrix has g = ", nrow(cteprf), " genes and n = ", ncol(cteprf), " cell types: ", paste0(colnames(cteprf), collapse = ", ")))
    
    this$refdat$cteprf <- cteprf
    
    return(this)
  }
  
  
  api$setRefDatCoexPrgms <- function(this, coexPrograms, level = c("sbj", "cel")) {
    
    # When two or more cell types share a coexProgram they are synchronized at the xSubject level
    # to achieve "independent coexpression" define separate coexPrograms
    # at the xCell level, this setting makes no differnece
    
    print(paste0("setting reference coexPrograms for level: ", level, "..."))
    
    # geneIds <- cormats[[1]] %>% rownames()
    
    print(paste0(length(coexPrograms), " coexPrograms for has been supplied."))
    
    coexPrograms %>% session$collectionUtils$lapplyWithName(function(id, program) {
      print(paste0("coexProgram: ", id, 
                   "; cell types: ", paste0(program$celltypes, collapse = ", "), 
                   "; n genes: ", length(program$genes), 
                   "; coex value = ", program$value))
    }, verbose = FALSE)
    
    this$refdat$coexPrograms[[level]] <- coexPrograms
    
    return(this)
  }
  
  
  api$fitMds <- function(this) {
    
    cteprf <- this$refdat$cteprf
    
    celltypes <- colnames(cteprf) %>% sort()
    names(celltypes) <- celltypes
    
    # fit mds per gene for subject level variability
    # use reference data, fit nb models
    # convert nb dispersion values into gamma ---- changed to directly fitting Gamma because the extra step just makes it too complicated. 
    # fit a trend to predict gamma variance from the mean
    # use the trend to generate estimates for each gene/cell type in cteprf to finalize a gamma distribution per gene
    
    sbjExprmat <- this$refdat$exprmat$sbj
    sbjGamMMs <- sbjExprmat %>% apply(1, function(vals) { c(mean = mean(vals), variance = var(vals)) }) %>% t() 
    sbjGamMMs <- log2(sbjGamMMs)
    sbjGamMMs <- sbjGamMMs %>% as_tibble(rownames = "gene_id")
    
    # sbjNbParams <- sbjExprmat %>% private$utils$estNbParams()
    # sbjGamParams <- sbjNbParams %>% private$utils$nbToGamma()
    # sbjGamMMs <- sbjGamParams %>% private$utils$computeGamMM()

    # sbjGamMMs %>%
    #   session$graphingUtils$ggplot(aes(x = mean, y = variance)) +
    #   geom_point()
    
    sbjMvModel <- lm(variance ~ poly(mean, 3, raw = TRUE), sbjGamMMs) # third degree polynomial
    # sbjMvModel <- loess(variance ~ mean, sbjGamMMs) # loess
    
    # let's check model fit? 
    # sbjGamMMs %>%
    #   mutate(fit = sbjMvModel$fitted) %>%
    #   session$graphingUtils$ggplot(aes(x = mean)) +
    #   geom_point(aes(y = variance)) +
    #   geom_line(aes(y = fit), color = "blue", size = 1)
    
    sbjErrorSd <- sbjMvModel$residuals %>% sd()
    
    print(paste0("fitting sbj gamma models for each of ", length(celltypes), " cell types: ", paste0(celltypes, collapse = ", "), "."))
    
    sbjGamMMsF <- celltypes %>% lapply(function(currCelltype) {
      
      params <- cteprf[, currCelltype] %>% as_tibble(rownames = "gene_id") %>% dplyr::select(gene_id, mean = value)
      params <- params %>% mutate(mean = log2(mean))
      
      params <- params %>% # generate variances for given genes
        mutate(variance = predict(sbjMvModel, params)) %>% 
        mutate(variance = variance + rnorm(length(variance), sd = sbjErrorSd))
      
      # params <- params %>% mutate(variance = pmax(variance, 2)) # enforce baseline variance
      # params <- params %>% mutate(mean = pmax(mean, 1)) # enforce CPM > 0
      paramsMat <- params %>% dplyr::select(mean, variance) %>% as.matrix()
      rownames(paramsMat) <- params$gene_id
      paramsMat <- (2^paramsMat)
      
      return(paramsMat)
    })
    
    # now convert the mean-var tibble to shape/scale paramtrization for gamma model definitions
    
    sbjGamParamsF <- sbjGamMMsF %>% lapply(private$utils$gamMMToParams)
    
    # COMMIT
    this$mdParams$sbj <- sbjGamParamsF
    
    # =====================
    # fit mds per gene for subject level variability; for each cell type
    
    celExprmat <- this$refdat$exprmat$cel
    
    celGamMMs <- celExprmat %>% apply(1, function(vals) { c(mean = mean(vals), variance = var(vals)) }) %>% t() 
    celGamMMs <- log2(celGamMMs)
    celGamMMs <- celGamMMs %>% as_tibble(rownames = "gene_id")
    
    # celNbParams <- celExprmat %>% private$utils$estNbParams()
    # celGamParams <- celNbParams %>% private$utils$nbToGamma()
    # celGamMMs <- celGamParams %>% private$utils$computeGamMM()
    
    # celGamMMs %>%
    #   session$graphingUtils$ggplot(aes(x = mean, y = variance)) +
    #   geom_point()
    
    celMvModel <- lm(variance ~ poly(mean, 3, raw = TRUE), celGamMMs)
    # celMvModel <- loess(variance ~ mean, celGamMMs) # loess
    
    # let's check model fit? 
    # celGamMMs %>%
    #   mutate(fit = celMvModel$fitted.values) %>%
    #   session$graphingUtils$ggplot(aes(x = mean)) +
    #   geom_point(aes(y = variance)) +
    #   geom_line(aes(y = fit), color = "blue", size = 1)
    
    celErrorSd <- celMvModel$residuals %>% sd()
    
    print(paste0("fitting cel gamma models for each of ", length(celltypes), " cell types: ", paste0(celltypes, collapse = ", "), "."))
    
    celGamMMsF <- celltypes %>% lapply(function(currCelltype) {
      
      params <- cteprf[, currCelltype] %>% as_tibble(rownames = "gene_id") %>% dplyr::select(gene_id, mean = value)
      params <- params %>% mutate(mean = log2(mean))
      
      params <- params %>% # generate variances for given genes
        mutate(variance = predict(celMvModel, params)) %>% 
        mutate(variance = variance + rnorm(length(variance), sd = celErrorSd))
      
      # params <- params %>% mutate(variance = pmax(variance, 0.1)) # enforce baseline variance
      # params <- params %>% mutate(mean = pmax(mean, 1))
      
      paramsMat <- params %>% dplyr::select(mean, variance) %>% as.matrix()
      rownames(paramsMat) <- params$gene_id
      paramsMat <- (2^paramsMat)
      
      return(paramsMat)
    })
    
    # checking the simulated version
    # celGamMMsF$microglia %>% as_tibble(rownames = "gene_id") %>% 
    #   session$graphingUtils$ggplot() + 
    #   geom_point(aes(x = log2(mean + 1), y = log2(variance + 1)))
    
    # now convert the mean-var tibble to shape/scale paramtrization for gamma model definitions
    
    celGamParamsF <- celGamMMsF %>% lapply(private$utils$gamMMToParams)
    
    # COMMIT
    this$mdParams$cel <- celGamParamsF
    
    return(this)
  }
  
  api$simBseLnCels <- function(this, nSubject, nCell) {
    # this function simulates cells from the single baseline marginal distribution per gene (same mean)
    # includes the cell level co-expresssion here
    
    celltypes <- this$refdat$cteprf %>% colnames() %>% sort()
    names(celltypes) <- celltypes
    
    subjects <- paste0("sbj_", formatC(1:nSubject, digits = ceiling(log10(nSubject)), flag = "0"))
    names(subjects) <- subjects
    
    celParams <- this$mdParams$cel
    
    print(paste0("Simulating N = ", nCell, " cells for ", "S = ", length(subjects), " subjects per C = ", length(celltypes), " cell types."))
    
    coexPrograms <- this$refdat$coexPrograms$cel
    
    # add additional coexProgram of 0 for any cell type not included
    prgmCelltypes <- coexPrograms %>% lapply(function(program) { program$celltypes }) %>% unlist() %>% unname()
    
    unspecCelltypes <- celltypes %>% setdiff(prgmCelltypes)
    
    if (length(unspecCelltypes) > 0) {
      coexPrograms <- coexPrograms %>% c(list(unspecified = list(celltypes = unspecCelltypes, genes = character(), values = 0)))
    }

    exprmat <- do.call(cbind, coexPrograms %>% lapply(function(coexProgram) {
      
      celltypes <- coexProgram$celltypes
      
      currCormat <- matrix(ncol = length(coexProgram$genes), nrow = length(coexProgram$genes))
      rownames(currCormat) <- coexProgram$genes
      colnames(currCormat) <- coexProgram$genes
      diag(currCormat) <- 1
      currCormat[which(is.na(currCormat))] <- coexProgram$value
      
      exprmat <- do.call(cbind, celltypes %>% lapply(function(currCelltype) {
        
        print(paste0("Simulating cell level expression values for cell type: ", currCelltype))
        
        currParams <- celParams[[currCelltype]]
        
        do.call(cbind, subjects %>% lapply(function(currSbj) { # separate cells by subject, though everything is actually coming from the same distribution (baseline)
          print(paste0("Simulating cell level expression values for N = ", nCell, " cells for cell type: ", currCelltype, " and subject: ", currSbj, "."))
          private$utils$generateVals(nCell, currParams, currCormat, idNames = paste0(currCelltype, ".", currSbj, ".cel"))
        }))
        
      })) 
      
    }))
    
    samples <- tibble(sample = colnames(exprmat)) %>% 
      mutate(subject = str_extract(sample, "sbj_[0-9]+"), 
             cell_type = str_extract(sample, "^[A-Za-z]+")) 
    
    return(list(simCells = list(exprmat = exprmat, samples = samples), 
                params = list(md = celParams, coexPrograms = coexPrograms)))
    
    
    # for debugging
    # n <- nCell
    # mdParams <- currParams 
    # corMat <- currCormat 
    # idNames <- paste0(currCelltype, ".", currSbj, ".cel")
  }
  
  api$simSbjLvMeans <- function(this, nSubject) {
    # generate subject level means based on gamma model
    # this subject level variability can then be combined with cell level models to generate cells with subject level biases
    
    celltypes <- this$refdat$cteprf %>% colnames() %>% sort()
    names(celltypes) <- celltypes
    
    print(paste0("Simulating mean expression values for S = ", nSubject, " subjects per C = ", length(celltypes), " cell types."))
    
    coexPrograms <- this$refdat$coexPrograms$sbj
    
    # add additional coexProgram of 0 for any cell type not included
    prgmCelltypes <- coexPrograms %>% lapply(function(program) { program$celltypes }) %>% unlist() %>% unname()
    
    unspecCelltypes <- celltypes %>% setdiff(prgmCelltypes)
    
    if (length(unspecCelltypes) > 0) {
      coexPrograms <- coexPrograms %>% c(list(unspecified = list(celltypes = unspecCelltypes, genes = character(), values = 0)))
    }
    
    exprmat <- do.call(cbind, coexPrograms %>% lapply(function(coexProgram) {
      
      currCelltypes <- coexProgram$celltypes
      
      print(paste0("Simulating subject level expression values for cell types: ", paste0(currCelltypes, collapse = ", ")))
      
      currCormat <- matrix(ncol = length(coexProgram$genes) * length(currCelltypes), nrow = length(coexProgram$genes) * length(currCelltypes))
      
      if (length(currCormat) > 0) {
        rownames(currCormat) <- currCelltypes %>% lapply(function(celltype) { paste0(celltype, "&", coexProgram$genes) }) %>% unlist()
        colnames(currCormat) <- rownames(currCormat)
        diag(currCormat) <- 1
        currCormat[which(is.na(currCormat))] <- coexProgram$value
      }

      currParams <- currCelltypes %>% lapply(function(celltype) {
        currMdParams <- this$mdParams$sbj[[celltype]]
        rownames(currMdParams) <- paste0(celltype, "&", rownames(currMdParams))
        return(currMdParams)
      }) %>% session$dataWrangler$rbind()
      
      currExprmat <- private$utils$generateVals(nSubject, currParams, currCormat, idNames = paste0("sbj")) 
      
      # now separate into cell types
      
      currExprmat <- do.call(cbind, currCelltypes %>% lapply(function(celltype) {
        exprmatCt <- currExprmat[grepl(celltype, rownames(currExprmat)), ]
        rownames(exprmatCt) <- rownames(exprmatCt) %>% strsplit("&") %>% sapply(function(strs) { strs[2] })
        colnames(exprmatCt) <- paste0(celltype, ".", colnames(exprmatCt))
        return(exprmatCt)
      }))
      
      return(currExprmat)
    }))
    
    samples <- tibble(sample = colnames(exprmat)) %>% 
      mutate(subject = str_extract(sample, "sbj_[0-9]+"), 
             cell_type = str_extract(sample, "^[A-Za-z]+")) 
    
    return(list(simSbjs = list(exprmat = exprmat, samples = samples), 
                params = list(md = this$mdParams$sbj, coexPrograms = coexPrograms)))
    
  }
  
  #excitatory.8739
  
  api$computeCelParams <- function(this, sbjLvMeans) {
    # given the subject level means, make cell level params that define a cell level distribution model per cell type per subject
    
    samples <- sbjLvMeans$samples
    exprmat <- sbjLvMeans$exprmat
    
    celltypes <- samples$cell_type %>% unique() %>% sort()
    names(celltypes) <- celltypes
    
    subjects <- samples$subject %>% unique() %>% sort()
    names(subjects) <- subjects
    
    print(paste0("Constructing final cell level marginal probability distribution models for C = ", 
                 length(celltypes), 
                 " cell types and S = ", 
                 length(subjects), 
                 " subjects."))
    
    # next, construct cell level models using sbjMus as input
    # basically, output an independent model per sbj/cell type combination
    celParams <- celltypes %>% lapply(function(currCelltype) {
      
      # convert gamma into mean var, move the mean, maintain rv, get the new variance, convert back to gamma 
      
      print(paste0("Constructing cell level models for cell type: ", currCelltype))
      
      mdParams <- this$mdParams$cel[[currCelltype]] 
      
      mdMM <- mdParams %>% private$utils$computeGamMM()
      
      rvs <- mdMM %>% apply(1, function(vals) { vals[2] / vals[1] })
      
      currSamples <- samples %>% filter(cell_type == currCelltype)
      currExprmat <- exprmat[, currSamples$sample]
      
      paramsF <- subjects %>% lapply(function(currSbj) {
        
        print(paste0("Constructing cell level models for cell type: ", currCelltype, " and subject: ", currSbj, "."))
        
        currSbjMu <- currExprmat[, grep(currSbj, colnames(currExprmat), value = TRUE)]
        currCelRv <- rvs[names(currSbjMu)]
        
        currVarAdj <- currSbjMu * currCelRv
        
        paramF <- cbind(mean = currSbjMu, variance = currVarAdj) %>% as.matrix()
        paramF <- paramF %>% apply(1, function(vals) { pmax(vals, 1e-10) }) %>% t() 
        
        paramF <- paramF %>% private$utils$gamMMToParams()
        
        return(paramF) #"excitatory.54997"
      })
      
      return(paramsF)
    }) 
    
    return(celParams)
  }
  
  
  api$convertCelLvDist <- function(simCells, celMdParamsOrig, celMdParamsNew) {
    #  convert the cell level expressions into the distributional models with sbj level variability added
    # basically, use the original params to convert the cells into pvalues and then convert those pvalues into corresponding values defined by the new params
    
    exprmat <- simCells$exprmat
    samples <- simCells$samples
    
    celltypes <- samples$cell_type %>% unique() %>% sort()
    names(celltypes) <- celltypes
    
    subjects <- samples$subject %>% unique() %>% sort()
    names(subjects) <- subjects
    
    print(paste0("Converting distributions for N = ", 
                 nrow(samples),
                 " cells across C = ", 
                 length(celltypes), 
                 " cell types and S = ", 
                 length(subjects), 
                 " subjects."))
    
    exprmatF <- do.call(cbind, celltypes %>% lapply(function(currCelltype) {
      
      print(paste0("Converting cells for cell type: ", currCelltype))
      
      result <- do.call(cbind, subjects %>% lapply(function(currSbj) {
        
        print(paste0("Converting cells for cell type: ", currCelltype, " and subject: ", currSbj, "."))
        
        currSamples <- samples %>% filter(cell_type == currCelltype, subject == currSbj)
        currExprmat <- exprmat[, currSamples$sample]
        
        # now, convert everything into p-values based on the original model
        curParamsOrig <- celMdParamsOrig[[currCelltype]]
        
        # then, convert the p-values into the new distribution
        curParamsNew <- celMdParamsNew[[currCelltype]][[currSbj]]
        
        result <- rownames(currExprmat) %>% sapply(function(geneId) {
          orig <- curParamsOrig[geneId,]
          new <- curParamsNew[geneId,]
          pvalues <- currExprmat[geneId, ] %>% pgamma(shape = orig[["shape"]], scale = orig[["scale"]])
          exprNew <- pvalues %>% qgamma(shape = new[["shape"]], scale = new[["scale"]])
          return(exprNew)
        }) %>% t() # transpose into rows are genes just before returning
        
        return(result)
      }))
      
      return(result)
      
    }))
    
    return(list(exprmat = exprmatF, samples = samples))
  }
  
  
  
  api$GENERATE_CC_SPECS <- function(nBkSamples, nTotalCells, baselineProps, sdFrac) { 
    # generate cellular composition specifications
    
    nCells <- baselineProps * nTotalCells
    sds <- sdFrac * nCells
    params <- cbind(nCells, sds)
    
    celltypes <- names(baselineProps)
    names(celltypes) <- celltypes
    
    result <- celltypes %>% sapply(function(currType) {
      currParams <- params[currType,]
      rnorm(nBkSamples, mean = currParams[1], sd = currParams[2]) %>% round()
    })
    
    rownames(result) <- paste0("sbj_", formatC(1:nBkSamples, digits = ceiling(log10(nBkSamples)), flag = "0"))
    
    return(result)
    
  }
  
  
  
  
  
  
  return(list(this = private$utils$initThis(), api = api))
  
}

