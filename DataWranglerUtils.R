constructDataWrangler <- function() {
  
  public <- list()
  
  public$readTable <- function(file, sep, header = TRUE, ...) {
    tibble <- read.table(file, 
                         stringsAsFactors = FALSE, 
                         header = header, 
                         quote = "", 
                         sep = sep, 
                         colClasses = "character", 
                         fill = TRUE,
                         ...) %>% 
      as_tibble()
    return(tibble)
  }
  
  
  public$readTsv <- function(file, ...) {
    tibble <- public$readTable(file, "\t", ...)
    return(tibble)
  }
  
  public$readCsv <- function(file, ...) {
    tibble <- public$readTable(file, ",", ...)
    return(tibble)
  }
  
  public$extractIndices <- function(vector, indices) {
    vector[indices]
  }
  
  # use this function only if the default functions don't work; probably inefficient
  # everything will be read as character here
  public$readCols <- function(file, indices, sep) {
    conn <- file(file, open = "r")
    lines <- readLines(conn)
    close(conn)
    
    header <- lines[1] %>% strsplit(sep) %>% unlist() %>% public$extractIndices(indices)
    
    resultTibble <- do.call(cbind, lines[2:length(lines)] %>% lapply(function(currLine) {
      currRow <- currLine %>% strsplit(sep) %>% unlist() %>% public$extractIndices(indices)
      tibble(currRow)
    })) %>% t() %>% as_tibble()
    
    names(resultTibble) <- header
    
    return(resultTibble)
  }
  
  public$vectorToTibble <- function(vector) {
    
    tibble(variable = names(vector), value = vector)
    
  }
  
  
  public$setColAsRownames <- function(tibble, colName) {
    tibble %>% 
      as.data.frame() %>% 
      remove_rownames() %>% 
      column_to_rownames(colName)
  }
  
  public$setRownameAsColumn <- function(df, colName) {
    df %>% 
      as.data.frame() %>% 
      rownames_to_column(colName) %>% 
      as_tibble()
  }
  
  public$extractColumn <- function(df, colName) {
    df[[colName]]
  }
  
  
  public$topN <- function(vector, n) {
    vector[1:n]
  }
  
  
  
  public$isInvalid <- function(value) {
    return(is.na(value) | tolower(value) == "" | tolower(value) == "na")
  }
  
  public$isNotApplicable <- function(value) {
    return(tolower(value) == "not_applicable")
  }
  
  public$isUnknown <- function(value) {
    return(tolower(value) == "unknown")
  }
  
  public$removeInvalidRows <- function(df, cols = character(0)) {
    if (length(cols) == 0) {
      cols <- names(df)
    }
    if (length(cols) >= 1) {
      if (length(cols) > 1) {
        df[which(!public$isInvalid(df[[cols[1]]])),] %>% public$removeInvalidRows(cols[2:length(cols)])
      } else {
        df[which(!public$isInvalid(df[[cols[1]]])),]
      }
    } else {
      df
    }
  }
  
  public$factorToCharacter <- function(df) {
    dfClone <- df
    for (col in names(dfClone)) {
      if (is.factor(dfClone[[col]])) {
        dfClone[[col]] <- as.character(dfClone[[col]])
      }
    }
    dfClone
  }
  
  public$colToFactor <- function(df, colName) {
    dfClone <- df
    dfClone[[colName]] <- factor(as.character(dfClone[[colName]]))
    dfClone
  }
  
  public$quantileNormalize <- function(df, colNames) {
    
    # https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
    quantile_normalisation <- function(df){
      df_rank <- apply(df,2,rank,ties.method="min")
      df_sorted <- data.frame(apply(df, 2, sort))
      df_mean <- apply(df_sorted, 1, mean)
      
      index_to_mean <- function(my_index, my_mean){
        return(my_mean[my_index])
      }
      
      df_final <- apply(df_rank, 2, index_to_mean, my_mean = df_mean)
      rownames(df_final) <- rownames(df)
      return(df_final)
    }
    
    dfClone <- df
    additionalCols <- names(dfClone)[which(!(names(dfClone) %in% colNames))]
    targetCols <- colNames
    
    normalized <- dfClone[targetCols] %>% as.data.frame() %>% quantile_normalisation() %>% as_tibble()
    
    result <- dfClone[additionalCols] %>% cbind(normalized) %>% as_tibble()
    
    return(result)
  }
    
  
  public$colsToCharacter <- function(targetTibble, colNames) {
    tibbleClone <- targetTibble
    for (currName in colNames) {
      tibbleClone[[currName]] <- as.character(tibbleClone[[currName]])
    }
    return(tibbleClone)
  }
  
  public$colsToNumeric <- function(targetTibble, colNames) {
    tibbleClone <- targetTibble
    for (currName in colNames) {
      tibbleClone[[currName]] <- as.numeric(tibbleClone[[currName]])
    }
    return(tibbleClone)
  }
  
  public$removeInvalidValues <- function(expressionData) {
    dataClone <- expressionData %>% na.omit()
    validRows <- which(is.finite(rowSums(dataClone[2:length(dataClone)])))
    dataClone <- dataClone[validRows,] %>% select(gene, everything())
    return(dataClone)
  }
  
  public$colsToZero <- function(targetTibble, colNames) {
    tibbleClone <- targetTibble
    for (currName in colNames) {
      tibbleClone[[currName]][which(is.na(targetTibble[[currName]]))] <- 0
    }
    return(tibbleClone)
  }
  
  public$zeroToNa <- function(targetTibble, colNames) {
    tibbleClone <- targetTibble
    for (currName in colNames) {
      tibbleClone[[currName]][which((targetTibble[[currName]]) == 0)] <- NA
    }
    return(tibbleClone)
  }
  
  public$fillNa <- function(targetTibble, colNames = NULL, value) {
    tibbleClone <- targetTibble
    if (is.null(colNames)) {
      allColNames <- names(tibbleClone)
    } else {
      allColNames <- colNames %>% intersect(names(tibbleClone))
    }
    for (currName in allColNames) {
      tibbleClone[[currName]][which(is.na(targetTibble[[currName]]))] <- value
    }
    return(tibbleClone)
  }
  
  public$fillEmpty <- function(targetTibble, colNames, value) {
    tibbleClone <- targetTibble
    allColNames <- colNames %>% intersect(names(tibbleClone))
    for (currName in allColNames) {
      tibbleClone[[currName]][which((targetTibble[[currName]] == ""))] <- value
    }
    return(tibbleClone)
  }
  
  public$replaceValue <- function(targetTibble, colNames, oldValue, newValue) {
    tibbleClone <- targetTibble
    allColNames <- colNames %>% intersect(names(tibbleClone))
    for (currName in allColNames) {
      tibbleClone[[currName]][which((targetTibble[[currName]] == oldValue))] <- newValue
    }
    return(tibbleClone)
  }
  
  public$colsToFalse <- function(targetTibble, colNames) {
    tibbleClone <- targetTibble
    for (currName in colNames) {
      tibbleClone[[currName]][which(is.na(targetTibble[[currName]]))] <- FALSE
    }
    return(tibbleClone)
  }
  
  public$rbind <- function(tbls) {
    do.call(rbind, tbls)
  } 
  
  public$cbind <- function(tbls) {
    do.call(cbind, tbls)
  } 
  
  public$mergeColumnsXy <- function(tibble, colName) {
    clone <- tibble
    cols <- names(clone)
    xCol <- paste0(colName, ".x")
    yCol <- paste0(colName, ".y")
    if ((xCol %in% cols) & (yCol %in% cols)) {
      clone[[colName]] <- clone[[xCol]]
      missingIndices <- which(public$isInvalid(clone[[colName]]))
      clone[[colName]][missingIndices] <- clone[[yCol]][missingIndices]
      clone <- clone[!(names(clone) %in% c(xCol, yCol))]
    }
    return(clone)
  }
  
  
  public$attachNames <- function(content, useNames = NULL) {
    if (is.null(useNames)) {
      names(content) <- content
      return(content)
    } else {
      names(content) <- useNames
      return(content)
    }
  }
  
  
  public$extractNestedField <- function(strs, field) {
    strs %>% str_extract(paste0("(?<=", field, "=)\\*"))
  }
  
  public$nestFields <- function(tibble, cols, new) {
    tibbleClone <- tibble
    for (col in cols) {
      tibbleClone[[col]] <- paste0(col, "=", tibbleClone[[col]])
    }
    tibbleClone[[new]] <- paste0("[", apply(tibbleClone[,cols], 1, paste0, collapse = ","), "]")
    tibbleClone <- tibbleClone[,!(names(tibbleClone) %in% cols)] 
    return(tibbleClone)
  }
  
  return(public)
}
