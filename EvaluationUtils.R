constructEvaluationUtils <- function() {
  
  public <- list()
  
  private <- list()
  
  private$computeRecall <- function(p, n, t) {
    tp <- p %>% intersect(t) %>% length()
    fn <- t %>% intersect(n) %>% length()
    return(tp / (tp + fn))
  }
  
  private$computeFPR <- function(p, n, t) {
    fp <- p %>% setdiff(t) %>% length()
    tn <- n %>% setdiff(t) %>% length()
    return(fp / (fp + tn))
  }
  
  private$computePrecision <- function(p, n, t) {
    tp <- p %>% intersect(t) %>% length()
    fp <- p %>% setdiff(t) %>% length() 
    return(tp /(tp + fp))
  }
  
  public$roc <- function(ranking, trueSet, step) {
    t <- trueSet %>% intersect(ranking)
    roc <- do.call(rbind, seq(1, length(ranking), step) %>% session$collectionUtils$lapply(function(i) {
      p <- ranking[1:i]
      n <- ranking %>% setdiff(p)
      tibble(recall = private$computeRecall(p, n, t), 
             false_positive_rate = private$computeFPR(p, n, t))
    }, reporter = "ROC"))
    return(roc)
  }
  
  public$precisionRecall <- function(ranking, trueSet, step) {
    t <- trueSet %>% intersect(ranking)
    precisionRecall <- do.call(rbind, seq(1, length(ranking), step) %>% session$collectionUtils$lapply(function(i) {
      p <- ranking[1:i]
      n <- ranking %>% setdiff(p)
      tibble(recall = private$computeRecall(p, n, t), 
             precision = private$computePrecision(p, n, t))
    }, reporter = "PR"))
    return(precisionRecall)
  }
  
  public$plotAuc <- function(x, y) {
    data <- tibble(x = x, y = y)
    data <- data %>% 
      mutate(x_diff = diff(c(0, x))) %>% 
      mutate(unit_area = x_diff * y) 
    data$unit_area %>% sum()
  }
  
  public$auroc <- function(ranking, trueSet) {
    setStatus <- tibble(item = ranking) %>% mutate(rank = 1:length(ranking), in_set = item %in% trueSet)
    result <- wilcox.test(rank ~ in_set, data = setStatus) 
    ns <- setStatus %>% group_by(in_set) %>% summarize(n = n()) %>% arrange(in_set)
    # auc
    auc <- result$statistic / (as.double(ns$n[1]) * as.double(ns$n[2]))
    return(c(n_true = length(trueSet), n_full = length(ranking), auroc = auc, pvalue = result$p.value))
  }
  
  public$aurocCI <- function(ranking, trueSet, iteration) {
    # f <- ranking %>% setdiff(trueSet)
    t <- trueSet
    sampleAuc <- 1:iteration %>% mclapply(function(i) {
      # currF <- f %>% sample(length(f), replace = TRUE) %>% unique()
      currT <- t %>% sample(length(t), replace = TRUE) %>% unique()
      # currRanking <- ranking[which(ranking %in% c(currF, currT))]
      public$auroc(ranking, currT)
    }, mc.cores = 20) %>% unlist()
    (absCI <- (sd(sampleAuc)*1.96))
    mean <- public$auroc(ranking, trueSet)
    return(tibble(mean = mean, lower = mean - absCI, upper = mean + absCI))
  }
  
  
  
  public$precision <- function(ranking, trueSet, recall) {
    
    recallN <- round((length(trueSet) * recall), 0)
    
    trueRanks <- tibble(item = ranking, rank = 1:length(ranking)) %>% 
      mutate(in_set = item %in% trueSet) %>% 
      filter(in_set)
    
    recallRank <- trueRanks[recallN,]$rank
    
    p <- ranking[1:recallRank]
    n <- ranking %>% setdiff(p)
    t <- trueSet
    
    private$computePrecision(p, n, t)
  }
  
  
  public$precisionCI <- function(ranking, trueSet, recall, iteration) {
    f <- ranking %>% setdiff(trueSet)
    t <- trueSet
    samplePrecision <- 1:iteration %>% sapply(function(i) {
      currF <- f %>% sample(length(f), replace = TRUE) %>% unique()
      currT <- t %>% sample(length(t), replace = TRUE) %>% unique()
      currRanking <- ranking[which(ranking %in% c(currF, currT))]
      public$precision(currRanking, currT, recall)
    })
    absCI <- (sd(samplePrecision)*1.96)
    mean <- public$precision(ranking, trueSet, recall)
    return(tibble(lower = mean - absCI, mean = mean, upper = mean + absCI))
  }
  
  
  public$computeOra <- function(hits, trueSet, background) {
    
    allPositives <- trueSet %>% intersect(background) %>% na.omit()
    allHits <- hits %>% intersect(background) %>% na.omit()
    
    truePositives <- allHits %>% intersect(allPositives)
    
    hits_frac = length(hits) / length(background)
    tp_frac = length(truePositives) / length(allPositives)
    
    pvalue = phyper(length(truePositives), length(allPositives), length(background) - length(allPositives), length(hits), lower.tail = FALSE)
    
    tibble(n_background = length(background), 
           n_sample = length(hits), 
           n_true = length(allPositives), 
           n_tp = length(truePositives), 
           sample_frac = hits_frac, 
           tp_frac = tp_frac,
           p_value = pvalue)
  }
  
  
  public$fisher <- function(predicted, notPredicted, trueSet, alternative) {
    # return 1. a fisher's test object, 2. contingency table
    
    predicted_false <- notPredicted
    predicted_true <- predicted %>% setdiff(predicted_false) # make sure the two are mutually exclusive
    background <- c(predicted_true, predicted_false)
    truth <- trueSet %>% intersect(background)
    
    a <- predicted_true %>% intersect(truth) 
    b <- predicted_true %>% setdiff(a)
    
    c <- predicted_false %>% intersect(truth)
    d <- predicted_false %>% setdiff(c)

    contingency <- data.frame(
      true = c(length(a), length(c)), 
      false = c(length(b), length(d)), 
      row.names = c("sampled", "not_sampled")
    )
    
    fisher <- fisher.test(contingency, alternative = alternative)
    
    # also return some basic stats
    # 1. or - odds ratio, 2. enrich, 3. depriv, 4. expected, ... etc
    
    contingency_null <- chisq.test(contingency)$expected
    
    stats <- c()
    stats["or"] <- (contingency[1, 1] / contingency[1, 2]) / (contingency[2, 1] / contingency[2, 2])
    stats["enrich"] <- (contingency[1, 1] / contingency[1, 2])
    stats["depriv"] <- (contingency[2, 1] / contingency[2, 2])
    stats["times_over_random"] <- (contingency[1, 1] / contingency_null[1, 1])
    stats["recovered_n"] <- contingency[1, 1]
    stats["recovered_frac"] <- contingency[1, 1] / length(trueSet) 
    
    return(list(tbl = contingency, tbl_null = contingency_null, test = fisher, stats = stats))
  }
  
  
  
  return(public)
}












  