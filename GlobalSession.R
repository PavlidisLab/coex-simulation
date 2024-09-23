globalSession <- list()
globalSession$static$HOME_DIR <- "/home/echu113/ws_rstudio/"

initGlobalSession <- function() {
  
  # LOAD INTERNAL DEPENDENCIES ===================================================
  library(tidyverse)
  library(limma)
  library(pheatmap)
  library(homologene)
  library(parallel)
  # ==============================================================================
  
  # -----------------------------
  # PUBLIC
  # -----------------------------
  
  public <- list()
  
  public$HOME_DIR <- globalSession$static$HOME_DIR
  public$WORKSPACES <- paste0(public$HOME_DIR, "_workspaces/")
  public$UTILS_DIR <- paste0(public$HOME_DIR, "utils/")
  public$EXPRESSION_DATA_DIR <- paste0(public$HOME_DIR, "expression_data/")
  public$BEDS_DATA_DIR <- paste0(public$HOME_DIR, "beds_data/")
  
  
  # LOAD INTERNAL DEPENDENCIES ===================================================
  source(paste0(public$UTILS_DIR, "CollectionUtils.R")); public$collectionUtils <- constructCollectionUtils()
  source(paste0(public$UTILS_DIR, "GraphingUtils.R")); public$graphingUtils <- constructGraphingUtils()
  source(paste0(public$UTILS_DIR, "DataWranglerUtils.R")); public$dataWrangler <- constructDataWrangler()
  source(paste0(public$UTILS_DIR, "EvaluationUtils.R")); public$evaluationUtils <- constructEvaluationUtils()
  # source(paste0(public$UTILS_DIR, "OntologyUtils.R")); public$ontologyUtils <- constructOntology(public)
  # source(paste0(public$UTILS_DIR, "LiftOverUtils.R")); public$liftOverUtils <- constructLiftOverUtils(public)
  # ==============================================================================
  
  return(public)
}

