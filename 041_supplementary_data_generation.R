# in this script generate all the supplementary data that need to be uploaded

source("/home/echu113/ws_rstudio/_workspaces/041_CH4_FINAL/041_SETUP.R")

# uncomment to run each numbered secion; this is so I can isolate each section to run in the terminal

# # ============================
# # 1. single cell co-expression matrices
# # - 39? one per dataset and cell type
# 
# fdata <- list()
# 
# fdata$files <- tibble(outfile = list.files(workspace$outputDir)) %>% 
#   filter(grepl("coexmats", outfile)) %>% 
#   mutate(outfile_name = gsub(".rds", "", outfile)) %>% 
#   mutate(strs = strsplit(outfile_name, "_")) %>% 
#   mutate(level = strs %>% sapply(function(currStr) { currStr[1] })) %>% 
#   mutate(type = strs %>% sapply(function(currStr) { currStr[2] })) %>% 
#   mutate(dataset = strs %>% sapply(function(currStr) { currStr[3] })) %>% 
#   mutate(cell_type = strs %>% sapply(function(currStr) { currStr[4] })) %>% 
#   dplyr::select(outfile, level, type, dataset, cell_type) %>% 
#   filter(level == "sc", dataset != "consensus")
# 
# # ramos doesn't have astrocyte and oligodendrocyte; lau doesn't have opc
# # that's why its missing 3, from 6 cell type * 7 datasets = 42 data groups, only 39 available
# fdata$files %>% session$graphingUtils$ggplot(aes(x = dataset, y = cell_type)) + geom_point()
# 
# # =================================================
# # OUTPUT S_FILE
# xdata <- list()
# 
# print("Writing xcell co-expression files...")
# 
# xdata$files <- fdata$files$outfile %>% session$dataWrangler$attachNames()
# 
# xdata$files %>% mclapply(function(currFile) {
#   mat <- read_rds(paste0(workspace$outputDir, currFile)) %>% session$dataWrangler$setRownameAsColumn("gene_id")
#   metadata <- fdata$files %>% filter(outfile == currFile)
#   filename <- paste0("xcell_coexmat_", metadata$cell_type, "_", metadata$dataset, ".tsv")
#   write.table(mat, file = paste0(workspace$sdataDir, filename), sep = "\t", row.names = FALSE)
# }, mc.cores = 40)
# # =================================================




# # ============================
# # 2. xSubject co-expression matrices
# # - 39? one per dataset and cell type
# 
# fdata <- list()
# 
# fdata$files <- tibble(outfile = list.files(workspace$outputDir)) %>% 
#   filter(grepl("coexmats", outfile)) %>% 
#   mutate(outfile_name = gsub(".rds", "", outfile)) %>% 
#   mutate(strs = strsplit(outfile_name, "_")) %>% 
#   mutate(level = strs %>% sapply(function(currStr) { currStr[1] })) %>% 
#   mutate(type = strs %>% sapply(function(currStr) { currStr[2] })) %>% 
#   mutate(dataset = strs %>% sapply(function(currStr) { currStr[3] })) %>% 
#   mutate(cell_type = strs %>% sapply(function(currStr) { currStr[4] })) %>% 
#   dplyr::select(outfile, level, type, dataset, cell_type) %>% 
#   filter(level == "sbj", dataset != "consensus")
# 
# # ramos doesn't have astrocyte and oligodendrocyte; lau doesn't have opc
# # that's why its missing 3, from 6 cell type * 7 datasets = 42 data groups, only 39 available
# fdata$files %>% session$graphingUtils$ggplot(aes(x = dataset, y = cell_type)) + geom_point()
# 
# # =================================================
# # OUTPUT S_FILE
# xdata <- list()
# 
# print("Writing xsbj co-expression files...")
# 
# xdata$files <- fdata$files$outfile %>% session$dataWrangler$attachNames()
# 
# xdata$files %>% mclapply(function(currFile) {
#   mat <- read_rds(paste0(workspace$outputDir, currFile)) %>% session$dataWrangler$setRownameAsColumn("gene_id")
#   metadata <- fdata$files %>% filter(outfile == currFile)
#   filename <- paste0("xsbj_coexmat_", metadata$cell_type, "_", metadata$dataset, ".tsv")
#   write.table(mat, file = paste0(workspace$sdataDir, filename), sep = "\t", row.names = FALSE)
# }, mc.cores = 40)
# # =================================================

# # ============================
# # 3. consensus xCell networks
# # - 6, one per cell type
# 
# fdata <- list()
# 
# fdata$file <- "sc_coexmats_consensus.rds"
# 
# print("Writing xcell consensus network...")
# 
# fdata$mat <- read_rds(paste0(workspace$outputDir, fdata$file)) 
# 
# names(fdata$mat) %>% mclapply(function(celltype) {
#   
#   coexmat <- fdata$mat[[celltype]] %>% session$dataWrangler$setRownameAsColumn("gene_id")
#   
#   filename <- paste0("xcell_coexmat_", celltype, "_consensus.tsv")
#   write.table(coexmat, file = paste0(workspace$sdataDir, filename), sep = "\t", row.names = FALSE)
#   
# }, mc.cores = 6)
# # =================================================


# # ==================================================
# # 4. bulk, before and after correction --- coex values --- which has been normalized
# # - 3, rosmap, rosmap_ihc, velmeshev
# 
# fdata <- list()
# 
# fdata$files <- c(rosmap = "bk_coexmats_rosmap.rds", rosmap_mgp_residual = "bk_coexmats_rosmap-mgpres.rds",
#                  rosmap_ihc = "bk_coexmats_rosmap-ihc.rds", rosmap_ihc_residual = "bk_coexmats_rosmap-ihcres.rds",
#                  velmeshev = "bk_coexmats_velmeshev.rds", velmeshev_mgp_residual = "bk_coexmats_velmeshev-mgpres.rds")
# 
# print("Writing bulk co-expression coex files...")
# 
# names(fdata$files) %>% mclapply(function(file) {
# 
#   coexmat <- read_rds(paste0(workspace$outputDir, fdata$files[[file]])) %>% session$dataWrangler$setRownameAsColumn("gene_id")
# 
#   filename <- paste0("bk_coexmat_", file, ".tsv")
#   write.table(coexmat, file = paste0(workspace$sdataDir, filename), sep = "\t", row.names = FALSE)
# 
# }, mc.cores = 6)
# 
# # =================================================

# # ==================================================
# # 5. bulk, before and after correction --- coefficient -- not normalized, just pearson's correlation coefficients
# # - 3, rosmap, rosmap_ihc, velmeshev
# 
# fdata <- list()
# 
# fdata$files <- c(rosmap = "bk_coefmats_rosmap.rds", rosmap_mgp_residual = "bk_coefmats_rosmap-mgpres.rds",
#                  rosmap_ihc = "bk_coefmats_rosmap-ihc.rds", rosmap_ihc_residual = "bk_coefmats_rosmap-ihcres.rds",
#                  velmeshev = "bk_coefmats_velmeshev.rds", velmeshev_mgp_residual = "bk_coefmats_velmeshev-mgpres.rds")
# 
# print("Writing bulk co-expression coef files...")
# 
# names(fdata$files) %>% mclapply(function(file) {
#   
#   coexmat <- read_rds(paste0(workspace$outputDir, fdata$files[[file]])) %>% session$dataWrangler$setRownameAsColumn("gene_id")
#   
#   filename <- paste0("bk_coefmat_", file, ".tsv")
#   write.table(coexmat, file = paste0(workspace$sdataDir, filename), sep = "\t", row.names = FALSE)
#   
# }, mc.cores = 6)
# 
# # =================================================


# # ==================================================
# # 6. CCV models for each residual? 
# # - 3, rosmap, rosmap_ihc, velmeshev
# 
# fdata <- list()
# 
# fdata$files <- list(rosmap_mgp = "bk_exprmats_rosmap-mgpres.rds", 
#                  velmeshev_mgp = "bk_exprmats_velmeshev-mgpres.rds",
#                  rosmap_ihc = "bk_exprmats_rosmap-ihcres.rds")
# 
# fdata$main <- fdata$files %>% lapply(function(currfile) {
#   read_rds(paste0(workspace$outputDir, currfile))
# })
# 
# # 1. ccv markers
# 
# xdata <- list()
# 
# xdata$main <- fdata$main$rosmap_mgp$mgps$mkrs$tbl %>% mutate(dataset = "rosmap") %>% 
#   rbind(fdata$main$velmeshev_mgp$mgps$mkrs$tbl %>% mutate(dataset = "velmeshev"))
# 
# write.table(xdata$main, file = paste0(workspace$sdataDir, "bk_ccv_mkrs.tsv"), sep = "\t", row.names = FALSE)
#  
# # 2. ccv cell type proportions
# 
# xdata <- list()
# 
# xdata$main <- fdata$main$rosmap_mgp$ccvModels$ctpMat %>% 
#   session$dataWrangler$setRownameAsColumn("sample") %>% 
#   gather(cell_type, estimate, -sample) %>% 
#   mutate(dataset = "rosmap", method = "mgp") %>% 
#   rbind(fdata$main$rosmap_ihc$ccvModel$ctpMat %>% 
#           session$dataWrangler$setRownameAsColumn("sample") %>%
#           gather(cell_type, estimate, -sample) %>% 
#           mutate(dataset = "rosmap", method = "ihc")) %>% 
#   rbind(fdata$main$velmeshev_mgp$ccvModels$ctpMat %>% 
#           session$dataWrangler$setRownameAsColumn("sample") %>% 
#           gather(cell_type, estimate, -sample) %>% 
#           mutate(dataset = "velmeshev", method = "mgp")) 
# 
# write.table(xdata$main, file = paste0(workspace$sdataDir, "bk_ccv_estimates.tsv"), sep = "\t", row.names = FALSE)
# 
# # 3. ccv statistics lm
# 
# xdata <- list()
# 
# xdata$main <- fdata$main$rosmap_mgp$ccvModels$stats$lmStats %>% mutate(dataset = "rosmap", method = "mgp") %>% 
#   rbind(fdata$main$rosmap_ihc$ccvModel$stats$lmStats %>% mutate(dataset = "rosmap", method = "ihc")) %>% 
#   rbind(fdata$main$rosmap_ihc$ccvModel$stats$lmStats %>% mutate(dataset = "velmeshev", method = "mgp"))
# 
# write.table(xdata$main, file = paste0(workspace$sdataDir, "bk_ccv_stat_lm.tsv"), sep = "\t", row.names = FALSE)
# 
# 
# # 4. ccv statistics lm
# 
# xdata <- list()
# 
# xdata$main <- fdata$main$rosmap_mgp$ccvModels$stats$coefStats %>% mutate(dataset = "rosmap", method = "mgp") %>% 
#   rbind(fdata$main$rosmap_ihc$ccvModel$stats$coefStats %>% mutate(dataset = "rosmap", method = "ihc")) %>% 
#   rbind(fdata$main$rosmap_ihc$ccvModel$stats$coefStats %>% mutate(dataset = "velmeshev", method = "mgp"))
# 
# write.table(xdata$main, file = paste0(workspace$sdataDir, "bk_ccv_stat_coef.tsv"), sep = "\t", row.names = FALSE)
# 
# ==========================================================


# # ==========================================================
# # 7. odds_ratio_repro (with a column for levels, ref_dataset, test_dataset)
# 
# fdata <- list(); gc()
# 
# fdata$main <- c(xcell_reproducibility =  "sc_repro_fisher.rds",
#                 xsbj_reproducibility = "sbj_repro_fisher.rds", 
#                 xbulk_reproducibility = "bk_repro_fisher.rds")
# 
# fdata$main <- fdata$main %>% 
#   session$collectionUtils$lapply(function(file) { read_rds(paste0(workspace$outputDir, file)) })
# 
# fdata$main <- fdata$main %>% session$collectionUtils$lapply(function(fisher) { fisher$tbl })
# 
# fdata$main <- fdata$main$xcell_reproducibility %>%
#   mutate(level = "xcell") %>% 
#   dplyr::select(level, 
#                 cell_type, 
#                 ref_dat = dat_truth, 
#                 test_dat = dat_sample, 
#                 odds_ratio = or, 
#                 n_pairs_recov = recovered_n, 
#                 frac_pairs_recov = recovered_frac, 
#                 pvalue) %>% 
#   rbind(fdata$main$xsbj_reproducibility %>%
#           mutate(level = "xsubject") %>% 
#           dplyr::select(level, 
#                         cell_type, 
#                         ref_dat = dat_truth, 
#                         test_dat = dat_sample, 
#                         odds_ratio = or, 
#                         n_pairs_recov = recovered_n, 
#                         frac_pairs_recov = recovered_frac, 
#                         pvalue)) %>% 
#   rbind(fdata$main$xbulk_reproducibility %>%
#           mutate(level = "xbulk", cell_type = NA) %>% 
#           dplyr::select(level, 
#                         cell_type, 
#                         ref_dat = dat_truth, 
#                         test_dat = dataset, 
#                         odds_ratio = or, 
#                         n_pairs_recov = recovered_n, 
#                         frac_pairs_recov = recovered_frac, 
#                         pvalue))
# 
# write.table(fdata$main, file = paste0(workspace$sdataDir, "odds_ratio_repro.tsv"), sep = "\t", row.names = FALSE)
# # ==========================================================


# # ==========================================================
# # 8. odds_ratio_preserv (dataset, ref_lvl, test_lvl)
# 
# fdata <- list(); gc()
# 
# fdata$main <- c(xsbj_preservation = "sc_sbj_preservation_fisher.rds", 
#                 xbulk_preservation = "bk_preservation_fisher.rds")
# 
# fdata$main <- fdata$main %>% 
#   session$collectionUtils$lapply(function(file) { read_rds(paste0(workspace$outputDir, file)) })
# 
# fdata$main <- fdata$main %>% session$collectionUtils$lapply(function(fisher) { fisher$tbl })
# 
# fdata$main <- fdata$main$xsbj_preservation %>%
#   mutate(ref_lvl = "xcell", test_lvl = "xsubject") %>% 
#   dplyr::select(ref_lvl, 
#                 test_lvl, 
#                 cell_type, 
#                 dataset = dat_truth, 
#                 odds_ratio = or, 
#                 n_pairs_recov = recovered_n, 
#                 frac_pairs_recov = recovered_frac, 
#                 pvalue) %>% 
#   rbind(fdata$main$xbulk_preservation %>%
#           mutate(ref_lvl = level, test_lvl = "xbulk") %>% 
#           dplyr::select(ref_lvl, 
#                         test_lvl, 
#                         cell_type, 
#                         dataset, 
#                         odds_ratio = or, 
#                         n_pairs_recov = recovered_n, 
#                         frac_pairs_recov = recovered_frac, 
#                         pvalue)) 
# 
# fdata$main$ref_lvl %>% unique()
# 
# fdata$main <- fdata$main %>% 
#   mutate(ref_lvl = ref_lvl %>% 
#            sapply(function(currlvl) { if (currlvl == "sc") { "xcell" } else if (currlvl == "sbj") { "xsubject" } else { currlvl }}))
# 
# write.table(fdata$main, file = paste0(workspace$sdataDir, "odds_ratio_preserv.tsv"), sep = "\t", row.names = FALSE)
#
# # ==========================================================


# ==========================================================
# 9. clusters information

# 1. cluster genes

fdata <- list(); gc()

fdata$main <- read_rds(paste0(workspace$outputDir, "sc_clusters_objects.rds"))

fdata$main <- fdata$main$clusterFlats %>% unname() %>% unlist(recursive = FALSE) %>% 
  session$collectionUtils$lapplyWithName(function(cluster, genes) { tibble(cluster_id = cluster, gene_id = genes) }) %>% 
  session$dataWrangler$rbind()

write.table(fdata$main, file = paste0(workspace$sdataDir, "clusters_genes.tsv"), sep = "\t", row.names = FALSE)

# 2. cluster integrity at different levels

fdata <- list(); gc()

fdata$main <- read_rds(paste0(workspace$outputDir, "sc_clusters_integrity_within_celltype.rds"))$egad

fdata$main <- fdata$main %>% filter(!grepl("psdbk", dataset), !grepl("generic", cluster)) 

fdata$main <- fdata$main %>% dplyr::select(cluster, cell_type = coex_ct, everything())

fdata$main <- fdata$main %>% dplyr::select(-cluster_ct)

fdata$main <- fdata$main %>% 
  mutate(level = level %>% 
           sapply(function(currLvl) { if (currLvl == "sc") { "xcell" } else if (currLvl == "sbj") { "xsubject" } else if (currLvl == "bk") { "xbulk" } }))

write.table(fdata$main, file = paste0(workspace$sdataDir, "clusters_integrity.tsv"), sep = "\t", row.names = FALSE)


# 3. cluster average co-expression in bulk

fdata <- list(); gc()

fdata$main <- read_rds(paste0(workspace$outputDir, "sc_clusters_integrity_in_bk_coefmats.rds"))

fdata$main <- fdata$main %>% filter(!grepl("psdbk", dataset), !grepl("generic", cluster)) 

fdata$main <- fdata$main %>% dplyr::select(cluster, cell_type = coex_ct, everything())

fdata$main <- fdata$main %>% dplyr::select(-cluster_ct)







