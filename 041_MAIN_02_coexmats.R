# in this script, compute the co-expression networks
# cross cell - one per cell type - elimiate cell types & samples where appropriate
# cross subjects - one per cell type - eliminate cell types & samples where appropriate 
# cross bulk samples - one 
# use all availage genes

# LOAD IN ALL THE PIECES OF DATA NEEDED
pdata <- list(); gc()

pdata$genes <- readRDS(paste0(workspace$outputDir, "genes_metadata.rds")) # genes

pdata$celltypes <- c("excitatory", "inhibitory", "opc", "oligodendrocyte", "astrocyte", "microglia")
names(pdata$celltypes) <- pdata$celltypes

pdata$files$exprmats <- tibble(outfile = list.files(workspace$outputDir)) %>% 
  filter(grepl("exprmats", outfile)) %>% 
  mutate(outfile_name = gsub(".rds", "", outfile)) %>% 
  filter(!grepl("RAW", outfile_name)) %>% 
  mutate(strs = strsplit(outfile_name, "_")) %>% 
  mutate(level = strs %>% sapply(function(currStr) { currStr[1] })) %>% 
  mutate(dataset = strs %>% sapply(function(currStr) { currStr[3] })) %>% 
  dplyr::select(outfile, level, dataset)

pdata$files$exprmats$dataset %>% unique() %>% sort()

pdata$files$coexmats <- tibble(outfile = list.files(workspace$outputDir)) %>% 
  filter(grepl("coexmats", outfile)) %>% 
  mutate(outfile_name = gsub(".rds", "", outfile)) %>% 
  mutate(strs = strsplit(outfile_name, "_")) %>% 
  mutate(level = strs %>% sapply(function(currStr) { currStr[1] })) %>% 
  mutate(type = strs %>% sapply(function(currStr) { currStr[2] })) %>% 
  mutate(dataset = strs %>% sapply(function(currStr) { currStr[3] })) %>% 
  mutate(cell_type = strs %>% sapply(function(currStr) { currStr[4] })) %>% 
  dplyr::select(outfile, level, type, dataset, cell_type)
  
# +++++++++++++++++++++
# compute cross cell co-expression tables
# velmeshev

fdata <- list(); gc()

fdata$exprdat <- readRDS(paste0(workspace$outputDir, "sc_exprmats_velmeshev.rds"))

fdata$exprmats <- fdata$exprdat$exprmats

fdata$samples <- fdata$exprdat$samples %>% lapply(function(currSamples) {
  currSamples %>% dplyr::select(sample, sample_bk, cell_type) %>% group_by(cell_type, sample_bk) %>% summarize(n_cell = n()) %>% ungroup()
}) %>% session$dataWrangler$rbind()

fdata$samples <- fdata$samples %>% filter(cell_type %in% pdata$celltypes)

fdata$samples %>% 
  session$graphingUtils$ggplot(aes(x = sample_bk, y = n_cell)) + 
  geom_bar(stat = "identity") + 
  geom_hline(yintercept = 20) +
  scale_y_log10() + 
  session$graphingUtils$tiltX(angle = 90) + 
  facet_wrap(~cell_type, ncol = 1) + 
  geom_bar(data = fdata$samples %>% filter(n_cell < 20), fill = "red", stat = "identity")

fdata$samples %>% 
  filter(n_cell >= 20) %>% 
  group_by(cell_type) %>% 
  summarize(n_sbj = n())

fdata$samples <- fdata$samples %>% filter(n_cell >= 20) # filter out any sample with less than 20 cells

# also filter out any cell type with < 20 samples
fdata$samples <- fdata$samples %>% filter(cell_type %in% (fdata$samples %>%
                           group_by(cell_type) %>% 
                           summarize(n_sample_bk = n()) %>% 
                           filter(n_sample_bk >= 20) %>% 
                           session$dataWrangler$extractColumn("cell_type")))

fdata$celltypes <- pdata$celltypes %>% intersect(fdata$samples$cell_type)

# process each cell type and write to disk at each step so I don't hold a lot of memory
fdata$celltypes %>% mclapply(function(currCelltype) { 
  print(paste0("Working on cell type: ", currCelltype, "..."))
  currSamples <- fdata$samples %>% filter(cell_type == currCelltype) 
  currCells <- fdata$exprdat$samples[[currCelltype]] %>% dplyr::select(sample, sample_bk) %>% unique()
  currCells <- currCells %>% filter(sample_bk %in% currSamples$sample_bk)
  currExprmat <- fdata$exprmats[[currCelltype]][, currCells$sample]
  currExprmat <- currExprmat[sort(rownames(currExprmat)), ]
  currCoexmat <- workspace$utils$computeAggScCoex(currCells, currExprmat)
  
  # ==========================================
  # COMMIT ===================================
  print("Saving coexmat to disk...")
  currCoexmat %>% saveRDS(paste0(workspace$outputDir, paste0("sc_coexmats_velmeshev_", currCelltype, ".rds")))
}, mc.cores = length(fdata$celltypes))


# +++++++++++++++++++++
# compute cross cell co-expression tables
# rosmap

fdata <- list(); gc()

fdata$exprdat <- readRDS(paste0(workspace$outputDir, "sc_exprmats_rosmap.rds"))

fdata$exprmats <- fdata$exprdat$exprmats

fdata$samples <- fdata$exprdat$samples %>% lapply(function(currSamples) {
  currSamples %>% dplyr::select(sample, sample_bk, cell_type) %>% group_by(cell_type, sample_bk) %>% summarize(n_cell = n()) %>% ungroup()
}) %>% session$dataWrangler$rbind()

fdata$samples <- fdata$samples %>% filter(cell_type %in% pdata$celltypes)

fdata$samples %>% 
  session$graphingUtils$ggplot(aes(x = sample_bk, y = n_cell)) + 
  geom_bar(stat = "identity") + 
  geom_hline(yintercept = 20) +
  scale_y_log10() + 
  session$graphingUtils$tiltX(angle = 90) + 
  facet_wrap(~cell_type, ncol = 1) + 
  geom_bar(data = fdata$samples %>% filter(n_cell < 20), fill = "red", stat = "identity")

fdata$samples %>% 
  filter(n_cell >= 20) %>% 
  group_by(cell_type) %>% 
  summarize(n_sbj = n())

fdata$samples <- fdata$samples %>% filter(n_cell >= 20) # filter out any sample with less than 20 cells

# also filter out any cell type with < 20 samples
fdata$samples <- fdata$samples %>% filter(cell_type %in% (fdata$samples %>%
                                                            group_by(cell_type) %>% 
                                                            summarize(n_sample_bk = n()) %>% 
                                                            filter(n_sample_bk >= 20) %>% 
                                                            session$dataWrangler$extractColumn("cell_type")))

fdata$celltypes <- pdata$celltypes %>% intersect(fdata$samples$cell_type)

# process each cell type and write to disk at each step so I don't hold a lot of memory
fdata$celltypes %>% mclapply(function(currCelltype) { 
  print(paste0("Working on cell type: ", currCelltype, "..."))
  currSamples <- fdata$samples %>% filter(cell_type == currCelltype) 
  currCells <- fdata$exprdat$samples[[currCelltype]] %>% dplyr::select(sample, sample_bk) %>% unique()
  currCells <- currCells %>% filter(sample_bk %in% currSamples$sample_bk)
  currExprmat <- fdata$exprmats[[currCelltype]][, currCells$sample]
  currExprmat <- currExprmat[sort(rownames(currExprmat)), ]
  currCoexmat <- workspace$utils$computeAggScCoex(currCells, currExprmat)
  
  # ==========================================
  # COMMIT ===================================
  print("Saving coexmat to disk...")
  currCoexmat %>% saveRDS(paste0(workspace$outputDir, paste0("sc_coexmats_rosmap_", currCelltype, ".rds")))
}, mc.cores = length(fdata$celltypes))


# +++++++++++++++++++++
# compute cross cell co-expression tables
# nagy

fdata <- list(); gc()

fdata$dataset <- "nagy"

fdata$exprdat <- readRDS(paste0(workspace$outputDir, "sc_exprmats_", fdata$dataset, ".rds"))

fdata$exprmats <- fdata$exprdat$exprmats

fdata$samples <- fdata$exprdat$samples %>% lapply(function(currSamples) {
  currSamples %>% dplyr::select(sample, sample_bk, cell_type) %>% group_by(cell_type, sample_bk) %>% summarize(n_cell = n()) %>% ungroup()
}) %>% session$dataWrangler$rbind()

fdata$samples <- fdata$samples %>% filter(cell_type %in% pdata$celltypes)

fdata$samples %>% 
  session$graphingUtils$ggplot(aes(x = sample_bk, y = n_cell)) + 
  geom_bar(stat = "identity") + 
  geom_hline(yintercept = 20) +
  scale_y_log10() + 
  session$graphingUtils$tiltX(angle = 90) + 
  facet_wrap(~cell_type, ncol = 1) + 
  geom_bar(data = fdata$samples %>% filter(n_cell < 20), fill = "red", stat = "identity")

fdata$samples %>% 
  filter(n_cell >= 20) %>% 
  group_by(cell_type) %>% 
  summarize(n_sbj = n())

fdata$samples <- fdata$samples %>% filter(n_cell >= 20) # filter out any sample with less than 20 cells

# also filter out any cell type with < 20 samples
fdata$samples <- fdata$samples %>% filter(cell_type %in% (fdata$samples %>%
                                                            group_by(cell_type) %>% 
                                                            summarize(n_sample_bk = n()) %>% 
                                                            filter(n_sample_bk >= 20) %>% 
                                                            session$dataWrangler$extractColumn("cell_type")))

fdata$celltypes <- pdata$celltypes %>% intersect(fdata$samples$cell_type)

# process each cell type and write to disk at each step so I don't hold a lot of memory
fdata$celltypes %>% mclapply(function(currCelltype) { 
  print(paste0("Working on cell type: ", currCelltype, "..."))
  currSamples <- fdata$samples %>% filter(cell_type == currCelltype) 
  currCells <- fdata$exprdat$samples[[currCelltype]] %>% dplyr::select(sample, sample_bk) %>% unique()
  currCells <- currCells %>% filter(sample_bk %in% currSamples$sample_bk)
  currExprmat <- fdata$exprmats[[currCelltype]][, currCells$sample]
  currExprmat <- currExprmat[sort(rownames(currExprmat)), ]
  currCoexmat <- workspace$utils$computeAggScCoex(currCells, currExprmat)
  
  # ==========================================
  # COMMIT ===================================
  print("Saving coexmat to disk...")
  currCoexmat %>% saveRDS(paste0(workspace$outputDir, paste0("sc_coexmats_", fdata$dataset, "_", currCelltype, ".rds")))
}, mc.cores = length(fdata$celltypes))


# +++++++++++++++++++++
# compute cross cell co-expression tables
# ramos

fdata <- list(); gc()

fdata$dataset <- "ramos"

fdata$exprdat <- readRDS(paste0(workspace$outputDir, "sc_exprmats_", fdata$dataset, ".rds"))

fdata$exprmats <- fdata$exprdat$exprmats

fdata$samples <- fdata$exprdat$samples %>% lapply(function(currSamples) {
  currSamples %>% dplyr::select(sample, sample_bk, cell_type) %>% group_by(cell_type, sample_bk) %>% summarize(n_cell = n()) %>% ungroup()
}) %>% session$dataWrangler$rbind()

fdata$samples <- fdata$samples %>% filter(cell_type %in% pdata$celltypes)

fdata$samples %>% 
  session$graphingUtils$ggplot(aes(x = sample_bk, y = n_cell)) + 
  geom_bar(stat = "identity") + 
  geom_hline(yintercept = 20) +
  scale_y_log10() + 
  session$graphingUtils$tiltX(angle = 90) + 
  facet_wrap(~cell_type, ncol = 1) + 
  geom_bar(data = fdata$samples %>% filter(n_cell < 20), fill = "red", stat = "identity")

fdata$samples %>% 
  filter(n_cell >= 20) %>% 
  group_by(cell_type) %>% 
  summarize(n_sbj = n())

fdata$samples <- fdata$samples %>% filter(n_cell >= 20) # filter out any sample with less than 20 cells

# also filter out any cell type with < 20 samples
fdata$samples <- fdata$samples %>% filter(cell_type %in% (fdata$samples %>%
                                                            group_by(cell_type) %>% 
                                                            summarize(n_sample_bk = n()) %>% 
                                                            filter(n_sample_bk >= 20) %>% 
                                                            session$dataWrangler$extractColumn("cell_type")))

fdata$celltypes <- pdata$celltypes %>% intersect(fdata$samples$cell_type)

# process each cell type and write to disk at each step so I don't hold a lot of memory
fdata$celltypes %>% mclapply(function(currCelltype) { 
  print(paste0("Working on cell type: ", currCelltype, "..."))
  currSamples <- fdata$samples %>% filter(cell_type == currCelltype) 
  currCells <- fdata$exprdat$samples[[currCelltype]] %>% dplyr::select(sample, sample_bk) %>% unique()
  currCells <- currCells %>% filter(sample_bk %in% currSamples$sample_bk)
  currExprmat <- fdata$exprmats[[currCelltype]][, currCells$sample]
  currExprmat <- currExprmat[sort(rownames(currExprmat)), ]
  currCoexmat <- workspace$utils$computeAggScCoex(currCells, currExprmat)
  
  # ==========================================
  # COMMIT ===================================
  print("Saving coexmat to disk...")
  currCoexmat %>% saveRDS(paste0(workspace$outputDir, paste0("sc_coexmats_", fdata$dataset, "_", currCelltype, ".rds")))
}, mc.cores = length(fdata$celltypes))

# +++++++++++++++++++++
# compute cross cell co-expression tables
# pineda

fdata <- list(); gc()

fdata$dataset <- "pineda"

fdata$exprdat <- readRDS(paste0(workspace$outputDir, "sc_exprmats_", fdata$dataset, ".rds"))

fdata$exprmats <- fdata$exprdat$exprmats

fdata$samples <- fdata$exprdat$samples %>% lapply(function(currSamples) {
  currSamples %>% dplyr::select(sample, sample_bk, cell_type) %>% group_by(cell_type, sample_bk) %>% summarize(n_cell = n()) %>% ungroup()
}) %>% session$dataWrangler$rbind()

fdata$samples <- fdata$samples %>% filter(cell_type %in% pdata$celltypes)

fdata$samples %>% 
  session$graphingUtils$ggplot(aes(x = sample_bk, y = n_cell)) + 
  geom_bar(stat = "identity") + 
  geom_hline(yintercept = 20) +
  scale_y_log10() + 
  session$graphingUtils$tiltX(angle = 90) + 
  facet_wrap(~cell_type, ncol = 1) + 
  geom_bar(data = fdata$samples %>% filter(n_cell < 20), fill = "red", stat = "identity")

fdata$samples %>% 
  filter(n_cell >= 20) %>% 
  group_by(cell_type) %>% 
  summarize(n_sbj = n())

fdata$samples <- fdata$samples %>% filter(n_cell >= 20) # filter out any sample with less than 20 cells

# also filter out any cell type with < 20 samples
fdata$samples <- fdata$samples %>% filter(cell_type %in% (fdata$samples %>%
                                                            group_by(cell_type) %>% 
                                                            summarize(n_sample_bk = n()) %>% 
                                                            filter(n_sample_bk >= 20) %>% 
                                                            session$dataWrangler$extractColumn("cell_type")))

fdata$celltypes <- pdata$celltypes %>% intersect(fdata$samples$cell_type)

# process each cell type and write to disk at each step so I don't hold a lot of memory
fdata$celltypes %>% mclapply(function(currCelltype) { 
  print(paste0("Working on cell type: ", currCelltype, "..."))
  currSamples <- fdata$samples %>% filter(cell_type == currCelltype) 
  currCells <- fdata$exprdat$samples[[currCelltype]] %>% dplyr::select(sample, sample_bk) %>% unique()
  currCells <- currCells %>% filter(sample_bk %in% currSamples$sample_bk)
  currExprmat <- fdata$exprmats[[currCelltype]][, currCells$sample]
  currExprmat <- currExprmat[sort(rownames(currExprmat)), ]
  currCoexmat <- workspace$utils$computeAggScCoex(currCells, currExprmat)
  
  # ==========================================
  # COMMIT ===================================
  print("Saving coexmat to disk...")
  currCoexmat %>% saveRDS(paste0(workspace$outputDir, paste0("sc_coexmats_", fdata$dataset, "_", currCelltype, ".rds")))
}, mc.cores = length(fdata$celltypes))



# +++++++++++++++++++++
# compute cross cell co-expression tables
# lim

fdata <- list(); gc()

fdata$dataset <- "lim"

fdata$exprdat <- readRDS(paste0(workspace$outputDir, "sc_exprmats_", fdata$dataset, ".rds"))

fdata$exprmats <- fdata$exprdat$exprmats

fdata$samples <- fdata$exprdat$samples %>% lapply(function(currSamples) {
  currSamples %>% dplyr::select(sample, sample_bk, cell_type) %>% group_by(cell_type, sample_bk) %>% summarize(n_cell = n()) %>% ungroup()
}) %>% session$dataWrangler$rbind()

fdata$samples <- fdata$samples %>% filter(cell_type %in% pdata$celltypes)

fdata$samples %>% 
  session$graphingUtils$ggplot(aes(x = sample_bk, y = n_cell)) + 
  geom_bar(stat = "identity") + 
  geom_hline(yintercept = 20) +
  scale_y_log10() + 
  session$graphingUtils$tiltX(angle = 90) + 
  facet_wrap(~cell_type, ncol = 1) + 
  geom_bar(data = fdata$samples %>% filter(n_cell < 20), fill = "red", stat = "identity")

fdata$samples %>% 
  filter(n_cell >= 20) %>% 
  group_by(cell_type) %>% 
  summarize(n_sbj = n())

fdata$samples <- fdata$samples %>% filter(n_cell >= 20) # filter out any sample with less than 20 cells

# also filter out any cell type with < 20 samples
fdata$samples <- fdata$samples %>% filter(cell_type %in% (fdata$samples %>%
                                                            group_by(cell_type) %>% 
                                                            summarize(n_sample_bk = n()) %>% 
                                                            filter(n_sample_bk >= 20) %>% 
                                                            session$dataWrangler$extractColumn("cell_type")))

fdata$celltypes <- pdata$celltypes %>% intersect(fdata$samples$cell_type)

# process each cell type and write to disk at each step so I don't hold a lot of memory
fdata$celltypes %>% mclapply(function(currCelltype) { 
  print(paste0("Working on cell type: ", currCelltype, "..."))
  currSamples <- fdata$samples %>% filter(cell_type == currCelltype) 
  currCells <- fdata$exprdat$samples[[currCelltype]] %>% dplyr::select(sample, sample_bk) %>% unique()
  currCells <- currCells %>% filter(sample_bk %in% currSamples$sample_bk)
  currExprmat <- fdata$exprmats[[currCelltype]][, currCells$sample]
  currExprmat <- currExprmat[sort(rownames(currExprmat)), ]
  currCoexmat <- workspace$utils$computeAggScCoex(currCells, currExprmat)
  
  # ==========================================
  # COMMIT ===================================
  print("Saving coexmat to disk...")
  currCoexmat %>% saveRDS(paste0(workspace$outputDir, paste0("sc_coexmats_", fdata$dataset, "_", currCelltype, ".rds")))
}, mc.cores = length(fdata$celltypes))

# +++++++++++++++++++++
# compute cross cell co-expression tables
# lau

fdata <- list(); gc()

fdata$dataset <- "lau"

fdata$exprdat <- readRDS(paste0(workspace$outputDir, "sc_exprmats_", fdata$dataset, ".rds"))

fdata$exprmats <- fdata$exprdat$exprmats

fdata$samples <- fdata$exprdat$samples %>% lapply(function(currSamples) {
  currSamples %>% dplyr::select(sample, sample_bk, cell_type) %>% group_by(cell_type, sample_bk) %>% summarize(n_cell = n()) %>% ungroup()
}) %>% session$dataWrangler$rbind()

fdata$samples <- fdata$samples %>% filter(cell_type %in% pdata$celltypes)

fdata$samples %>% 
  session$graphingUtils$ggplot(aes(x = sample_bk, y = n_cell)) + 
  geom_bar(stat = "identity") + 
  geom_hline(yintercept = 20) +
  scale_y_log10() + 
  session$graphingUtils$tiltX(angle = 90) + 
  facet_wrap(~cell_type, ncol = 1) + 
  geom_bar(data = fdata$samples %>% filter(n_cell < 20), fill = "red", stat = "identity")

fdata$samples %>% 
  filter(n_cell >= 20) %>% 
  group_by(cell_type) %>% 
  summarize(n_sbj = n())

fdata$samples <- fdata$samples %>% filter(n_cell >= 20) # filter out any sample with less than 20 cells

# also filter out any cell type with < 20 samples
fdata$samples <- fdata$samples %>% filter(cell_type %in% (fdata$samples %>%
                                                            group_by(cell_type) %>% 
                                                            summarize(n_sample_bk = n()) %>% 
                                                            filter(n_sample_bk >= 20) %>% 
                                                            session$dataWrangler$extractColumn("cell_type")))

fdata$celltypes <- pdata$celltypes %>% intersect(fdata$samples$cell_type)

# process each cell type and write to disk at each step so I don't hold a lot of memory
fdata$celltypes %>% mclapply(function(currCelltype) { 
  print(paste0("Working on cell type: ", currCelltype, "..."))
  currSamples <- fdata$samples %>% filter(cell_type == currCelltype) 
  currCells <- fdata$exprdat$samples[[currCelltype]] %>% dplyr::select(sample, sample_bk) %>% unique()
  currCells <- currCells %>% filter(sample_bk %in% currSamples$sample_bk)
  currExprmat <- fdata$exprmats[[currCelltype]][, currCells$sample]
  currExprmat <- currExprmat[sort(rownames(currExprmat)), ]
  currCoexmat <- workspace$utils$computeAggScCoex(currCells, currExprmat)
  
  # ==========================================
  # COMMIT ===================================
  print("Saving coexmat to disk...")
  currCoexmat %>% saveRDS(paste0(workspace$outputDir, paste0("sc_coexmats_", fdata$dataset, "_", currCelltype, ".rds")))
}, mc.cores = length(fdata$celltypes))


# ==============================
# ==============================
# ==============================
# SUBJECT LEVEL COEXPRESSION NETWORKS NOW!
# ==============================
# ==============================
# ==============================

pdata$files$coexmats <- tibble(outfile = list.files(workspace$outputDir)) %>% 
  filter(grepl("coexmats", outfile)) %>% 
  mutate(outfile_name = gsub(".rds", "", outfile)) %>% 
  mutate(strs = strsplit(outfile_name, "_")) %>% 
  mutate(level = strs %>% sapply(function(currStr) { currStr[1] })) %>% 
  mutate(type = strs %>% sapply(function(currStr) { currStr[2] })) %>% 
  mutate(dataset = strs %>% sapply(function(currStr) { currStr[3] })) %>% 
  mutate(cell_type = strs %>% sapply(function(currStr) { currStr[4] })) %>% 
  dplyr::select(outfile, level, type, dataset, cell_type)

# pdata$files$coexmats %>% filter(level == "sbj") %>% View()

pdata$files$exprmats %>% filter(level == "sbj")

# +++++++++++++++++++++
# compute cross subject co-expression tables
# velmeshev first

fdata <- list(); gc()

fdata$dataset <- "velmeshev"

fdata$exprdat <-  readRDS(paste0(workspace$outputDir, "sbj_exprmats_", fdata$dataset, ".rds"))

fdata$exprmats <- fdata$exprdat$exprmats

fdata$samples <- fdata$exprdat$samples %>% filter(n_cell > 20)

fdata$celltypes <- fdata$samples %>% 
  dplyr::select(cell_type, sample_bk) %>% unique() %>% group_by(cell_type) %>% summarize(n_sbj = n()) %>% 
  filter(n_sbj >= 20)

fdata$celltypes <- pdata$celltypes %>% intersect(fdata$celltypes$cell_type)  
names(fdata$celltypes) <- fdata$celltypes

fdata$celltypes %>% mclapply(function(currCelltype) {
  print(paste0("Working on cell type: ", currCelltype, "..."))
  currSamples <- fdata$samples %>% filter(cell_type == currCelltype) %>% arrange(sample)
  currExprmat <- fdata$exprmats[[currCelltype]][, currSamples$sample]
  currExprmat <- currExprmat[sort(rownames(currExprmat)), ]
  currCoexmat <- currExprmat %>% workspace$utils$computeCoexmat()
  currCoexmat <- currCoexmat %>% workspace$utils$normalizeCoexmat()
  
  # ==========================================
  # COMMIT ===================================
  print("Saving coexmat to disk...")
  currCoexmat %>% saveRDS(paste0(workspace$outputDir, paste0("sbj_coexmats_", fdata$dataset, "_", currCelltype, ".rds")))
}, mc.cores = length(fdata$celltypes))


# +++++++++++++++++++++
# compute cross subject co-expression tables
# rosmap

fdata <- list(); gc()

fdata$dataset <- "rosmap"

fdata$exprdat <-  readRDS(paste0(workspace$outputDir, "sbj_exprmats_", fdata$dataset, ".rds"))

fdata$exprmats <- fdata$exprdat$exprmats

fdata$samples <- fdata$exprdat$samples %>% filter(n_cell > 20)

(fdata$celltypes <- fdata$samples %>% 
  dplyr::select(cell_type, sample_bk) %>% unique() %>% group_by(cell_type) %>% summarize(n_sbj = n()) %>% 
  filter(n_sbj >= 20))

fdata$celltypes <- pdata$celltypes %>% intersect(fdata$celltypes$cell_type)  
names(fdata$celltypes) <- fdata$celltypes

fdata$celltypes %>% mclapply(function(currCelltype) {
  print(paste0("Working on cell type: ", currCelltype, "..."))
  currSamples <- fdata$samples %>% filter(cell_type == currCelltype) %>% arrange(sample)
  currExprmat <- fdata$exprmats[[currCelltype]][, currSamples$sample]
  currExprmat <- currExprmat[sort(rownames(currExprmat)), ]
  currCoexmat <- currExprmat %>% workspace$utils$computeCoexmat()
  currCoexmat <- currCoexmat %>% workspace$utils$normalizeCoexmat()
  
  # ==========================================
  # COMMIT ===================================
  print("Saving coexmat to disk...")
  currCoexmat %>% saveRDS(paste0(workspace$outputDir, paste0("sbj_coexmats_", fdata$dataset, "_", currCelltype, ".rds")))
}, mc.cores = length(fdata$celltypes))


# +++++++++++++++++++++
# compute cross subject co-expression tables
# ramos

fdata <- list(); gc()

fdata$dataset <- "ramos"

fdata$exprdat <-  readRDS(paste0(workspace$outputDir, "sbj_exprmats_", fdata$dataset, ".rds"))

fdata$exprmats <- fdata$exprdat$exprmats

fdata$samples <- fdata$exprdat$samples %>% filter(n_cell > 20)

(fdata$celltypes <- fdata$samples %>% 
    dplyr::select(cell_type, sample_bk) %>% unique() %>% group_by(cell_type) %>% summarize(n_sbj = n()) %>% 
    filter(n_sbj >= 20))

fdata$celltypes <- pdata$celltypes %>% intersect(fdata$celltypes$cell_type)  
names(fdata$celltypes) <- fdata$celltypes

fdata$celltypes %>% mclapply(function(currCelltype) {
  print(paste0("Working on cell type: ", currCelltype, "..."))
  currSamples <- fdata$samples %>% filter(cell_type == currCelltype) %>% arrange(sample)
  currExprmat <- fdata$exprmats[[currCelltype]][, currSamples$sample]
  currExprmat <- currExprmat[sort(rownames(currExprmat)), ]
  currCoexmat <- currExprmat %>% workspace$utils$computeCoexmat()
  currCoexmat <- currCoexmat %>% workspace$utils$normalizeCoexmat()
  
  # ==========================================
  # COMMIT ===================================
  print("Saving coexmat to disk...")
  currCoexmat %>% saveRDS(paste0(workspace$outputDir, paste0("sbj_coexmats_", fdata$dataset, "_", currCelltype, ".rds")))
}, mc.cores = length(fdata$celltypes))


# +++++++++++++++++++++
# compute cross subject co-expression tables
# pineda

fdata <- list(); gc()

fdata$dataset <- "pineda"

fdata$exprdat <-  readRDS(paste0(workspace$outputDir, "sbj_exprmats_", fdata$dataset, ".rds"))

fdata$exprmats <- fdata$exprdat$exprmats

fdata$samples <- fdata$exprdat$samples %>% filter(n_cell > 20)

(fdata$celltypes <- fdata$samples %>% 
    dplyr::select(cell_type, sample_bk) %>% unique() %>% group_by(cell_type) %>% summarize(n_sbj = n()) %>% 
    filter(n_sbj >= 20))

fdata$celltypes <- pdata$celltypes %>% intersect(fdata$celltypes$cell_type)  
names(fdata$celltypes) <- fdata$celltypes

fdata$celltypes %>% mclapply(function(currCelltype) {
  print(paste0("Working on cell type: ", currCelltype, "..."))
  currSamples <- fdata$samples %>% filter(cell_type == currCelltype) %>% arrange(sample)
  currExprmat <- fdata$exprmats[[currCelltype]][, currSamples$sample]
  currExprmat <- currExprmat[sort(rownames(currExprmat)), ]
  currCoexmat <- currExprmat %>% workspace$utils$computeCoexmat()
  currCoexmat <- currCoexmat %>% workspace$utils$normalizeCoexmat()
  
  # ==========================================
  # COMMIT ===================================
  print("Saving coexmat to disk...")
  currCoexmat %>% saveRDS(paste0(workspace$outputDir, paste0("sbj_coexmats_", fdata$dataset, "_", currCelltype, ".rds")))
}, mc.cores = length(fdata$celltypes))


# +++++++++++++++++++++
# compute cross subject co-expression tables
# nagy

fdata <- list(); gc()

fdata$dataset <- "nagy"

fdata$exprdat <-  readRDS(paste0(workspace$outputDir, "sbj_exprmats_", fdata$dataset, ".rds"))

fdata$exprmats <- fdata$exprdat$exprmats

fdata$samples <- fdata$exprdat$samples %>% filter(n_cell > 20)

(fdata$celltypes <- fdata$samples %>% 
    dplyr::select(cell_type, sample_bk) %>% unique() %>% group_by(cell_type) %>% summarize(n_sbj = n()) %>% 
    filter(n_sbj >= 20))

fdata$celltypes <- pdata$celltypes %>% intersect(fdata$celltypes$cell_type)  
names(fdata$celltypes) <- fdata$celltypes

fdata$celltypes %>% mclapply(function(currCelltype) {
  print(paste0("Working on cell type: ", currCelltype, "..."))
  currSamples <- fdata$samples %>% filter(cell_type == currCelltype) %>% arrange(sample)
  currExprmat <- fdata$exprmats[[currCelltype]][, currSamples$sample]
  currExprmat <- currExprmat[sort(rownames(currExprmat)), ]
  currCoexmat <- currExprmat %>% workspace$utils$computeCoexmat()
  currCoexmat <- currCoexmat %>% workspace$utils$normalizeCoexmat()
  
  # ==========================================
  # COMMIT ===================================
  print("Saving coexmat to disk...")
  currCoexmat %>% saveRDS(paste0(workspace$outputDir, paste0("sbj_coexmats_", fdata$dataset, "_", currCelltype, ".rds")))
}, mc.cores = length(fdata$celltypes))


# +++++++++++++++++++++
# compute cross subject co-expression tables
# lim

fdata <- list(); gc()

fdata$dataset <- "lim"

fdata$exprdat <-  readRDS(paste0(workspace$outputDir, "sbj_exprmats_", fdata$dataset, ".rds"))

fdata$exprmats <- fdata$exprdat$exprmats

fdata$samples <- fdata$exprdat$samples %>% filter(n_cell > 20)

(fdata$celltypes <- fdata$samples %>% 
    dplyr::select(cell_type, sample_bk) %>% unique() %>% group_by(cell_type) %>% summarize(n_sbj = n()) %>% 
    filter(n_sbj >= 20))

fdata$celltypes <- pdata$celltypes %>% intersect(fdata$celltypes$cell_type)  
names(fdata$celltypes) <- fdata$celltypes

fdata$celltypes %>% mclapply(function(currCelltype) {
  print(paste0("Working on cell type: ", currCelltype, "..."))
  currSamples <- fdata$samples %>% filter(cell_type == currCelltype) %>% arrange(sample)
  currExprmat <- fdata$exprmats[[currCelltype]][, currSamples$sample]
  currExprmat <- currExprmat[sort(rownames(currExprmat)), ]
  currCoexmat <- currExprmat %>% workspace$utils$computeCoexmat()
  currCoexmat <- currCoexmat %>% workspace$utils$normalizeCoexmat()
  
  # ==========================================
  # COMMIT ===================================
  print("Saving coexmat to disk...")
  currCoexmat %>% saveRDS(paste0(workspace$outputDir, paste0("sbj_coexmats_", fdata$dataset, "_", currCelltype, ".rds")))
}, mc.cores = length(fdata$celltypes))



# +++++++++++++++++++++
# compute cross subject co-expression tables
# lau

fdata <- list(); gc()

fdata$dataset <- "lau"

fdata$exprdat <-  readRDS(paste0(workspace$outputDir, "sbj_exprmats_", fdata$dataset, ".rds"))

fdata$exprmats <- fdata$exprdat$exprmats

fdata$samples <- fdata$exprdat$samples %>% filter(n_cell > 20)

(fdata$celltypes <- fdata$samples %>% 
    dplyr::select(cell_type, sample_bk) %>% unique() %>% group_by(cell_type) %>% summarize(n_sbj = n()) %>% 
    filter(n_sbj >= 20))

fdata$celltypes <- pdata$celltypes %>% intersect(fdata$celltypes$cell_type)  
names(fdata$celltypes) <- fdata$celltypes

fdata$celltypes %>% mclapply(function(currCelltype) {
  print(paste0("Working on cell type: ", currCelltype, "..."))
  currSamples <- fdata$samples %>% filter(cell_type == currCelltype) %>% arrange(sample)
  currExprmat <- fdata$exprmats[[currCelltype]][, currSamples$sample]
  currExprmat <- currExprmat[sort(rownames(currExprmat)), ]
  currCoexmat <- currExprmat %>% workspace$utils$computeCoexmat()
  currCoexmat <- currCoexmat %>% workspace$utils$normalizeCoexmat()
  
  # ==========================================
  # COMMIT ===================================
  print("Saving coexmat to disk...")
  currCoexmat %>% saveRDS(paste0(workspace$outputDir, paste0("sbj_coexmats_", fdata$dataset, "_", currCelltype, ".rds")))
}, mc.cores = length(fdata$celltypes))


# +++++++++++++++++++++
# compute cross bulk sample co-expression tables
# velmeshev

fdata <- list(); gc()

fdata$exprdat <-  readRDS(paste0(workspace$outputDir, "bk_exprmats_velmeshev.rds"))

fdata$exprmat <- fdata$exprdat$exprmat
fdata$exprmat <- fdata$exprmat[sort(rownames(fdata$exprmat)), ]

fdata$coefmat <- fdata$exprmat %>% t() %>% cor() # coefmat is the raw pearson's correlation version; used to compare changes in absolute correlation, disregarding relativity with other genes
fdata$coexmat <- fdata$coefmat %>% workspace$utils$normalizeCoexmat() # coexmat is the normalized version

# =================================  
# COMMIT ==========================
fdata$coefmat %>% saveRDS(paste0(workspace$outputDir, "bk_coefmats_velmeshev.rds"))
fdata$coexmat %>% saveRDS(paste0(workspace$outputDir, "bk_coexmats_velmeshev.rds"))
# =================================

# +++++++++++++++++++++
# compute cross bulk sample co-expression tables
# ROSMAP

fdata <- list(); gc()

fdata$exprdat <-  readRDS(paste0(workspace$outputDir, "bk_exprmats_rosmap.rds"))

fdata$exprmat <- fdata$exprdat$exprmat
fdata$exprmat <- fdata$exprmat[sort(rownames(fdata$exprmat)), ]

fdata$coefmat <- fdata$exprmat %>% t() %>% cor() # coefmat is the raw pearson's correlation version; used to compare changes in absolute correlation, disregarding relativity with other genes
fdata$coexmat <- fdata$coefmat %>% workspace$utils$normalizeCoexmat()

# =================================  
# COMMIT ==========================
fdata$coefmat %>% saveRDS(paste0(workspace$outputDir, "bk_coefmats_rosmap.rds"))
fdata$coexmat %>% saveRDS(paste0(workspace$outputDir, "bk_coexmats_rosmap.rds"))
# =================================

# +++++++++++++++++++++
# compute cross bulk sample co-expression tables
# ROSMAP subset for those with ihc cell type proportion estimates

fdata <- list(); gc()

fdata$exprdat <-  readRDS(paste0(workspace$outputDir, "bk_exprmats_rosmap-ihc.rds"))

fdata$exprmat <- fdata$exprdat$exprmat
fdata$exprmat <- fdata$exprmat[sort(rownames(fdata$exprmat)), ]

fdata$coefmat <- fdata$exprmat %>% t() %>% cor() # coefmat is the raw pearson's correlation version; used to compare changes in absolute correlation, disregarding relativity with other genes
fdata$coexmat <- fdata$coefmat %>% workspace$utils$normalizeCoexmat()

# =================================  
# COMMIT ==========================
fdata$coefmat %>% saveRDS(paste0(workspace$outputDir, "bk_coefmats_rosmap-ihc.rds"))
fdata$coexmat %>% saveRDS(paste0(workspace$outputDir, "bk_coexmats_rosmap-ihc.rds"))
# =================================

# +++++++++++++++++++++
# compute cross bulk sample co-expression tables
# ROSMAP with IHC cell type proprtions regressed out

fdata <- list(); gc()

fdata$exprdat <-  readRDS(paste0(workspace$outputDir, "bk_exprmats_rosmap-ihcres.rds"))

fdata$exprmat <- fdata$exprdat$exprmat
fdata$exprmat <- fdata$exprmat[sort(rownames(fdata$exprmat)), ]

fdata$coefmat <- fdata$exprmat %>% t() %>% cor() # coefmat is the raw pearson's correlation version; used to compare changes in absolute correlation, disregarding relativity with other genes
fdata$coexmat <- fdata$coefmat %>% workspace$utils$normalizeCoexmat()

# =================================  
# COMMIT ==========================
fdata$coefmat %>% saveRDS(paste0(workspace$outputDir, "bk_coefmats_rosmap-ihcres.rds"))
fdata$coexmat %>% saveRDS(paste0(workspace$outputDir, "bk_coexmats_rosmap-ihcres.rds"))
# =================================


# +++++++++++++++++++++
# compute cross bulk sample co-expression tables
# ROSMAP with MGPs regressed out

fdata <- list(); gc()

fdata$exprdat <-  readRDS(paste0(workspace$outputDir, "bk_exprmats_rosmap-mgpres.rds"))$exprdat

fdata$exprmat <- fdata$exprdat$exprmat
fdata$exprmat <- fdata$exprmat[sort(rownames(fdata$exprmat)), ]

fdata$coefmat <- fdata$exprmat %>% t() %>% cor() # coefmat is the raw pearson's correlation version; used to compare changes in absolute correlation, disregarding relativity with other genes
fdata$coexmat <- fdata$coefmat %>% workspace$utils$normalizeCoexmat()

# =================================  
# COMMIT ==========================
fdata$coefmat %>% saveRDS(paste0(workspace$outputDir, "bk_coefmats_rosmap-mgppres.rds"))
fdata$coexmat %>% saveRDS(paste0(workspace$outputDir, "bk_coexmats_rosmap-mgppres.rds"))
# =================================


# +++++++++++++++++++++
# compute cross bulk sample co-expression tables
# Velmeshev with MGPs regressed out

fdata <- list(); gc()

fdata$exprdat <-  readRDS(paste0(workspace$outputDir, "bk_exprmats_velmeshev-mgpres.rds"))$exprdat

fdata$exprmat <- fdata$exprdat$exprmat
fdata$exprmat <- fdata$exprmat[sort(rownames(fdata$exprmat)), ]

fdata$coefmat <- fdata$exprmat %>% t() %>% cor() # coefmat is the raw pearson's correlation version; used to compare changes in absolute correlation, disregarding relativity with other genes
fdata$coexmat <- fdata$coefmat %>% workspace$utils$normalizeCoexmat()

# =================================  
# COMMIT ==========================
fdata$coefmat %>% saveRDS(paste0(workspace$outputDir, "bk_coefmats_velmeshev-mgpres.rds"))
fdata$coexmat %>% saveRDS(paste0(workspace$outputDir, "bk_coexmats_velmeshev-mgpres.rds"))
# =================================

# compute co-expression networks for the two psdkbkres networks... 
# psdkbk vs psdbkres in both datasets

# +++++++++++++++++++++
# Velmeshev psdbk vs psdbkres

fdata <- list(); gc()

fdata$exprdat <-  readRDS(paste0(workspace$outputDir, "bk_exprmats_velmeshev-psdbkres.rds"))

fdata$coexmats <- fdata$exprdat$exprmat[c("bkOrig", "bkRes")] %>% session$collectionUtils$lapplyWithName(function(id, exprmat) {
  exprmatCopy <- exprmat[sort(rownames(exprmat)), ]
  print(paste0("computing coexpression network for ", id))
  coefmat <- exprmatCopy %>% t() %>% cor() # coefmat is the raw pearson's correlation version; used to compare changes in absolute correlation, disregarding relativity with other genes
  print(paste0("normalizing coexpression network for ", id))
  coexmat <- coefmat %>% workspace$utils$normalizeCoexmat()
  return(list(coefmat = coefmat, coexmat = coexmat))
})

# =================================  
# COMMIT ==========================
fdata$coexmats$bkOrig$coefmat %>% saveRDS(paste0(workspace$outputDir, "bk_coefmats_velmeshev-psdbk.rds"))
fdata$coexmats$bkOrig$coexmat %>% saveRDS(paste0(workspace$outputDir, "bk_coexmats_velmeshev-psdbk.rds"))
fdata$coexmats$bkRes$coefmat %>% saveRDS(paste0(workspace$outputDir, "bk_coefmats_velmeshev-psdbkres.rds"))
fdata$coexmats$bkRes$coexmat %>% saveRDS(paste0(workspace$outputDir, "bk_coexmats_velmeshev-psdbkres.rds"))
# =================================

# +++++++++++++++++++++
# ROSMAP psdbk vs psdbkres

fdata <- list(); gc()

fdata$exprdat <-  readRDS(paste0(workspace$outputDir, "bk_exprmats_rosmap-psdbkres.rds"))

fdata$coexmats <- fdata$exprdat$exprmat[c("bkOrig", "bkRes")] %>% session$collectionUtils$lapplyWithName(function(id, exprmat) {
  exprmatCopy <- exprmat[sort(rownames(exprmat)), ]
  print(paste0("computing coexpression network for ", id))
  coefmat <- exprmatCopy %>% t() %>% cor() # coefmat is the raw pearson's correlation version; used to compare changes in absolute correlation, disregarding relativity with other genes
  print(paste0("normalizing coexpression network for ", id))
  coexmat <- coefmat %>% workspace$utils$normalizeCoexmat()
  return(list(coefmat = coefmat, coexmat = coexmat))
})

# =================================  
# COMMIT ==========================
fdata$coexmats$bkOrig$coefmat %>% saveRDS(paste0(workspace$outputDir, "bk_coefmats_rosmap-psdbk.rds"))
fdata$coexmats$bkOrig$coexmat %>% saveRDS(paste0(workspace$outputDir, "bk_coexmats_rosmap-psdbk.rds"))
fdata$coexmats$bkRes$coefmat %>% saveRDS(paste0(workspace$outputDir, "bk_coefmats_rosmap-psdbkres.rds"))
fdata$coexmats$bkRes$coexmat %>% saveRDS(paste0(workspace$outputDir, "bk_coexmats_rosmap-psdbkres.rds"))
# =================================
























