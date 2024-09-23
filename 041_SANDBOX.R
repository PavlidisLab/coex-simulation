
# August 8, 2024
# let's examine the genes in each cluster
# maybe sort them by node degree in each cluster? As a method for measuring "importance"

source("/home/echu113/ws_rstudio/_workspaces/041_CH4_FINAL/041_SETUP.R")

pdata <- list()

pdata$genes <- read_rds(paste0(workspace$outputDir, "genes_metadata.rds"))


# read in the clusters & annotate with genes

pdata$clusters <- read.table(paste0(workspace$sdataDir, "clusters_genes.tsv"), sep = "\t", header = TRUE) %>% as_tibble()
pdata$clusters <- pdata$clusters %>% left_join(pdata$genes %>% dplyr::select(gene_id, gene))
pdata$clusters <- pdata$clusters %>% 
  mutate(ct_origin = strsplit(cluster_id, "\\.")) %>% 
  mutate(ct_origin = ct_origin %>% sapply(function(str) { str[1] }))
pdata$clusters <- pdata$clusters %>% dplyr::select(ct_origin, cluster_id, gene_id, gene)


# load cluster integrity

pdata$integrity <- read.table(paste0(workspace$sdataDir, "clusters_integrity.tsv"), sep = "\t", header = TRUE) %>% as_tibble()
pdata$clustersRnk <- pdata$integrity %>% filter(level == "xbulk", dataset %in% c("rosmap", "velmeshev")) %>% group_by(cluster) %>% summarize(auc = max(auc)) %>% arrange(desc(auc))
pdata$clustersTop <- pdata$clustersRnk$cluster[1:20]
pdata$clustersBtm <- pdata$clustersRnk$cluster[96:115]

# now, let's work out the node degree for each cluster of genes, in each consensus network

fdata <- list()

fdata$clustersFlat <- pdata$clusters$cluster_id %>% unique()

fdata$celltypes <- pdata$clusters$ct_origin %>% unique() %>% sort()
names(fdata$celltypes) <- fdata$celltypes

fdata$filenames <- paste0("xcell_coexmat_", fdata$celltype, "_consensus.tsv")
names(fdata$filenames) <- fdata$celltype

fdata$main <- fdata$filenames %>% session$collectionUtils$lapplyWithName(function(celltype, file) {
 
  coexmat <- read.table(paste0(workspace$sdataDir, file), sep = "\t", header = TRUE)
  coexmat <- coexmat %>% session$dataWrangler$setColAsRownames("gene_id") %>% as.matrix()
  
  fdata$clustersFlat %>% session$collectionUtils$lapply(function(clusterId) {
    clusterGenes <- pdata$clusters %>% filter(cluster_id == clusterId) %>% session$dataWrangler$extractColumn("gene_id")
    clusterGenes <- clusterGenes %>% intersect(rownames(coexmat))
    
    if (length(clusterGenes) < 2) {
      return(tibble(ct_origin = c(), cluster_id = c(), gene_id = c(), gene = c(), degree = c(), ct_coex = c()))
    }
    
    result <- coexmat[clusterGenes, clusterGenes] %>%
      apply(1, sum) %>% 
      session$dataWrangler$vectorToTibble() %>% 
      dplyr::select(gene_id = variable, degree = value)
    
    pdata$clusters %>% 
      filter(cluster_id == clusterId) %>% 
      inner_join(result, by = "gene_id") %>% 
      mutate(ct_coex = celltype)
    
  }) %>% session$dataWrangler$rbind()
  
}) %>% session$dataWrangler$rbind()


fdata$main %>% filter(cluster_id == "astrocyte.id_0007") %>% 
  group_by(ct_coex) %>% 
  mutate(score = degree / n()) %>% 
  session$graphingUtils$ggplot(aes(x = ct_coex, y = score)) + 
  geom_violin() + 
  geom_point(position = "jitter") +
  session$graphingUtils$tiltX(angle = 90)


xdata <- list()

xdata$main <- fdata$main %>% filter(cluster_id %in% "excitatory.id_0020") %>% 
  group_by(cluster_id, ct_coex) %>% 
  mutate(score = degree / n()) %>% 
  ungroup()

xdata$geneLabs <- xdata$main %>% group_by(ct_coex) %>% mutate(score_rank = rank(-score)) %>% filter(score_rank <= 5) %>% ungroup()

xdata$main %>% 
  session$graphingUtils$ggplot(aes(x = ct_coex, y = score)) + 
  geom_violin() + 
  geom_point(position = "jitter") + 
  geom_label_repel(data = xdata$geneLabs, aes(label = gene), ) +
  session$graphingUtils$tiltX(angle = 90) + 
  facet_wrap(~cluster_id, ncol = 3) + 
  geom_hline(yintercept = 0.5, linetype = "dashed")

(xdata$geneLabs %>% filter(ct_coex == "excitatory") %>% arrange(desc(score)) -> x)
x$gene


fdata$main %>% filter(cluster_id == "excitatory.id_0028") %>% arrange(desc(degree))

fdata$main %>% filter(cluster_id == "excitatory.id_0028") %>% filter(gene == "KHDRBS3")

# generate a table with 5 genes!

fdata$main %>% 
  filter(ct_coex == ct_origin) %>% 
  group_by(cluster_id) %>% 
  mutate(degree_rnk = rank(-degree)) %>% 
  ungroup() %>% 
  filter(degree_rnk <= 5) %>% 
  filter(cluster_id %in% pdata$clustersTop) %>% 
  group_by(cluster_id) %>% 
  summarize(genes = str_c(gene, collapse = ", ")) %>% 
  left_join(tibble(cluster_id = pdata$clustersTop, i = 1:length(pdata$clustersTop))) %>% 
  arrange(i) %>% 
  dplyr::select(-i) %>% View()

fdata$main %>% 
  filter(ct_coex == ct_origin) %>% 
  group_by(cluster_id) %>% 
  mutate(degree_rnk = rank(-degree)) %>% 
  ungroup() %>% 
  filter(degree_rnk <= 5) %>% 
  filter(cluster_id %in% pdata$clustersBtm) %>% 
  group_by(cluster_id) %>% 
  summarize(genes = str_c(gene, collapse = ", ")) %>% 
  left_join(tibble(cluster_id = pdata$clustersBtm, i = 1:length(pdata$clustersBtm))) %>% 
  arrange(i) %>% 
  dplyr::select(-i) %>% View()


# get cluster sizes

pdata$clusters %>% group_by(cluster_id) %>% summarize(n_gene = n()) %>% 
  filter(cluster_id %in% pdata$clustersTop) %>% 
  group_by(cluster_id) %>%
  left_join(tibble(cluster_id = pdata$clustersTop, i = 1:length(pdata$clustersTop))) %>% 
  arrange(i) %>% 
  dplyr::select(-i) %>% View()

pdata$clusters %>% group_by(cluster_id) %>% summarize(n_gene = n()) %>% 
  filter(cluster_id %in% pdata$clustersBtm) %>% 
  group_by(cluster_id) %>%
  left_join(tibble(cluster_id = pdata$clustersBtm, i = 1:length(pdata$clustersBtm))) %>% 
  arrange(i) %>% 
  dplyr::select(-i) %>% View()
  




coexmat <- read.table(paste0(workspace$sdataDir, fdata$filenames[["inhibitory"]]), sep = "\t", header = TRUE)
coexmat <- coexmat %>% session$dataWrangler$setColAsRownames("gene_id") %>% as.matrix()

currCluster <- pdata$clusters %>% filter(cluster_id == "inhibitory.id_0012") %>% session$dataWrangler$extractColumn("gene_id")

currCluster <- currCluster %>% intersect(rownames(coexmat))
coexmat[currCluster, currCluster] %>% apply(1, sum) %>% session$dataWrangler$vectorToTibble() %>% dplyr::select(gene_id = variable, degree = value) -> x
pdata$clusters %>% filter(cluster_id == "inhibitory.id_0012") %>% inner_join(x) %>% arrange(desc(degree)) -> x





# August 31 , 2023
# effect of log norm on simulated null 

# TODO - clean this up and talk about it over meeting once back ... 
# show a couple of examples as to why this occured. 

pdata <- list()

set.seed(2); pdata$exprmat$real <- rnorm(100000, 1000, 1000) %>% abs() %>% round() %>% matrix(nrow = 100)
rownames(pdata$exprmat$real) <- paste0("gene.", formatC(1:nrow(pdata$exprmat$real), digits = 2, flag = "0"))
colnames(pdata$exprmat$real) <- paste0("sample.", formatC(1:ncol(pdata$exprmat$real), digits = 3, flag = "0"))

pdata$exprmat$realfrac <- pdata$exprmat$real %>% apply(2, function(vals) { vals / sum(vals) })

set.seed(3); pdata$libSizes <- rnorm(1000, 1000, 1000) %>% abs() %>% round()
names(pdata$libSizes) <- colnames(pdata$exprmat$realfrac)

pdata$exprmat$ctmat <- 1:1000 %>% sapply(function(currI) { pdata$exprmat$realfrac[, currI] * pdata$libSizes[currI] }) %>% round()
colnames(pdata$exprmat$ctmat) <- names(pdata$libSizes)

# remove those with less than 30 reads in total
pdata$samplesV <- which((pdata$exprmat$ctmat %>% colSums() >= 30))

pdata$exprmat <- pdata$exprmat %>% lapply(function(exprmat) { exprmat[, pdata$samplesV] })

pdata$exprmat$cpm <- pdata$exprmat$ctmat %>% workspace$utils$ctmatToExprmat2()
pdata$exprmat$logcpm <- pdata$exprmat$ctmat %>% workspace$utils$ctmatToExprmat2() %>% workspace$utils$logTrans()

# remove really bad samples; those with library size 0
pdata$coexmat$pearson <- pdata$exprmat %>% session$collectionUtils$lapply(function(exprmat) { exprmat %>% t() %>% cor(method = "pearson") })
pdata$coexmat$spearman <- pdata$exprmat %>% session$collectionUtils$lapply(function(exprmat) { exprmat %>% t() %>% cor(method = "spearman") })

pdata$coextbl$pearson <- do.call(cbind, (pdata$coexmat$pearson %>% lapply(function(coexmat) { coexmat %>% workspace$utils$vectorize() }))) 
rownames(pdata$coextbl$pearson) <- pdata$coexmat$pearson$real %>% workspace$utils$getPairIds()

pdata$coextbl$spearman <- do.call(cbind, (pdata$coexmat$spearman %>% lapply(function(coexmat) { coexmat %>% workspace$utils$vectorize() }))) 
rownames(pdata$coextbl$spearman) <- pdata$coexmat$spearman$real %>% workspace$utils$getPairIds()

pdata$coextbl$pearson %>% as.data.frame() %>% GGally::ggpairs()
pdata$coextbl$spearman %>% as.data.frame() %>% GGally::ggpairs()

pdata$coextbl$pearson %>% session$dataWrangler$setRownameAsColumn("pair_id") %>% arrange(desc(logcpm)) # gene.093.gene.005
pdata$coextbl$pearson %>% session$dataWrangler$setRownameAsColumn("pair_id") %>% arrange(desc(real)) 

pdata$coextbl$pearson %>% session$dataWrangler$setRownameAsColumn("pair_id") %>% 
  session$graphingUtils$ggplot(aes(x = logcpm, y = ctmat)) + geom_point(shape = 1)

# change the matrix to see different effects
pdata$exprmat$ctmat[c("gene.073", "gene.025"),] %>% t() %>% session$dataWrangler$setRownameAsColumn("sample") %>% 
  setnames(c("sample", "gene_a", "gene_b")) %>% 
  session$graphingUtils$ggplot(aes(x = gene_a, y = gene_b)) + geom_point() + 
  geom_smooth(method = "lm", se = FALSE)

pdata$coextbl$pearson %>% session$dataWrangler$setRownameAsColumn("pair_id") %>% 
  gather("mat", "coex", -pair_id) %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + geom_density(aes(color = mat))

pdata$coextbl$pearson %>% session$dataWrangler$setRownameAsColumn("pair_id") %>% arrange(desc(logcpm)) -> x
pdata$exprmat$real %>% rowMeans() %>% session$dataWrangler$vectorToTibble() %>% dplyr::select(gene_id = variable, expr = value) -> y

x$pair_id %>% strsplit("\\.")


# let's examine coincidence with library size
xdata <- list()
xdata$samples <- names(pdata$libSizes) %>% intersect(colnames(pdata$exprmat$logcpm))
xdata$main <- pdata$exprmat$logcpm[, xdata$samples] %>% apply(1, function(vals) { cor(vals, pdata$libSizes[xdata$samples]) })
pdata$exprmat$real[, xdata$samples] %>% apply(1, mean) %>% session$dataWrangler$vectorToTibble() %>% 
  dplyr::select(gene = variable, expr = value) %>% 
  left_join(xdata$main %>% session$dataWrangler$vectorToTibble() %>% dplyr::select(gene = variable, libsize_bias = value)) %>% 
  session$graphingUtils$ggplot(aes(x = expr, y = libsize_bias)) + geom_point() + geom_smooth(method = "lm") + 
  ylim(-0.1, 0.3)



# August 2, 2023
# exploration of the clusters in velmeshev and rosmap post correction

pdata <- list()

pdata$coexmat <- readRDS(paste0(workspace$outputDir, "sc_coexmats_rosmap_inhibitory.rds"))
# pdata$coexmat <- readRDS(paste0(workspace$spaceDir, "output_20230731/", "sc_coexmats_velmeshev_inhibitory.rds"))

pdata$clusters <- pdata$coexmat %>% workspace$utils$computeClusters(function(dendro, distMat) {
  dynamicTreeCut::cutreeDynamic(dendro, distM = distMat, deepSplit = FALSE)
})

pdata$clusters %>% workspace$utils$plotDendro()

# July 26, 2023
# exploration of normalization effects on co-expression
# in this script, let's pick one sample in Velmeshev to perform co-expression analysis
# see if the "general co-expression" is reproduced
# see if gene expression profiles correlate with sequencing depth
# see if different normalization techniques get rid of it

pdata <- list()

pdata$exprdat <- readRDS(paste0(workspace$outputDir, "sc_exprmats_velmeshev.rds"))
pdata$samples <- pdata$exprdat$samples$inhibitory %>% filter(sample_bk == "5387_PFC")

pdata$ctmat <- pdata$exprdat$raw$ctmat[, pdata$samples$sample] %>% as.matrix()
pdata$ctmat <- pdata$ctmat %>% workspace$utils$cleanCtmat2()

# check the distribution of correlation of gene expression profiles with library size
pdata$genes <- rownames(pdata$ctmat)
names(pdata$genes) <- pdata$genes
pdata$genes <- read_rds(paste0(workspace$outputDir, "genes_metadata.rds")) %>% filter(gene_id %in% pdata$genes)
pdata$genes <- pdata$genes %>% arrange(gene_id)

# filter ctmat with the genes that are present
pdata$ctmat <- pdata$ctmat[pdata$genes$gene_id, ]

# work out libsizes
pdata$libsizes <- pdata$ctmat %>% apply(2, sum)

# here - try different versions of normalization
# 1. cpm
# 2. logcpm
# 3. logcpmres
# 4. frac
# 5. logfrac 
# 6. pf
# 7. logpf 
# 8. pflogpf

pdata$exprmat <- list()

# 1. cpm 

fdata <- list()

fdata$normalize <- function(ctmat) {
  libSizes <- ctmat %>% apply(2, sum)
  result <- t( t(ctmat) / libSizes) * 1e6 # convert the unit into CPM
  return(result)
}

fdata$exprmat <- pdata$ctmat %>% fdata$normalize()
pdata$exprmat$cpm <- fdata$exprmat 

# 2. logcpm

fdata <- list()

fdata$normalize <- function(ctmat) {
  libSizes <- ctmat %>% apply(2, sum)
  result <- t( t(ctmat) / libSizes) * 1e6 # convert the unit into CPM
  result <- log2(result + 1)
  return(result)
}

fdata$exprmat <- pdata$ctmat %>% fdata$normalize()
pdata$exprmat$logcpm <- fdata$exprmat 

# 3. logcpmres

fdata <- list()

fdata$normalize <- function(ctmat) {
  libSizes <- ctmat %>% apply(2, sum)
  result <- t( t(ctmat) / libSizes) * 1e6 # convert the unit into CPM
  result <- log2(result + 1)
  
  # fit it to library size
  genes <- pdata$genes$gene_id
  names(genes) <- genes
  
  result <- genes %>% sapply(function(currGene) {
    expr <- result[currGene, ]
    lm <- lm(expr ~ libSizes)
    return(lm$residuals)
  }) %>% t()
  
  return(result)
}

fdata$exprmat <- pdata$ctmat %>% fdata$normalize()
pdata$exprmat$logcpmres <- fdata$exprmat 

# 4. frac

fdata <- list()

fdata$normalize <- function(ctmat) {
  libSizes <- ctmat %>% apply(2, sum)
  result <- t( t(ctmat) / libSizes)
  return(result)
}

fdata$exprmat <- pdata$ctmat %>% fdata$normalize()
pdata$exprmat$frac <- fdata$exprmat 

# 5. logfrac 

fdata <- list()

fdata$normalize <- function(ctmat) {
  libSizes <- ctmat %>% apply(2, sum)
  result <- t( t(ctmat) / libSizes)
  result <- log2(result + 1)
  return(result)
}

fdata$exprmat <- pdata$ctmat %>% fdata$normalize()
pdata$exprmat$logfrac <- fdata$exprmat 

# 6. pf

fdata <- list()

fdata$normalize <- function(ctmat) {
  mtx <- t(ctmat)
  pf <- rowSums(mtx)
  mtx_pf <- mtx / (pf / mean(pf))
  mtx_pf <- t(mtx_pf)
  return(mtx_pf)
}

fdata$exprmat <- pdata$ctmat %>% fdata$normalize()
pdata$exprmat$pf <- fdata$exprmat 

# 7. logpf 

fdata <- list()

fdata$normalize <- function(ctmat) {
  mtx <- t(ctmat)
  pf <- rowSums(mtx)
  mtx_pf <- mtx / (pf / mean(pf))
  mtx_pf <- t(mtx_pf)
  mtx_pf <- log2(mtx_pf + 1)
  return(mtx_pf)
}

fdata$exprmat <- pdata$ctmat %>% fdata$normalize()
pdata$exprmat$logpf <- fdata$exprmat 

# 8. pflogpf

fdata <- list()

fdata$normalize <- function(ctmat) {
  
  doPf <- function(mtx) {
    pf <- rowSums(mtx)
    mtx_pf <- mtx / (pf / mean(pf))
    return(mtx_pf)
  }
  
  result <- doPf(t(ctmat))
  result <- log2(result + 1)
  result <- doPf(result)

  result <- t(result)
  
  return(result)
}

fdata$exprmat <- pdata$ctmat %>% fdata$normalize()
pdata$exprmat$pflogpf <- fdata$exprmat 


# ++++ first thing: check whether the gene expressions are correlated with library size
pdata$bias <- pdata$exprmat %>% session$collectionUtils$lapplyWithName(function(method, currExprmat) {
  tbl <- pdata$genes$gene_id %>% lapply(function(currGene) {
    tibble(gene_id = currGene, cor_coef = cor(currExprmat[currGene, ], pdata$libsizes, method = "pearson")[1, 1])
  }) %>% 
    session$dataWrangler$rbind() %>%
    mutate(method = method)
}) %>% session$dataWrangler$rbind()

pdata$bias %>% 
  filter(method != "logcpmres") %>% 
  session$graphingUtils$ggplot(aes(x = cor_coef, fill = method, color = method)) + 
  geom_density(alpha = 0.3) + 
  geom_point(y = 0, size = 3, alpha = 0.5) + 
  facet_wrap(~method, ncol = 1) + 
  geom_vline(xintercept = 0, linetype = "dashed")

# examine more closely how the expression level plots against the library size in different cases

fdata <- list()

fdata$genes <- pdata$bias %>% filter(method == "logcpm") %>% arrange(desc(cor_coef)) ; 
pdata$bias %>% filter(method == "logcpm") %>% arrange(desc(cor_coef))

# c("ENSG00000164742", "ENSG00000248049"); lots of drop outs
# c("ENSG00000251562", "ENSG00000224078"); no drop outs
# c("ENSG00000112293", "ENSG00000001460"); high co-expression

fdata <- list()

fdata$main <- pdata$exprmat$logcpm[c("ENSG00000146677", "ENSG00000144713"), ] %>% 
  t() %>% cbind(pdata$libsizes) %>% 
  session$dataWrangler$setRownameAsColumn("sample") %>% 
  dplyr::select(everything(), lib_size = V3)

fdata$main %>% 
  session$graphingUtils$ggplot(aes(x = lib_size, y = ENSG00000146677)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)

fdata$main %>% 
  session$graphingUtils$ggplot(aes(x = lib_size, y = ENSG00000144713)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)

fdata$main %>% 
  session$graphingUtils$ggplot(aes(x = ENSG00000146677, y = ENSG00000144713)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)

fdata$main %>% session$dataWrangler$setColAsRownames("sample") %>% cor(method = "pearson")


# ++++ second: compute and compare the co-expression networks

pdata$coexmat <- pdata$exprmat %>% mclapply(workspace$utils$computeCoexmat, mc.cores = length(pdata$exprmat)) 
pdata$coexmatRnk <- pdata$coexmat %>% mclapply(workspace$utils$normalizeCoexmat, mc.cores = length(pdata$exprmat))

# first, select 1000 random genes and look at the overall distribution of correlations

fdata <- list(); gc()

fdata$bias <- pdata$bias %>% filter(method == "logcpm") %>% arrange(desc(cor_coef))

set.seed(1); fdata$genes$rand <- pdata$genes$gene_id %>% sample(1000)
fdata$genes$rpgs <- pdata$genes %>% filter(grepl("^RP[S|L]*", gene)) %>% session$dataWrangler$extractColumn("gene_id")
fdata$genes$biased <- fdata$bias$gene_id[1:100]

fdata$coextbl <- pdata$coexmat %>% session$collectionUtils$lapplyWithName(function(method, currCoexmat) {
  fdata$genes %>% session$collectionUtils$lapplyWithName(function(setName, set) {
    currCoexmat[set, set] %>% workspace$utils$getCoords() %>% mutate(coex = currCoexmat[set, set] %>% workspace$utils$vectorize(), method = method, set = setName) 
  }, verbose = FALSE) %>% 
    session$dataWrangler$rbind()
}) %>% session$dataWrangler$rbind()

fdata$coextbl %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_density(aes(color = set, fill = set), alpha = 0.2) + 
  facet_wrap(~method, ncol = 1)

fdata$coextbl %>% 
  filter(set == "rpgs") %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_density(aes(color = method), alpha = 0.2) + 
  geom_vline(xintercept = 0, linetype = "dashed")

fdata$coexRnktbl <- pdata$coexmatRnk %>% session$collectionUtils$lapplyWithName(function(method, currCoexmat) {
  fdata$genes %>% session$collectionUtils$lapplyWithName(function(setName, set) {
    currCoexmat[set, set] %>% workspace$utils$getCoords() %>% mutate(coex = currCoexmat[set, set] %>% workspace$utils$vectorize(), method = method, set = setName) 
  }, verbose = FALSE) %>% 
    session$dataWrangler$rbind()
}) %>% session$dataWrangler$rbind()

fdata$coexRnktbl %>% 
  filter(set == "biased") %>% 
  session$graphingUtils$ggplot(aes(x = coex)) + 
  geom_density(aes(color = method), alpha = 0.2) + 
  facet_wrap(~method, nrow = 1)

# correlations of correlations !

xdata <- list()

xdata$main <- fdata$coextbl %>% filter(set == "rand") %>% 
  mutate(pair_id = paste0(gene_a, ".", gene_b)) %>%
  dplyr::select(pair_id, method, coex) %>% 
  spread(method, coex) %>% 
  session$dataWrangler$setColAsRownames("pair_id")

# do another round of subsampling
set.seed(0); xdata$indices <- 1:nrow(xdata$main) %>% sample(5000)

xdata$main[xdata$indices, ] %>% GGally::ggpairs()


# let's take a look at how consistent the top 1% are.... 

fdata <- list(); gc()

fdata$lnkmats <- pdata$coexmatRnk %>% session$collectionUtils$lapply(workspace$utils$computeLnksMat)

fdata$lnks <- fdata$lnkmats %>% session$collectionUtils$lapply(function(lnkmat) {
  lnks <- lnkmat %>% workspace$utils$getCoords() %>% mutate(lnk = lnkmat %>% workspace$utils$vectorize())
  lnks <- lnks %>% filter(lnk > 0)
  lnks <- lnks %>% mutate(pair_id = paste0(gene_a, ".", gene_b))
  return(sort(lnks$pair_id))
})

fdata$lnks %>% lapply(length)

fdata$overlaps <- fdata$lnks %>% session$collectionUtils$lapplyWithName(function(method_a, lnks_a) {
  fdata$lnks %>% session$collectionUtils$lapplyWithName(function(method_b, lnks_b) {
    nOverlap <- lnks_a %>% intersect(lnks_b) %>% length()
    fracOverlap <- nOverlap / length(lnks_a)
    tibble(method_a = method_a, method_b = method_b, n_intersect = nOverlap, frac_intersect = fracOverlap)
  }, verbose = FALSE) %>% session$dataWrangler$rbind()
}) %>% session$dataWrangler$rbind()

fdata$overlaps %>% 
  dplyr::select(method_a, method_b, frac_intersect) %>% 
  spread(method_b, frac_intersect) %>% 
  session$dataWrangler$setColAsRownames("method_a") %>%
  session$graphingUtils$heatmap(decimalNums = 2)



# next, let's look at the difference in correlation of the most biased genes

fdata <- list()


# let's check the top 10








# ok try 2 different normalization methods
pdata$logcpm <- pdata$ctmat %>% workspace$utils$ctmatToExprmat2(logtrans = TRUE)
pdata$logcpm2 <- pdata$logcpm %>% workspace$utils$ctmatToExprmat2(logtrans = FALSE)

pdata$logcpmDepthBias <- pdata$genes$gene_id %>% session$collectionUtils$lapply(function(currGene) { # ENSG00000164742, ENSG00000022355, ENSG00000114948
  tibble(gene_id = currGene, cor_coef = cor(fdata$exprmat[currGene, ], pdata$libsizes, method = "pearson")[1, 1])
}) %>% session$dataWrangler$rbind()
pdata$logcpmDepthBias %>% session$graphingUtils$ggplot(aes(x = cor_coef)) + geom_density()

pdata$logcpm2DepthBias <- pdata$genes$gene_id %>% session$collectionUtils$lapply(function(currGene) {
  tibble(gene_id = currGene, cor_coef = cor(pdata$logcpm2[currGene, ], pdata$libsizes, method = "pearson")[1, 1])
}) %>% session$dataWrangler$rbind()
pdata$logcpm2DepthBias %>% session$graphingUtils$ggplot(aes(x = cor_coef)) + geom_density()

# try regressing it out? 
pdata$logcpmRes <- fdata$genes %>% sapply(function(currGene) {
  expr <- pdata$logcpm[currGene, ]
  lm <- lm(expr ~ pdata$libsizes)
  return(lm$residuals)
}) %>% t()

pdata$logcpmResDepthBias <- fdata$genes %>% session$collectionUtils$lapply(function(currGene) {
  tibble(gene_id = currGene, cor_coef = cor(pdata$logcpmRes[currGene, ], pdata$libsizes)[1, 1])
}) %>% session$dataWrangler$rbind()
pdata$logcpmResDepthBias %>% session$graphingUtils$ggplot(aes(x = cor_coef)) + geom_histogram()


# TEMP CODE

pdata$logcpm["ENSG00000164742",] %>% session$dataWrangler$vectorToTibble() %>% dplyr::select(sample = variable, expr = value) %>% 
  left_join(pdata$logcpm %>% apply(2, sum) %>% session$dataWrangler$vectorToTibble() %>% dplyr::select(sample = variable, lib_size = value), by = "sample") %>% 
  session$graphingUtils$ggplot(aes(x = lib_size, y = expr)) + 
  geom_point(alpha = 0.1) +
  geom_smooth(method = "lm", se = FALSE)


pdata$logcpm %>% apply(2, sum) %>% session$dataWrangler$vectorToTibble() %>% dplyr::select(sample = variable, lib_size_logcpm = value) %>% 
  left_join(pdata$libsizes %>% session$dataWrangler$vectorToTibble() %>% dplyr::select(sample = variable, lib_size_orig = value)) -> x
  
x %>% session$graphingUtils$ggplot(aes(x = lib_size_orig, y = lib_size_logcpm)) + geom_point()



# make the co-expression matrics
pdata$logcpmCoexmat <- pdata$logcpm %>% workspace$utils$computeCoexmat() %>% workspace$utils$normalizeCoexmat()
pdata$logcpm2Coexmat <- pdata$logcpm2 %>% workspace$utils$computeCoexmat() %>% workspace$utils$normalizeCoexmat()
pdata$logcpmResCoexmat <- pdata$logcpmRes %>% workspace$utils$computeCoexmat() %>% workspace$utils$normalizeCoexmat()

pdata$logcpmClusters <- pdata$logcpmCoexmat %>% workspace$utils$computeClusters(function(dendro, distMat) {
  dynamicTreeCut::cutreeDynamic(dendro, distM = distMat, deepSplit = FALSE)
})

pdata$logcpm2Clusters <- pdata$logcpm2Coexmat %>% workspace$utils$computeClusters(function(dendro, distMat) {
  dynamicTreeCut::cutreeDynamic(dendro, distM = distMat, deepSplit = 0, cutHeight = 0.2)
})

pdata$logcpmResClusters <- pdata$logcpmResCoexmat %>% workspace$utils$computeClusters(function(dendro, distMat) {
  dynamicTreeCut::cutreeDynamic(dendro, distM = distMat, deepSplit = 0, cutHeight = 0.2)
})
  
# inspect the dendrograms

pdata$logcpmClusters %>% workspace$utils$plotDendro()
pdata$logcpm2Clusters %>% workspace$utils$plotDendro()
pdata$logcpmResClusters %>% workspace$utils$plotDendro()


workspace$utils$computeEgad(pdata$logcpmCoexmat, clusterMat = pdata$logcpm2Clusters$clusterMat)
workspace$utils$computeEgad(pdata$logcpmResCoexmat, clusterMat = pdata$logcpm2Clusters$clusterMat)
workspace$utils$computeEgad(pdata$logcpm2Coexmat, clusterMat = pdata$logcpm2Clusters$clusterMat)




plotDendroAndColors(pdata$logcpmClusters$dendro, 
                    colors = labels2colors(pdata$logcpmClusters$clusters), 
                    rowText = pdata$logcpmClusters$clusters, 
                    dendroLabels = FALSE)


pdata$ctmat[1:10, 1:10]





fdata <- list()
  
fdata$genes$rand <- pdata$genes$gene_id %>% sample(1000)
fdata$genes$rpgs <- pdata$genes %>% filter(grepl("^RP[S|L]*", gene)) %>% session$dataWrangler$extractColumn("gene_id")

fdata$coexmat$rand <- pdata$logcpmCoexmat[fdata$genes$rand, fdata$genes$rand]
fdata$coexmat$rpgs <- pdata$logcpmCoexmat[fdata$genes$rpgs, fdata$genes$rpgs]

pdata$logcpmCoextbl <- fdata$coexmat$rand %>% workspace$utils$getCoords() %>% mutate(coex = fdata$coexmat$rand %>% workspace$utils$vectorize(), type = "rand") %>% 
  rbind(fdata$coexmat$rpgs %>% workspace$utils$getCoords() %>% mutate(coex = fdata$coexmat$rpgs %>% workspace$utils$vectorize(), type = "rpgs"))

pdata$logcpmCoextbl %>% session$graphingUtils$ggplot(aes(x = coex)) + geom_density(aes(fill = type, alpha = 0.2))



pdata$logcpm[c("ENSG00000146677", "ENSG00000144713"), ] %>% t() %>% as.data.frame() %>% GGally::ggpairs()


# for each normalization method
# 1. how does each gene correlate with the library size - (original)?; both the distribution and a few examples
# 2. range of co-expression values (before & after rank transformation)
# 3. dendrogram + cluster distributions
# 4. EGAD of cluster qualities
# 5. Co-expression of ribosomal protein genes... 




pdata <- list()
pdata$genes <- c("a", "b", "c")
pdata$sigma <- rbind(c(1, 0.9, 0), c(0.9, 1, 0), c(0,0,1))
rownames(pdata$sigma) <- pdata$genes
colnames(pdata$sigma) <- pdata$genes

set.seed(0)
pdata$sim <- mvtnorm::rmvnorm(30, mean = c(10, 8, 3), sigma = pdata$sigma)
colnames(pdata$sim) <- pdata$genes

pdata$sim %>% cor()

pdata$sim <- pdata$sim %>%
  as_tibble() %>% 
  mutate(sample = formatC(1:nrow(pdata$sim), digits = 3, flag = "0"))

pdata$sim <- pdata$sim %>% gather(gene, expr_sim, -sample)

pdata$plot <- pdata$sim %>% 
  session$graphingUtils$ggplot(aes(x = sample, y = expr_sim)) + 
  geom_line(aes(group = gene), size = 1) + 
  theme(axis.text.x = element_blank()) + 
  theme(axis.text.y = element_blank()) +
  ylab("Expression") + 
  xlab("Sample")

ggsave(filename = paste0(workspace$workspaceDir, "figures_sandbox/", "schematic_01.png"), 
       plot = pdata$plot, device = "png", units = "cm", width = 20, height = 15) 









