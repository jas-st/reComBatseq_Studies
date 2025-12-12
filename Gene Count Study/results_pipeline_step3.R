library(edgeR)
library(limma)
source("../helpers.R")

experiment <- "2(14)"
folder <- "Gene Count Study"
gene_count = round(2^14) # amount of genes

N_total_sample <- round(2^10) # total number of samples
n_batches <- 2 # amount of batches
n_groups <- 2 # amount of biological groups

batch <- rep(1:n_batches, each=N_total_sample/n_batches)
group <- rep(rep(0:(n_groups-1), n_batches), each=N_total_sample/(n_batches*n_groups))

stats_recombat <- vector("list", 5)
stats_pycombat <- vector("list", 5)

for(iter in 1:5){
  cat(paste("Iteration", iter, "\n"))

  pycombat_df <- read.csv(paste0(folder, "/DATA/", experiment, "_iter", iter, "_pycombat_df.csv"), row.names=1)
  recombat_df <- read.csv(paste0(folder, "/DATA/", experiment, "_iter", iter, "_recombat_df.csv"), row.names=1)


  # Differential expression
  de_genes <- readRDS(paste0(folder, "/DATA/", experiment, "_iter", iter, "_DEgenes.rds"))

  # reComBat
  vfit <- lmFit(scale(recombat_df), model.matrix(~as.factor(group)))
  efit <- eBayes(vfit)
  tests <- decideTests(efit)
  recombat_cor <- which(tests[,2]!=0)

  # pyComBat
  pycombat_cor <- edgeR_DEpipe(pycombat_df, batch=batch, group=group,
                                  include.batch=FALSE, alpha.unadj=0.05, alpha.fdr=0.01)[["unadj"]]
  pycombat_cor <- as.numeric(unlist(stringr::str_extract_all(pycombat_cor, "\\d+")))

  # Stats
  stats_recombat[[iter]] <- perfStats(recombat_cor, de_genes, gene_count)
  stats_pycombat[[iter]] <- perfStats(pycombat_cor, de_genes, gene_count)


  rm(de_genes, vfit, efit, tests, recombat_cor, pycombat_cor)
}

colMeans(do.call(rbind, stats_recombat))
colMeans(do.call(rbind, stats_pycombat))
