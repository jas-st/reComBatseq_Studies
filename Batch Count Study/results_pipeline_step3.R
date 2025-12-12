experiment <- "2"
folder <- "Batch Count Study"

gene_count = round(2^10) # amount of genes
N_total_sample <- 6000 # total number of samples

stats_recombat <- vector("list", 5)

for(iter in 1:5){
  cat(paste("Iteration", iter, "\n"))

  recombat_df <- read.csv(paste0(folder, "/DATA/", experiment, "_iter", iter, "_recombat_df.csv"), row.names=1)

  metadata <- read.csv(paste0(folder, "/DATA/", experiment, "_iter", iter, "_metadata.csv"), row.names=1)
  group <- as.factor(metadata[["group"]])


  # Differential expression
  de_genes <- readRDS(paste0(folder, "/DATA/", experiment, "_iter", iter, "_DEgenes.rds"))
  de_genes <- as.numeric(unlist(str_extract_all(de_genes, "\\d+")))

  # reComBat
  vfit <- lmFit(recombat_df, model.matrix(~as.factor(group)))
  efit <- eBayes(vfit)
  tests <- decideTests(efit)
  recombat_cor <- which(tests[,2]!=0)

  # Stats
  stats_recombat[[iter]] <- perfStats(recombat_cor, de_genes, gene_count)

  rm(de_genes, vfit, efit, tests, recombat_cor, pycombat_cor)
}

colMeans(do.call(rbind, stats_recombat))
