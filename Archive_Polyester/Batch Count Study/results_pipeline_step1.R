library(airway)
library(recombatseqv2)
library(splatter)
library(stringr)
source("../helpers.R")

data(airway)
airway_counts <- as.matrix(assay(airway))
params <- splatEstimate(airway_counts)
params <- newSplatParams()
folder <- "Batch Count Study"

n_batches = 40
experiment <- "40"

gene_count = round(2^10)
sample_count = 6000
lambda_reg <- 0.8
alpha_reg <- 0.3

### SIMULATE
res_list <- vector("list", 5)

for(iter in 1:5){
  cat(paste("Iteration", iter, "\n"))
  seed <- sample.int(1e6, 1)
  facLoc = abs(rnorm(n_batches, mean=2.5, sd=2.5))
  facScale = abs(rnorm(n_batches, mean=2.5, sd=2.5))

  sim_batch <- splatSimulate(nGenes=gene_count, batchCells=rep(sample_count/n_batches, n_batches),
                       mean.rate = params@mean.rate, mean.shape = params@mean.shape,
                       bcv.common = params@bcv.common, bcv.df = params@bcv.df,
                       group.prob=c(0.5,0.5), de.prob=0.1,
                       batch.facLoc=2, batch.facScale=0.8,
                       method='groups',
                       verbose=F,
                       batch.rmEffect=FALSE, seed=seed)


  sim_noBatch <- splatSimulate(nGenes=gene_count, batchCells=rep(sample_count/n_batches,n_batches),
                               mean.rate = params@mean.rate, mean.shape = params@mean.shape,
                               bcv.common = params@bcv.common, bcv.df = params@bcv.df,
                               group.prob=c(0.5,0.5), de.prob=0.1,
                               batch.facLoc=2, batch.facScale=0.8,
                               method='groups',
                               verbose=F,
                               batch.rmEffect=TRUE, seed=seed)

  batch <- colData(sim_batch)@listData[["Batch"]]
  group <- colData(sim_batch)@listData[["Group"]]
  batch <- as.factor(as.numeric(unlist(str_extract_all(batch, "\\d+"))))
  group <- as.factor(as.numeric(unlist(str_extract_all(group, "\\d+"))))

  metadata <- cbind(batch, group)
  write.csv(metadata,
            paste0(folder, "/DATA/", experiment, "/iter",iter,"_metadata.csv"))

  batch_df <- as.matrix(counts(sim_batch))
  nobatch_df <- as.matrix(counts(sim_noBatch))

  write.csv(batch_df,
            paste0(folder, "/DATA/", experiment, "/iter",iter,"_batch_df.csv"))
  write.csv(nobatch_df,
            paste0(folder, "/DATA/", experiment, "/iter",iter,"_nobatch_df.csv"))

  count_batch_transformed <- DGEList(counts=batch_df)
  count_batch_transformed <- edgeR::calcNormFactors(count_batch_transformed, method="TMM")
  count_batch_transformed <- voom(count_batch_transformed, model.matrix(~as.factor(group)))

  write.csv(count_batch_transformed$E,
            paste0(folder, "/DATA/", experiment, "/iter",iter,"_countmat_batch_transformed.csv"))

  # Batch correction
  covmats <- createConfoundedDesign(sample_count)
  covmat <- covmats[[2]]
  #qr(covmat)$rank

  write.csv(covmats[[1]],
            paste0(folder, "/DATA/", experiment, "/covmat.csv"))

  start.time <- Sys.time()
  recombatseq_df <- reComBat_seq(batch_df, batch = as.factor(batch), group = as.factor(group),
                               covar_mod = covmat, lambda_reg=lambda_reg, alpha_reg=alpha_reg)
  end.time <- Sys.time()
  time_reg <- as.numeric(difftime(end.time,start.time, units="mins"))

  write.csv(recombatseq_df,
            paste0(folder, "/DATA/", experiment, "/iter",iter,"_recombatseq_df.csv"))
  
  ## DE Analysis
  de_cols <- grep("^DEFacGroup", colnames(rowData(sim_batch)), value = TRUE)
  de_sums <- rowSums(as.matrix(rowData(sim_batch)[, de_cols]))
  de_ground_truth <- rownames(rowData(sim_batch))[de_sums != length(de_cols)]

  saveRDS(de_ground_truth, file = paste0(folder, "/DATA/", experiment, "/iter",iter,"_DEgenes.rds"))

  nobatch_cor <- edgeR_DEpipe(nobatch_df, batch=batch, group=group,
                              include.batch=FALSE, alpha.unadj=0.05, alpha.fdr=0.01)[["fdr"]]
  batch_cor <- edgeR_DEpipe(batch_df, batch=batch, group=group,
                            include.batch=FALSE, alpha.unadj=0.05, alpha.fdr=0.01)[["fdr"]]
  recombatseq_cor <- edgeR_DEpipe(recombatseq_df, batch=batch, group=group,
                                  include.batch=FALSE, alpha.unadj=0.05, alpha.fdr=0.01)[["fdr"]]

  stats_nobatch <- perfStats(nobatch_cor, de_ground_truth, gene_count)
  stats_batch <- perfStats(batch_cor, de_ground_truth, gene_count)
  stats_recombatseq <- perfStats(recombatseq_cor, de_ground_truth, gene_count)

  # RESULTS DATA FRAME
  res_df <- data.frame(
    time=c(0, 0, time_reg, 0),
    tpr=c(stats_batch[["tpr"]], stats_nobatch[["tpr"]], stats_recombatseq[["tpr"]], 0),
    fpr=c(stats_batch[["fpr"]], stats_nobatch[["fpr"]], stats_recombatseq[["fpr"]], 0),
    prec=c(stats_batch[["prec"]], stats_nobatch[["prec"]], stats_recombatseq[["prec"]], 0),
    lda=c(0, 0, 0, 0)
  )

  res_df <- t(res_df)
  colnames(res_df) <- c("withBatch", "withoutBatch", "reComBat-seq", "reComBat")

  res_list[[iter]] <- res_df


  rm(sim_batch, sim_noBatch, batch, group, metadata, batch_df, recombatseq_df,
     nobatch_df, covmat, covmats, res_df, batch_cor, de_ground_truth, end.time,
     nobatch_cor, count_batch_transformed, recombatseq_cor, start.time, stats_batch,
     stats_nobatch, stats_recombatseq, time_reg)
}

res_df <- Reduce(`+`, res_list) / length(res_list)

#### save results
write.csv(res_df,
          paste0(folder, "/DATA/", experiment, "/results.csv"))
