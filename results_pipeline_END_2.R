library(polyester)
library(recombatseqv2)
library(ggplot2)
library(Rtsne)
library(scales)
library(DESeq2)
library(Biostrings)
library(MASS)
library(caret)
library(SingleCellExperiment)
library(scater)
library(ggpubr)
library(stringr)
library(dplyr)

# remove all loaded variables except functions
rm(list = setdiff(ls(), lsf.str()))

# gene lengths + experiment settings
gene_lengths <- t(read.csv("/Users/zhasmina/Downloads/gene_lengths_ch38.csv", row.names=1))
experiment <- "2(14)"
folder <- "SAMPLES"
N_total_sample <- round(2^12) # total number of samples

gene_count = round(2^10) # amount of genes
n_batches <- 2 # amount of batches
n_groups <- 2 # amount of biological groups
N_samples <- rep(N_total_sample/(n_batches*n_groups), n_groups*n_batches)
bio_fold <- 1.3 # biological signal
size_1 <- 1/0.15 # 1/dispersion batch1

#batch & biological vectors
batch <- rep(1:n_batches, each=N_total_sample/n_batches)
group <- rep(rep(0:(n_groups-1), n_batches), each=N_total_sample/(n_batches*n_groups))

batch_fold <- 1.5 # mean batch effect: mean(batch2) = 1.5*mean(batch1) 1,1.5,2,3
disp_fold_level <- 4 # dispersion batch effect: disp(batch2) = 5*disp(batch1) 1,2,3,4
size_2 <- 1/(0.15*disp_fold_level) # 1/dispersion batch2

lambda_reg <- 0.8
alpha_reg <- 0.3

res_list <- vector("list", 5)

for(iter in 1:5){
  cat(paste("Iteration", iter, "\n"))
  # set gene widths
  gene_widths <- sample(gene_lengths, gene_count)
  gene_names <- paste0("gene", 1:gene_count)

  ## true DE genes
  de_ground_truth_ind <- sample(1:gene_count, round(0.1*gene_count), replace=FALSE)
  G_ups <- de_ground_truth_ind[1:round(length(de_ground_truth_ind)/2)]
  G_downs <- de_ground_truth_ind[(round(length(de_ground_truth_ind)/2)+1):length(de_ground_truth_ind)]
  G_nulls <- setdiff(1:gene_count, c(G_ups, G_downs))
  #ground truth vectors
  de_ground_truth <- gene_names[de_ground_truth_ind]
  true_ups <- gene_names[G_ups]
  true_downs <- gene_names[G_downs]
  true_nulls <- gene_names[G_nulls]


  saveRDS(de_ground_truth_ind, file = paste0("/Users/zhasmina/Desktop/EXPERIMENTS/",folder,
                                   "/experiment_",experiment,
                                   "/iter",iter,"_DEgenes.rds"))

  ## Fold change matrix and size matrix
  #for baseline datasets without batch effect
  fold_changes_base <- constructFCMatrix_Comp(G=gene_count, FC_group=c(0,1,0,1), G_ups=G_ups,
                                              G_downs=G_downs, bioFC=bio_fold, batchFC=1)
  size_mat_base <- constructSizeMatrix(G=gene_count, size_vec=c(size_1, size_1, size_1, size_1))

  #for data with batch effect
  fold_changes <- constructFCMatrix_Comp(G=gene_count, FC_group=c(0,1,0,1), G_ups=G_ups,
                                         G_downs=G_downs, bioFC=bio_fold, batchFC=batch_fold)
  size_mat <- constructSizeMatrix(G=gene_count, size_vec=c(size_1, size_1, size_2, size_2))

  ## Generate count matrices
  seed <- sample.int(1e6, 1)
  nobatch_df <- simulate_experiment2(gene_widths, size=size_mat_base, num_reps=N_samples,
                                        fold_changes=fold_changes_base, seed=seed)
  batch_df <- simulate_experiment2(gene_widths, size=size_mat, num_reps=N_samples,
                                      fold_changes=fold_changes, seed=seed)

  rownames(batch_df) <- gene_names
  colnames(batch_df) <- paste0("Sample", seq_len(N_total_sample))
  rownames(nobatch_df) <- gene_names
  colnames(nobatch_df) <- paste0("Sample", seq_len(N_total_sample))

  write.csv(batch_df,
            paste0("/Users/zhasmina/Desktop/EXPERIMENTS/",folder,
                   "/experiment_",experiment,
                   "/iter",iter,"_batch_df.csv"))
  write.csv(nobatch_df,
            paste0("/Users/zhasmina/Desktop/EXPERIMENTS/",folder,
                   "/experiment_",experiment,
                   "/iter",iter,"_nobatch_df.csv"))

  count_batch_transformed <- DGEList(counts=batch_df)
  count_batch_transformed <- edgeR::calcNormFactors(count_batch_transformed, method="TMM")
  count_batch_transformed <- voom(count_batch_transformed, model.matrix(~as.factor(group)))

  write.csv(count_batch_transformed$E,
            paste0("/Users/zhasmina/Desktop/EXPERIMENTS/",folder,
                   "/experiment_",experiment,
                   "/iter",iter,"_countmat_batch_transformed.csv"))

  # Batch correction
  start.time <- Sys.time()
  combatseq_df <- ComBat_seq(batch_df, batch = as.factor(batch), group = as.factor(group))
  end.time <- Sys.time()
  time_noreg <- as.numeric(difftime(end.time,start.time, units="mins"))

  covmats <- createConfoundedDesign(N_total_sample)
  covmat <- covmats[[2]]
  #qr(covmat)$rank

  write.csv(covmats[[1]],
            paste0("/Users/zhasmina/Desktop/EXPERIMENTS/",folder,
                   "/experiment_",experiment,"/covmat.csv"))

  start.time <- Sys.time()
  recombatseq_df <- ComBat_seq(batch_df, batch = as.factor(batch), group = as.factor(group),
                                       covar_mod = covmat, lambda_reg=lambda_reg, alpha_reg=alpha_reg)
  end.time <- Sys.time()
  time_reg <- as.numeric(difftime(end.time,start.time, units="mins"))

  write.csv(combatseq_df,
            paste0("/Users/zhasmina/Desktop/EXPERIMENTS/",folder,
                   "/experiment_",experiment,
                   "/iter",iter,"_combatseq_df.csv"))
  write.csv(recombatseq_df,
            paste0("/Users/zhasmina/Desktop/EXPERIMENTS/",folder,
                   "/experiment_",experiment,
                   "/iter",iter,"_recombatseq_df.csv"))

  ## DE Analysis
  nobatch_cor <- edgeR_DEpipe(nobatch_df, batch=batch, group=group,
                              include.batch=FALSE, alpha.unadj=0.05, alpha.fdr=0.01)[["unadj"]]
  batch_cor <- edgeR_DEpipe(batch_df, batch=batch, group=group,
                            include.batch=FALSE, alpha.unadj=0.05, alpha.fdr=0.01)[["unadj"]]
  combatseq_cor <- edgeR_DEpipe(combatseq_df, batch=batch, group=group,
                              include.batch=FALSE, alpha.unadj=0.05, alpha.fdr=0.01)[["unadj"]]
  recombatseq_cor <- edgeR_DEpipe(recombatseq_df, batch=batch, group=group,
                                  include.batch=FALSE, alpha.unadj=0.05, alpha.fdr=0.01)[["unadj"]]

  stats_nobatch <- perfStats(nobatch_cor, de_ground_truth, gene_count)
  stats_batch <- perfStats(batch_cor, de_ground_truth, gene_count)
  stats_combatseq <- perfStats(combatseq_cor, de_ground_truth, gene_count)
  stats_recombatseq <- perfStats(recombatseq_cor, de_ground_truth, gene_count)

  # RESULTS DATA FRAME
  res_df <- data.frame(
    time=c(0, 0, time_noreg, time_reg, 0),
    tpr=c(stats_batch[["tpr"]], stats_nobatch[["tpr"]], stats_combatseq[["tpr"]], stats_recombatseq[["tpr"]], 0),
    fpr=c(stats_batch[["fpr"]], stats_nobatch[["fpr"]], stats_combatseq[["fpr"]], stats_recombatseq[["fpr"]], 0),
    prec=c(stats_batch[["prec"]], stats_nobatch[["prec"]], stats_combatseq[["prec"]], stats_recombatseq[["prec"]], 0),
    lda=c(0, 0, 0, 0, 0)
  )

  res_df <- t(res_df)
  colnames(res_df) <- c("withBatch", "withoutBatch", "ComBat-seq", "reComBat-seq", "reComBat")

  res_list[[iter]] <- res_df


  rm(batch_df, combatseq_df, recombatseq_df, nobatch_df, covmat,
     fold_changes, fold_changes_base, res_df,size_mat, size_mat_base,batch_cor,
     de_ground_truth, de_ground_truth_ind, end.time, G_downs, G_ups, G_nulls, gene_names,
     gene_widths, nobatch_cor, combatseq_cor, count_batch_transformed,
     recombatseq_cor, start.time,stats_batch, stats_nobatch, stats_combatseq, stats_recombatseq,
     time_noreg, time_reg, true_downs, true_nulls, true_ups)
}

res_df <- Reduce(`+`, res_list) / length(res_list)

#### save results
write.csv(res_df,
          paste0("/Users/zhasmina/Desktop/EXPERIMENTS/",folder,
                 "/experiment_",experiment,"/results.csv"))
