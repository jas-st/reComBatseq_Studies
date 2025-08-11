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

# gene lengths + experiment settings
gene_lengths <- t(read.csv("/Users/zhasmina/Downloads/gene_lengths_ch38.csv", row.names=1))
experiment <- "2(7)"
folder <- "HYPERPARAMETERS"

mean_batchfolds <- c(1,1.5,2,3)
disp_batchfolds <- c(1,2,3,4)
alpha_regfolgs <- c(0,0.0001,0.001,0.1,0.3,0.5,0.7,0.9)
lambda_regfolds <- c(0.0001,0.001,0.1,0.3,0.5,0.7,0.9,1)

N_total_sample <- round(2^7) # total number of samples
n_batches <- 2 # amount of batches
n_groups <- 2 # amount of biological groups
N_samples <- rep(N_total_sample/(n_batches*n_groups), n_groups*n_batches)
gene_count = round(2^7) # amount of genes
bio_fold <- 1.3 # biological signal
size_1 <- 1/0.15 # 1/dispersion batch1

#batch & biological vectors
batch <- rep(1:n_batches, each=N_total_sample/n_batches)
group <- rep(rep(0:(n_groups-1), n_batches), each=N_total_sample/(n_batches*n_groups))

i_no = 1
experiment_list <- vector("list", length(alpha_regfolgs)*length(lambda_regfolds))

batch_fold <- 3 # mean batch effect: mean(batch2) = 1.5*mean(batch1) 1,1.5,2,3
disp_fold_level <- 4 # dispersion batch effect: disp(batch2) = 5*disp(batch1) 1,2,3,4

size_2 <- 1/(0.15*disp_fold_level) # 1/dispersion batch2

for (areg in alpha_regfolgs) {
  for (lreg in lambda_regfolds) {

    lambda_reg=lreg
    alpha_reg=areg

    res_list <- vector("list", 5)

    for(iter in 1:5){
      cat(paste("Iteration", iter, lambda_reg,alpha_reg,"\n"))
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
      count_nobatch <- simulate_experiment2(gene_widths, size=size_mat_base, num_reps=N_samples,
                                            fold_changes=fold_changes_base, seed=seed)
      count_batch <- simulate_experiment2(gene_widths, size=size_mat, num_reps=N_samples,
                                          fold_changes=fold_changes, seed=seed)

      rownames(count_batch) <- gene_names
      colnames(count_batch) <- paste0("Sample", seq_len(N_total_sample))
      rownames(count_nobatch) <- gene_names
      colnames(count_nobatch) <- paste0("Sample", seq_len(N_total_sample))

      # Batch correction
      start.time <- Sys.time()
      count_batchremoved <- ComBat_seq(count_batch, batch = as.factor(batch), group = as.factor(group))
      end.time <- Sys.time()
      time_noreg <- as.numeric(difftime(end.time,start.time, units="mins"))

      covmat <- createConfoundedDesign(N_total_sample)
      qr(covmat)$rank

      start.time <- Sys.time()
      count_batchremoved_reg <- ComBat_seq(count_batch, batch = as.factor(batch), group = as.factor(group),
                                           covar_mod = covmat, lambda_reg=lambda_reg, alpha_reg=alpha_reg)
      end.time <- Sys.time()
      time_reg <- as.numeric(difftime(end.time,start.time, units="mins"))

      ## DE Analysis
      nobatch_cor <- edgeR_DEpipe(count_nobatch, batch=batch, group=group,
                                  include.batch=FALSE, alpha.unadj=0.05, alpha.fdr=0.01)[["unadj"]]
      batch_cor <- edgeR_DEpipe(count_batch, batch=batch, group=group,
                                include.batch=FALSE, alpha.unadj=0.05, alpha.fdr=0.01)[["unadj"]]
      removed_cor <- edgeR_DEpipe(count_batchremoved, batch=batch, group=group,
                                  include.batch=FALSE, alpha.unadj=0.05, alpha.fdr=0.01)[["unadj"]]
      removed_reg_cor <- edgeR_DEpipe(count_batchremoved_reg, batch=batch, group=group,
                                      include.batch=FALSE, alpha.unadj=0.05, alpha.fdr=0.01)[["unadj"]]

      stats_nobatch <- perfStats(nobatch_cor, de_ground_truth, gene_count)
      stats_batch <- perfStats(batch_cor, de_ground_truth, gene_count)
      stats_removed <- perfStats(removed_cor, de_ground_truth, gene_count)
      stats_removed_reg <- perfStats(removed_reg_cor, de_ground_truth, gene_count)

      res_df <- data.frame(
        time=c(0,0,time_noreg, time_reg),
        tpr=c(stats_batch[["tpr"]], stats_nobatch[["tpr"]], stats_removed[["tpr"]], stats_removed_reg[["tpr"]]),
        fpr=c(stats_batch[["fpr"]], stats_nobatch[["fpr"]], stats_removed[["fpr"]], stats_removed_reg[["fpr"]]),
        prec=c(stats_batch[["prec"]], stats_nobatch[["prec"]], stats_removed[["prec"]], stats_removed_reg[["prec"]])
      )
      res_df <- t(res_df)
      colnames(res_df) <- c("withBatch", "withoutBatch", "notRegularized", "Regularized")

      res_list[[iter]] <- res_df


      rm(count_batch, count_batchremoved, count_batchremoved_reg, count_nobatch, covmat,
         fold_changes, fold_changes_base, res_df,size_mat, size_mat_base,batch_cor,
         de_ground_truth, de_ground_truth_ind, end.time, G_downs, G_ups, gene_names,
         gene_widths, nobatch_cor, removed_cor,
         removed_reg_cor, start.time,stats_batch, stats_nobatch, stats_removed, stats_removed_reg,
         time_noreg, time_reg, true_downs, true_nulls, true_ups)
    }
    mean_precision = mean(sapply(res_list, function(a) a["prec","Regularized"]))
    df_row <- c(batch_fold, disp_fold_level, alpha_reg, lambda_reg, mean_precision)
    experiment_list[[i_no]] <- df_row
    i_no <- i_no + 1
  }
}

for (fix in 57:64) {
  experiment_list[[fix]] <- c(3, 4, 0.9, lambda_regfolds[fix-56], 0.1)
}

hyperparams_df <- t(as.data.frame(experiment_list))
colnames(hyperparams_df) <- c("bmean","bdisp","alpha","lambd","prec")
rownames(hyperparams_df) <- seq_len(nrow(hyperparams_df))
hyperparams_df <- as.data.frame(hyperparams_df)
hyperparams_df$alpha <- gsub("1e-04","0.0001",as.factor(hyperparams_df$alpha))
hyperparams_df$lambd <- gsub("1e-04","0.0001",as.factor(hyperparams_df$lambd))


write.csv(hyperparams_df,
          paste0("/Users/zhasmina/Desktop/EXPERIMENTS/",folder,"/mbfold",
          as.character(batch_fold),"_dbfold",as.character(disp_fold_level),".csv"))


ggplot(hyperparams_df, aes(x=alpha, y=lambd, fill=prec))+
  geom_tile() +
  facet_grid(bmean~bdisp)


## PCA
plot_group1 <- PCAplotter_sce(count_batch, as.factor(batch), as.factor(group), "Group")
plot_batch1 <- PCAplotter_sce(count_batch, as.factor(batch), as.factor(group), "Batch")

plot_group2 <- PCAplotter_sce(count_nobatch, as.factor(batch), as.factor(group), "Group")
plot_batch2 <- PCAplotter_sce(count_nobatch, as.factor(batch), as.factor(group), "Batch")

ggarrange(plot_group1, plot_batch1, plot_group2, plot_batch2, ncol=2, nrow=2,
          labels=c("A", "", "B"))
ggsave(paste0("/Users/zhasmina/Desktop/EXPERIMENTS/",folder,
              "/experiment_",experiment,
              "/bfold",as.character(batch_fold),
              "_dfold",as.character(disp_fold_level),"_batchnobatch.png"))

plot_group1 <- PCAplotter_sce(count_batchremoved, as.factor(batch), as.factor(group), "Group")
plot_batch1 <- PCAplotter_sce(count_batchremoved, as.factor(batch), as.factor(group), "Batch")

plot_group2 <- PCAplotter_sce(count_batchremoved_reg, as.factor(batch), as.factor(group), "Group")
plot_batch2 <- PCAplotter_sce(count_batchremoved_reg, as.factor(batch), as.factor(group), "Batch")

ggarrange(plot_group1, plot_batch1, plot_group2, plot_batch2, ncol=2, nrow=2,
          labels=c("A", "", "B"))
ggsave(paste0("/Users/zhasmina/Desktop/EXPERIMENTS/",folder,
              "/experiment_",experiment,
              "/bfold",as.character(batch_fold),
              "_dfold",as.character(disp_fold_level),"_regnoreg.png"))

#### save results
write.csv(res_df,
          paste0("/Users/zhasmina/Desktop/EXPERIMENTS/",folder,
                 #"/experiment_",experiment,
                 "/",as.character(format(lambda_reg, scientific=FALSE)),
                 "_",as.character(format(alpha_reg, scientific=FALSE)),
                 "/bfold",as.character(batch_fold),
                 "_dfold",as.character(disp_fold_level),"_res_df.csv"))
