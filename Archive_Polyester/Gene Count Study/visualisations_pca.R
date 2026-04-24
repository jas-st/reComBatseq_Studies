source("helpers.R")
library(ggplot2)
library(ggpubr)
library(scater)

experiment <- "2(14)"
folder <- "Gene Count Study"
gene_count = round(2^14) # amount of genes

N_total_sample <- round(2^10) # total number of samples

n_batches <- 2 # amount of batches
n_groups <- 2 # amount of biological groups

batch <- rep(1:n_batches, each=N_total_sample/n_batches)
group <- rep(rep(0:(n_groups-1), n_batches), each=N_total_sample/(n_batches*n_groups))

for(iter in 1:5){
  cat(paste("Iteration", iter, "\n"))

  # read in count matrices
  recombat_df <- read.csv(paste0(folder, "/DATA/", experiment, "_iter", iter, "_recombat_df.csv"), row.names=1)
  batch_df <- read.csv(paste0(folder, "/DATA/", experiment, "_iter", iter, "_batch_df.csv"), row.names=1)
  nobatch_df <- read.csv(paste0(folder, "/DATA/", experiment, "_iter", iter, "_nobatch_df.csv"), row.names=1)
  combatseq_df <- read.csv(paste0(folder, "/DATA/", experiment, "_iter", iter, "_combatseq_df.csv"), row.names=1)
  recombatseq_df <- read.csv(paste0(folder, "/DATA/", experiment, "_iter", iter, "_recombatseq_df.csv"), row.names=1)
  pycombat_df <- read.csv(paste0(folder, "/DATA/", experiment, "_iter", iter, "_pycombat_df.csv"), row.names=1)

  ### PCA
  group_legend <- as_ggplot(get_legend(
    PCAplotter(PCAplotter_sce(batch_df, as.factor(batch), as.factor(group)), "Group") +
      guides(shape = guide_legend(override.aes = list(size = 10)),
             color = guide_legend(override.aes = list(size = 10))) +
      theme(legend.title = element_text(size = 20),
            legend.text  = element_text(size = 20))
  ))

  batch_legend <- as_ggplot(get_legend(
    PCAplotter(PCAplotter_sce(batch_df, as.factor(batch), as.factor(group)), "Batch") +
      guides(shape = guide_legend(override.aes = list(size = 10)),
             color = guide_legend(override.aes = list(size = 10))) +
      theme(legend.title = element_text(size = 20),
            legend.text  = element_text(size = 20))
  ))

  blank <- ggplot() + geom_blank() + theme_void()
  leg <- ggarrange(blank, group_legend, batch_legend, blank,
                   ncol=4, nrow=1, widths=c(1,2,2,1))

  batch_PCA <- PCAplotter_sce(batch_df, as.factor(batch), as.factor(group))
  nobatch_PCA <- PCAplotter_sce(nobatch_df, as.factor(batch), as.factor(group))
  combatseq_PCA <- PCAplotter_sce(combatseq_df, as.factor(batch), as.factor(group))
  recombatseq_PCA <- PCAplotter_sce(recombatseq_df, as.factor(batch), as.factor(group))
  pycombat_PCA <- PCAplotter_sce(pycombat_df, as.factor(batch), as.factor(group))

  recombat_sce <- SingleCellExperiment(assays=list(counts=recombat_df),
                                       colData=list(Batch=as.factor(batch), Group=as.factor(group)))
  recombat_sce@assays@data[["logcounts"]] <- scale(recombat_df)
  recombat_PCA <- runPCA(recombat_sce)

  plot_group1 <- PCAplotter(batch_PCA, "Group") + ggtitle("Unadjusted")
  plot_group2 <- PCAplotter(nobatch_PCA, "Group") + ggtitle("No batch effect")
  plot_group3 <- PCAplotter(combatseq_PCA, "Group") + ggtitle("ComBat-seq (non-singular design)")
  plot_group4 <- PCAplotter(recombatseq_PCA, "Group") + ggtitle("reComBat-seq (singular design)")
  plot_group5 <- PCAplotter(recombat_PCA, "Group") + ggtitle("reComBat (singular design)")
  plot_group6 <- PCAplotter(pycombat_PCA, "Group") + ggtitle("pyComBat-seq (singular design)")

  ggar1_1 <- ggarrange(plot_group1, plot_group2, plot_group3, nrow=1, ncol=3,
                       common.legend = TRUE, legend = "right")

  ggar1_2 <- ggarrange(plot_group4, plot_group5, plot_group6, nrow=1, ncol=3,
                       common.legend = TRUE, legend = "right")


  plot_batch1 <- PCAplotter(batch_PCA, "Batch") + ggtitle("")
  plot_batch2 <- PCAplotter(nobatch_PCA, "Batch") + ggtitle("")
  plot_batch3 <- PCAplotter(combatseq_PCA, "Batch") + ggtitle("")
  plot_batch4 <- PCAplotter(recombatseq_PCA, "Batch") + ggtitle("")
  plot_batch5 <- PCAplotter(recombat_PCA, "Batch") + ggtitle("")
  plot_batch6 <- PCAplotter(pycombat_PCA, "Batch") + ggtitle("")

  ggar2_1 <- ggarrange(plot_batch1, plot_batch2, plot_batch3, nrow=1, ncol=3,
                       common.legend = TRUE, legend = "right")
  ggar2_2 <- ggarrange(plot_batch4, plot_batch5, plot_batch6, nrow=1, ncol=3,
                       common.legend = TRUE, legend = "right")


  ggarrange(ggar1_1, ggar2_1, ggar1_2, ggar2_2, nrow=4,
            labels=c("A","A","B","B"))

  ggsave(paste0(folder, "/PCA Plots/", experiment, "_iter", iter, "_PCAplots.png"), height=10, width=10, dpi=300)


  rm(recombat_df, batch_df, nobatch_df, combatseq_df, recombatseq_df, group_legend,
     batch_legend, blank, leg, batch_PCA, nobatch_PCA, combatseq_PCA, recombatseq_PCA,
     recombat_sce, recombat_PCA, plot_group1, plot_group2, plot_group3, plot_group4,
     plot_group5, ggar2_1, ggar2_2, ggar1_1, ggar1_2, plot_batch1, plot_batch2, plot_batch3,
     plot_batch4, plot_batch5, plot_batch6, plot_group6, pycombat_df, pycombat_PCA)

}
