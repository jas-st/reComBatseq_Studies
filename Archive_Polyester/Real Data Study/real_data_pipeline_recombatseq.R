library(anndataR)
library(recombatseqv2)


h5ad_pfad <- "DATA/Adipose_tissue_exp.h5ad"

dat <- read_h5ad(h5ad_pfad)
tissue_name <- tools::file_path_sans_ext(basename(h5ad_pfad))
sce <- dat$as_SingleCellExperiment()

counts_df <- sce@assays@data@listData[["X"]]
colnames(counts_df) <- sce@colData@rownames

## Remove genes with near 0 counts
keep_lst_genes <- which(apply(counts_df, 1, function(x){sum(x)>ncol(counts_df)}))
rm_genes <- setdiff(1:nrow(counts_df), keep_lst_genes)
counts_df_reduced <- counts_df[keep_lst_genes, ]

cat("Gene Amount after Removal: ", nrow(counts_df_reduced), "\n")

batches <- sce[["batch"]]
group <- sce[["disease"]]

# which batches contain more than 1 sample
valid_batches = names(table(batches)[table(batches) > 1])
keep_lst_batches = which(batches %in% valid_batches)
rm_batches <- setdiff(1:ncol(counts_df_reduced), keep_lst_batches)
counts_df_reduced <- counts_df_reduced[,keep_lst_batches]

batches <- droplevels(batches[keep_lst_batches])
group <- droplevels(group[keep_lst_batches])

# batch correction
recombatseq_df <- ComBat_seq(counts_df_reduced, batch = batches, group = group,
                             lambda_reg=0.8, alpha_reg=0.3)


## Add back samples (so that dimensions won't change)
adjust_counts_whole <- matrix(NA, nrow=nrow(counts_df), ncol=ncol(counts_df))
dimnames(adjust_counts_whole) <- dimnames(counts_df)
adjust_counts_whole[keep_lst_genes, keep_lst_batches] <- recombatseq_df
adjust_counts_whole[, rm_batches] <- counts_df[, rm_batches]
adjust_counts_whole[rm_genes, ] <- counts_df[rm_genes, ]

write.csv(adjust_counts_whole, paste0("DATA/tissue_results/", tissue_name, "_corr.csv"))

