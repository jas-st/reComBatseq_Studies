#getJaccard <- function(de_genes, truegenes) {
#  intersection <- length(intersect(de_genes, truegenes))
#  union <- length(de_genes) + length(truegenes) - intersection
#  return(intersection/union)
#}

#PCAplotter_matrix <- function(countdf, batchvec, groupvec) {
#  coldat <- cbind(batchvec, groupvec)
#  colnames(coldat) <- c("batch", "group")

  #counts_norm <- as.matrix(apply(countdf, 2, function(x){x/sum(x)}))
#  counts_norm <- scale(normalizeCounts(countdf))

#  seobj <- SummarizedExperiment(assays=counts_norm, colData=coldat)
#  pca_obj <- plotPCA(DESeqTransform(seobj), intgroup=c("batch","group"))
#  pca_obj$data$group.1 <- as.factor(pca_obj$data$group.1)
#  pca_obj$data$batch <- as.factor(pca_obj$data$batch)
#  plt <- ggplot(pca_obj$data, aes(x=PC1, y=PC2, shape=group.1, colour=batch)) +
#    geom_point(size=3) + theme_bw() +
#    labs(x=sprintf("PC1: %s Variance", percent(pca_obj$plot_env$percentVar[1])),
#         y=sprintf("PC2: %s Variance", percent(pca_obj$plot_env$percentVar[2])))
#  plt
#}

#PCAplotter_group <- function(countdf, batchvec, groupvec) {
#  coldat <- cbind(batchvec, groupvec)
#  colnames(coldat) <- c("batch", "group")

  #counts_norm <- as.matrix(apply(countdf, 2, function(x){x/sum(x)}))
#  counts_norm <- scale(normalizeCounts(countdf))

#  seobj <- SummarizedExperiment(assays=counts_norm, colData=coldat)
#  pca_obj <- plotPCA(DESeqTransform(seobj), intgroup=c("batch","group"))
#  pca_obj$data$group.1 <- as.factor(pca_obj$data$group.1)
#  pca_obj$data$batch <- as.factor(pca_obj$data$batch)
#  plt <- ggplot(pca_obj$data, aes(x=PC1, y=PC2, shape=batch, colour=group.1)) +
#    geom_point(size=3) + theme_bw() +
#    labs(x=sprintf("PC1: %s Variance", percent(pca_obj$plot_env$percentVar[1])),
#         y=sprintf("PC2: %s Variance", percent(pca_obj$plot_env$percentVar[2])))
#  plt
#}

PCAplotter_sce <- function(countdf, batchvec, groupvec) {
  sce <- SingleCellExperiment(assays=list(counts=countdf),
                              colData=list(Batch=batchvec, Group=groupvec))
  #cnorm <- scale(countdf)
  #cnorm <- as.matrix(apply(countdf, 2, function(x){x/sum(x)}))
  cnorm <- scale(normalizeCounts(countdf))

  sce@assays@data[["logcounts"]] <- cnorm
  sim_pca <- runPCA(sce)
  return(sim_pca)
}

#LDAscore <- function(countdf, batchdf) {
#  X <- t(countdf)

#  filter1 <- colSums(X) > 1
#  X <- X[,filter1]
  #batchdf <- batchdf[filter1]

#  filter2 <- apply(X, 2, function(col) length(unique(col)))>2
#  X <- X[,filter2]
  #batchdf <- batchdf[filter2]

#  folds <- createFolds(batchdf, k = 10, list = TRUE)
#  scores <- numeric(10)

 # for (i in seq_along(folds)) {
#    val_idx <- folds[[i]]
#    train_idx <- setdiff(seq_along(batchdf), val_idx)

#    X_train <- X[train_idx, , drop = FALSE]
#    y_train <- batchdf[train_idx]

#    X_val <- X[val_idx, , drop = FALSE]
#    y_val <- batchdf[val_idx]

 #   lda_model <- lda(X_train, grouping = y_train)
#    preds <- predict(lda_model, X_val)$class

 #   acc <- mean(preds == y_val)
#    scores[i] <- acc
#  }
#  return(median(scores))
#}

PCAplotter <- function(sim_pca, colorby) {
  if (colorby=="Group") {
    plotPCA(sim_pca, colour_by = colorby) + scale_color_viridis_d(option='D') +
      geom_point(shape=21, stroke=0) + theme_bw() + labs(color=colorby)
  }
  else {
    plotPCA(sim_pca, colour_by = colorby) + scale_color_brewer(palette = "Set2") +
      geom_point(shape=21, stroke=0) + theme_bw() + labs(color=colorby)
  }
}

#Rtsneplotter_batch <- function(countlog, batchvec, groupvec) {
  # Calculate tSNE using Rtsne(0 function)
#  tsne_out <- Rtsne(t(countlog), check_duplicates = FALSE,
#                    perplexity = floor((ncol(countlog) - 1) / 3))

  # Conversion of matrix to dataframe
#  tsne_plot <- data.frame(x = tsne_out$Y[,1],
#                          y = tsne_out$Y[,2],
#                          batches = as.factor(batchvec),
#                          group = as.factor(groupvec))

  # Plotting the plot using ggplot() function
#  ggplot2::ggplot(tsne_plot) + geom_point(aes(x=x,y=y,colour=batches, shape=group), size=2) +
#    theme_bw()
#}

#Rtsneplotter_group <- function(countlog, batchvec, groupvec) {
  # Calculate tSNE using Rtsne(0 function)
#  tsne_out <- Rtsne(t(countlog), check_duplicates = FALSE,
 #                   perplexity = floor((ncol(countlog) - 1) / 3))

  # Conversion of matrix to dataframe
#  tsne_plot <- data.frame(x = tsne_out$Y[,1],
#                          y = tsne_out$Y[,2],
#                          batches = as.factor(batchvec),
#                          group = as.factor(groupvec))

  # Plotting the plot using ggplot() function
#  ggplot2::ggplot(tsne_plot) + geom_point(aes(x=x,y=y,colour=group, shape=batches), size=2) +
#    theme_bw()
#}

#trueDE <- function(object){
  #Output DE genes, defined by > 1.5 logFC, for the two target groups

#  output <- list()
#  rowdata <- rowData(object)
#  fc <- rowdata[,"DEFacGroup1"]/rowdata[,"DEFacGroup2"]
#  fc[fc < 1] <- 1/fc[fc < 1]
#  output$DEidx <- as.numeric(fc > 1.5)
#  output$DE <- as.character(rowdata$Gene)[output$DEidx==1]
#  return(output)
#}

edgeR_DEpipe <- function(counts_mat, batch, group, include.batch, alpha.unadj, alpha.fdr, covar=NULL){
  cat("DE tool: edgeR\n")
  y <- DGEList(counts=counts_mat)
  y <- edgeR::calcNormFactors(y, method="TMM")
  if(include.batch){
    cat("Including batch as covariate\n")
    design <- model.matrix(~ as.factor(group) + as.factor(batch))
  }else{
    cat("Default group as model matrix\n")
    design <- model.matrix(~as.factor(group))
  }
  if(!is.null(covar)){
    cat("Including surrogate variables or unwanted variation variables\n")
    design <- cbind(design, covar)
  }
  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmQLFit(y, design)
  qlf <- edgeR::glmQLFTest(fit, coef=2)
  de_res <- edgeR::topTags(qlf, n=nrow(counts_mat))$table

  de_called <- rownames(de_res)[de_res$PValue < alpha.unadj]
  de_called_fdr <- rownames(de_res)[de_res$FDR < alpha.fdr]
  return(list(unadj=de_called, fdr=de_called_fdr, de_res=de_res, design=design))
}

perfStats <- function(called_vec, ground_truth_vec, N_genes){
  if(length(called_vec)==0){
    tpr <- fpr <- 0
    prec <- NA
  }else{
    tp <- length(intersect(called_vec, ground_truth_vec))
    fp <- length(setdiff(called_vec, ground_truth_vec))
    N_DE <- length(ground_truth_vec)
    N_nonDE <- N_genes - N_DE

    tpr <- tp / N_DE
    fpr <- fp / N_nonDE
    prec <- tp / length(called_vec)
  }
  return(c(tpr=tpr, fpr=fpr, prec=prec))
}

createConfoundedDesign <- function(n_samples) {
  new_group <- replicate(3, sample(0:2, n_samples, replace=TRUE))
  last_col <- (new_group[,1]+1)%%3
  covmat <- data.frame()
  covmat <- cbind(new_group, last_col)
  covmat <- as.data.frame(covmat)
  covmat[] <- lapply(covmat, as.factor)
  covmat_model <- model.matrix(~., data = covmat)[, -1]

  return(list(covmat, covmat_model))
}
