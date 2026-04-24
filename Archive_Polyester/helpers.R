### Simulates RNA-seq read counts under a negative binomial model.
### Input:
### gene_width – transcript lengths or baseline expression levels
### num_reps – number of replicates per group (default: two groups of 10)
### fold_changes – multiplicative effect sizes between groups.
### size – dispersion parameter(s), numeric or matrix (default: basemeans / 3).
### Constructs baseline means per gene and group.Generates counts with the negative binomial (NB).
### Output: A matrix of simulated read counts (genes × samples).
simulate_experiment2 <- function (gene_width, num_reps = c(10, 10), size = NULL,
                                  fold_changes, paired = TRUE,...)
{
  extras = list(...)
  extras = .check_extras(extras, paired, total.n = sum(num_reps))
  if ("seed" %in% names(extras)) {
    set.seed(extras$seed)
  }

  #b0 = -3.0158
  #b1 = 0.8688
  #sigma = 4.152
  #logmus = b0 + b1 * log2(gene_width) + rnorm(length(gene_width), 0, sigma)
  #reads_per_transcript = 2^logmus - 1
  #reads_per_transcript = pmax(reads_per_transcript, 1e-06)
  reads_per_transcript = gene_width

  if (length(num_reps) == 1) {
    fold_changes = matrix(rep(1, length(gene_width)))
    basemeans = reads_per_transcript * fold_changes
  }
  else if (length(num_reps) == 2) {
    if (length(reads_per_transcript) == 1) {
      basemeans = matrix(reads_per_transcript, ncol = 2,
                         nrow = length(gene_width))
      basemeans[, 2] = fold_changes * basemeans[, 1]
    }
    else if (length(reads_per_transcript) == length(gene_width)) {
      basemeans = matrix(c(reads_per_transcript, reads_per_transcript),
                         nrow = length(reads_per_transcript))
      basemeans = basemeans * fold_changes
    }
    else {
      stop("reads_per_transcript is the wrong length.")
    }
  }
  else {
    basemeans = reads_per_transcript * fold_changes
  }
  if (is.null(size)) {
    size = basemeans/3
  }
  else if (is.numeric(size)) {
    if (is.matrix(basemeans)) {
      num_rows = nrow(basemeans)
      num_cols = ncol(basemeans)
    }
    else {
      num_rows = length(basemeans)
      num_cols = 1
    }
    size = matrix(size, nrow = num_rows, ncol = num_cols)
  }
  else if (class(size) == "matrix") {
    if (!is.matrix(basemeans)) {
      stop("If you provide a matrix for size, you also need a matrix for reads_per_transcript.")
    }
    stopifnot(nrow(size) == nrow(basemeans))
    stopifnot(ncol(size) == ncol(basemeans))
  }
  else {
    stop("size must be a number, numeric vector, or matrix.")
  }
  if ("seed" %in% names(extras)) {
    set.seed(extras$seed)
  }
  group_ids = rep(1:length(num_reps), times = num_reps)
  numreadsList = vector("list", sum(num_reps))
  numreadsList = lapply(1:sum(num_reps), function(i) {
    group_id = group_ids[i]
    NB(as.matrix(basemeans)[, group_id], as.matrix(size)[,
                                                         group_id])
  })
  readmat = matrix(unlist(numreadsList), ncol = sum(num_reps))
  readmat = t(extras$lib_sizes * t(readmat))

  return(readmat)
}

### Builds a fold-change matrix that encodes both biological and batch effects across genes and sample groups.
### Output: A gene × group matrix of fold changes combining batch and biological effects.
constructFCMatrix_Comp <- function(G, FC_group, G_ups, G_downs, bioFC, batchFC){
  fold_changes <- matrix(1, nrow=G, ncol=length(FC_group))

  # randomly split all genes into 2 groups - one increased in batch 2 vs 1, the other decreased
  # samples 1 and 2 918 times for the categorization of the genes
  # in the batch free case - the biological signal
  batch_genes_split <- sample(c(1,2), G, replace=TRUE)
  batch_up_genes <- which(batch_genes_split==1)
  batch_down_genes <- which(batch_genes_split==2)


  # batch fold changes
  # in the case of biological signal it doesnt matter its all 1,1,1,1 matrix
  # in the case of batch signal - batch up are the genes with higher signal in batch 2 etc
  fold_changes[batch_up_genes, ] <- matrix(rep(c(1, 1, batchFC, batchFC), length(batch_up_genes)),
                                           nrow=length(batch_up_genes), byrow=TRUE)
  fold_changes[batch_down_genes, ] <- matrix(rep(c(batchFC, batchFC, 1, 1), length(batch_down_genes)),
                                             nrow=length(batch_down_genes), byrow=TRUE)


  # add biological fold changes
  fold_changes[G_ups, FC_group==1] <- fold_changes[G_ups, FC_group==1] * bioFC
  fold_changes[G_downs, FC_group==0] <- fold_changes[G_downs, FC_group==0] * bioFC

  return(fold_changes)
}

### Constructs a dispersion (size) matrix for negative binomial simulations.
### Output: A gene × group matrix of size (dispersion) values.
constructSizeMatrix <- function(G, size_vec){
  size_mat <- matrix(rep(size_vec, G), ncol=length(size_vec), byrow=TRUE)
  return(size_mat)
}


# function to make sure everything is compatible and assign sane defaults
# internal

.check_extras = function(extras, paired, total.n){

  if(!('readlen' %in% names(extras))){
    extras$readlen = 100
  }

  if(!('error_model' %in% names(extras))){
    extras$error_model = 'uniform'
  }

  if(!('error_rate' %in% names(extras))){
    extras$error_rate = 0.005
  }
  if(extras$error_model == 'custom'){
    extras$path = paste0(extras$model_path, '/', extras$model_prefix)
  }#this should work beause we already checked stuff.

  if(!('bias' %in% names(extras))){
    extras$bias = 'none'
  }else{
    extras$bias = match.arg(extras$bias, c('none', 'rnaf', 'cdnaf'))
  }

  if(!('lib_sizes' %in% names(extras))){
    extras$lib_sizes = rep(1, total.n)
  }else{
    stopifnot(is.numeric(extras$lib_sizes))
    stopifnot(length(extras$lib_sizes) == total.n)
  }

  return(extras)

}

### Creates a SingleCellExperiment object from count data, normalizes it, runs PCA, and returns the PCA-transformed object
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

### Generates a PCA plot, coloring points by the given variable
### the two functions were separated for more flexibility and to not repeat the pca every time we want a different color
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

### edgeR-based differential expression pipeline, taken from ComBat-seq
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

### Computes performance metrics (TPR, FPR, precision) for detected DE genes against a ground truth set
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

### Generates a simulated experimental design matrix with confounded variables
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
