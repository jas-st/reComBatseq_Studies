library(dplyr)
library(viridis)
library(stringr)
library(tidyverse)

# remove all loaded variables except functions
rm(list = setdiff(ls(), lsf.str()))

#### HYPERPARAMETERS

temp = list.files(path="/Users/zhasmina/Desktop/EXPERIMENTS/HYPERPARAMETERS",
                  pattern="\\.csv$", full.names=TRUE)
myfiles = lapply(temp, read.csv)

full_df = bind_rows(myfiles)
full_df$alpha <- gsub("1e-04","0.0001",as.factor(full_df$alpha))
full_df$lambd <- gsub("1e-04","0.0001",as.factor(full_df$lambd))

max_points <- full_df %>% filter(prec > 0.85)
lmean = c("1"="No mean BE", "1.5"="Mean BE 1.5",
          "2"="Mean BE 2", "3"="Mean BE 3")
ldisp = c("1"="No dispersion BE", "2"="Dispersion BE 2",
          "3"="Dispersion BE 3", "4"="Dispersion BE 4")

ggplot(full_df, aes(x=alpha, y=lambd, fill=prec))+
  geom_tile(colour="white", linewidth=0.1) +
  geom_point(data = max_points, aes(x = alpha, y = lambd), shape = 8, color = "darkred", size = 2) +
  facet_grid(bmean~bdisp, labeller = labeller(bmean = lmean, bdisp = ldisp)) +
  labs(x="alpha", y="lambda") +
  scale_fill_viridis(option = "D", direction = 1, oob=scales::squish) +
  guides(fill = guide_colourbar(title = "Precision")) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave("/Users/zhasmina/Desktop/EXPERIMENTS/HYPERPARAMETERS/hyperparams_heatmap2.png", scale=2)


#### GENES

temp = list.files(path="/Users/zhasmina/Desktop/EXPERIMENTS/GENES",
                  pattern="results\\.csv$", recursive=T, full.names=TRUE)
myfiles = lapply(temp, read.csv)
identifier <- as.numeric(sapply(temp, str_extract,pattern="experiment_2\\((.*)\\)/",group=1))


#### SAMPLES

temp = list.files(path="/Users/zhasmina/Desktop/EXPERIMENTS/SAMPLES",
                  pattern="results\\.csv$", recursive=T, full.names=TRUE)
myfiles = lapply(temp, read.csv)
identifier <- as.numeric(sapply(temp, str_extract,pattern="experiment_2\\((.*)\\)/",group=1))

#### BATCHES

temp = list.files(path="/Users/zhasmina/Desktop/EXPERIMENTS/BATCHES_bigeffect",
                  pattern="results\\.csv$", recursive=T, full.names=TRUE)
myfiles = lapply(temp, read.csv)
identifier <- as.numeric(sapply(temp, str_extract,pattern="experiment_(.*)/",group=1))


## get full dataframe

full_df <- bind_rows(myfiles)
full_df["sample_length"] = rep(2^identifier,each=5) # gene_length otherwise
#full_df["batch_amount"] = rep(identifier,each=5)
full_df <- full_df %>%
  pivot_longer(cols = c(ComBat.seq, reComBat.seq, reComBat, pyComBat.seq, withBatch, withoutBatch),
  #pivot_longer(cols = c(reComBat.seq, reComBat, withBatch, withoutBatch),
               names_to = "Method")

# Plot time/prec/LDA

full_df_reduced <- full_df %>% filter((X %in% c("time", "prec", "lda")) &
                                     !(X == "prec" & Method == "withoutBatch") &
                                     !(X == "time" & Method %in% c("withoutBatch", "withBatch")) &
                                       (sample_length<500))

ggplot(full_df_reduced, aes(x = sample_length, y = value, color = Method)) + #x=gene_length,sample_length, batch_amount
  geom_line(show.legend = TRUE, linewidth = 1.2) +
  geom_point() +
  scale_color_brewer(palette="Set2") +
  facet_wrap(~X, scales = "free",
             labeller = labeller(X = c("lda" = "LDA","prec" = "Precision","time" = "Time (minutes)"))) +
  ylab(NULL) +
  xlab("Sample Amount") +  # Gene/Sample/Batch
  theme_minimal() +
  theme(strip.background = element_rect(fill = "grey90", color = NA),
        strip.text.x = element_text(size = 15), legend.position = "bottom")


ggsave("/Users/zhasmina/Desktop/EXPERIMENTS/SAMPLES/sample_plot_3.png",
       height=5, width=12, dpi=300)


# Plot tpr/fpr

full_df_reduced <- full_df %>% filter((X %in% c("tpr", "fpr") & Method != "withoutBatch"))

ggplot(full_df_reduced, aes(x = batch_amount, y = value, color = Method)) +
  geom_line(show.legend = TRUE, linewidth = 1.2) +
  geom_point() +
  scale_color_brewer(palette="Set2") +
  facet_wrap(~X, scales = "free",
             labeller = labeller(X = c("tpr" = "True Positive Rate", "fpr" = "False Positive Rate"))) +
  ylab(NULL) +
  xlab("Batch Amount") +
  theme_minimal() +
  theme(strip.background = element_rect(fill = "grey90", color = NA),
        strip.text.x = element_text(size = 15), legend.position = "bottom")


ggsave("/Users/zhasmina/Desktop/EXPERIMENTS/BATCHES_bigeffect/batch_plot_2.png",
       height=4, width=8, dpi=300)


#### COMPARISON

d <- DGEList(counts)
d <- edgeR::estimateCommonDisp(d)
plotMeanVar(d, show.raw.vars = TRUE, show.ave.raw.vars = FALSE, NBline = TRUE)
legend("topleft", legend = c("Poisson", "Negative Binomial"), fill = c("darkred", "darkblue"),
       border = c("darkred", "darkblue"), bty = "n", cex=0.7)
text(locator(), labels = c("Poisson", "Negative Binomial"))


