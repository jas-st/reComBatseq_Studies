library(dplyr)
library(viridis)
library(stringr)
library(tidyverse)

#### BATCHES

# path containing the results dataframes
temp = list.files(path="/Batch Count Study/DATA",
                  pattern="results\\.csv$", recursive=T, full.names=TRUE)
myfiles = lapply(temp, read.csv)
identifier <- as.numeric(sapply(temp, str_extract,pattern="experiment_(.*)/",group=1))

## get full dataframe

full_df <- bind_rows(myfiles)
full_df["batch_amount"] = rep(identifier,each=5)
full_df <- full_df %>%
  pivot_longer(cols = c(reComBat.seq, reComBat, withBatch, withoutBatch),
               names_to = "Method")

# Plot time/prec/LDA

full_df_reduced <- full_df %>% filter((X %in% c("time", "prec", "lda")) &
                                     !(X == "prec" & Method == "withoutBatch") &
                                     !(X == "time" & Method %in% c("withoutBatch", "withBatch")) &
                                       (sample_length<500))

ggplot(full_df_reduced, aes(x = batch_amount, y = value, color = Method)) +
  geom_line(show.legend = TRUE, linewidth = 1.2) +
  geom_point() +
  scale_color_brewer(palette="Set2") +
  facet_wrap(~X, scales = "free",
             labeller = labeller(X = c("lda" = "LDA","prec" = "Precision","time" = "Time (minutes)"))) +
  ylab(NULL) +
  xlab("Batch Amount") + 
  theme_minimal() +
  theme(strip.background = element_rect(fill = "grey90", color = NA),
        strip.text.x = element_text(size = 15), legend.position = "bottom")


ggsave("batch_plot.png", height=5, width=12, dpi=300)


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


ggsave("batch_plot_2.png",height=4, width=8, dpi=300)


