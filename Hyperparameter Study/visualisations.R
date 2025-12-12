library(dplyr)
library(viridis)
library(tidyverse)

#### HYPERPARAMETERS

## path to the csvs created with the pipeline (uncompressed)
temp = list.files(path="Hyperparameter Study",
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

#ggsave("Hyperparameter Study/hyperparams_heatmap2.png", scale=2)


