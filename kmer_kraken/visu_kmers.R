library(ggplot2)
library(tidyverse)

data <- read.table("./kmers_kraken/result.txt", header = TRUE)

# proportion de reads alignés
ggplot(data, aes(x = condition, y = nb_aligned_kraken/4611205)) +
  geom_col(fill = "#8B978A") +
  geom_text(aes(label = round(nb_aligned_kraken/4611205, 3),
                vjust = -0.5)) +
  geom_hline(aes(yintercept = nb_aligned_bowtie/4611205, 
                 linetype = "Bowtie"),
             color = "darkred") +
  scale_x_continuous(labels=as.character(data$condition), breaks=data$condition) +
  scale_linetype_manual(values = 2) +
  labs(linetype = NULL, x = "Kraken2 k-mer length", y = "Proportion of aligned reads") +
  scale_y_continuous(limits=c(0, 1))


# concordance bowtie-kraken
ggplot(data, aes(x = condition, y = nb_concordant/nb_aligned_kraken)) +
  geom_col(fill = "#8A8B97") +
  geom_text(aes(label = round(nb_concordant/nb_aligned_kraken, 3),
                vjust = -0.5)) +
  scale_x_continuous(labels=as.character(data$condition), breaks=data$condition) +
  labs(linetype = NULL, x = "Kraken2 k-mer length", y = "Proportion of reads attributed\n to the same species") +
  scale_y_continuous(limits=c(0, 1))
