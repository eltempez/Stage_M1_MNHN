library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(gridExtra)

# Récupération des données
data <- read.table("./test_kraken_bowtie/results.txt", header = TRUE, dec = ".", sep = " ")
program <- c(rep("kr", 6), rep("bt_kr", 6))
database <- rep(c(rep("lu", 2), rep("standard", 2), rep("kefir", 2)), 2)
ech <- rep(c("mcf", "kA"), 6)
data_large <- cbind(ech, program, database, data[, -1])

############################
####### MAPPED READS #######
############################
# Tableau de répartition des reads
mapped_reads <- data_large[, 1:5]
p_unaligned <- 100 - (mapped_reads$p_aligned_reads_bowtie + mapped_reads$p_aligned_reads_kraken)
mapped_reads <- cbind(mapped_reads, p_unaligned)
# Transformation des données en format long
reads_long <- mapped_reads %>%
  pivot_longer(cols = -c(ech, program, database), names_to = "mapped", values_to = "percents")

# Visualisation
ggplot(reads_long, aes(x = database, y = percents, fill = factor(mapped, c("p_unaligned", "p_aligned_reads_kraken", "p_aligned_reads_bowtie")))) +
  geom_bar(stat = "identity") + 
  facet_grid(ech ~ program, labeller = label_both) +
  labs(x = "Library", y = "% of reads", fill = "Mapped") +
  scale_fill_manual(values=c('#868686', '#894c9f', '#4c799f'), 
                    name = "", labels = c("unmapped", "bowtie", "kraken"))





