library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(gridExtra)
library(ggpattern)

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


########################################
####### LACTO, ACETO, AND YEASTS #######
########################################
# Tableau de distribution
distrib <- data_large[, c(1:3, 6:11)]
p_undistrib_kraken <- 100 - (distrib$p_lacto_kraken + distrib$p_aceto_kraken + distrib$p_yeasts_kraken)
p_undistrib_bowtie <- 100 - (distrib$p_lacto_bowtie + distrib$p_aceto_bowtie + distrib$p_yeasts_bowtie)
distrib <- cbind(distrib, p_undistrib_bowtie, p_undistrib_kraken)
# Transformation des pourcentages relatifs en pourcentages globaux
equiv_bowtie <- data_large$p_aligned_reads_bowtie / (data_large$p_aligned_reads_bowtie + data_large$p_aligned_reads_kraken) * 100
equiv_kraken <- data_large$p_aligned_reads_kraken / (data_large$p_aligned_reads_bowtie + data_large$p_aligned_reads_kraken) * 100
distrib[, c(4, 6, 8, 10)] <- distrib[, c(4, 6, 8, 10)] * equiv_bowtie/100
distrib[, c(5, 7, 9, 11)] <- distrib[, c(5, 7, 9, 11)] * equiv_kraken/100
# Transformation des données en format long
distrib_long <- distrib %>%
  pivot_longer(cols = -c(ech, program, database), 
               names_to = c(".value", "from"),
               names_pattern = "(.+)_(.+)") %>%
  pivot_longer(cols = -c(ech, program, database, from), 
               names_to = "classif", values_to = "percents")


# Visualisation
ggplot(distrib_long, aes(x = database, y = percents, 
                         fill = factor(classif, c("p_undistrib", "p_yeasts",
                                                  "p_aceto", "p_lacto")),
                         pattern = from)) +
  geom_col() +
  scale_fill_manual(values=c('#9B4E89', '#BC843B', '#3F5198', '#3F984E'), 
                    name = "type of micro-organism", 
                    labels = c("other", "yeasts", "acetobacter", "lactobacillus"),
                    guide = guide_legend(override.aes = list(pattern = 'none'))) +
  geom_col_pattern(pattern_density = 0.1, 
                   pattern_colour = rgb(230/255, 230/255, 230/255, alpha = 0.2),
                   pattern_fill = rgb(230/255, 230/255, 230/255, alpha = 0.2),
                   pattern_spacing = 0.02,
                   width = 0.88,
                   ) +
  facet_grid(ech ~ program, labeller = label_both) +
  labs(x = "Library", y = "% of reads", fill = "Type of micro-organism") +
  scale_pattern_manual(name = "obtained from", 
                       values = c('none', 'crosshatch')) 



#################################
####### NUMBER OF SPECIES #######
#################################
ggplot(data_large, aes(x = database, y = nb_species_kraken, group = program)) +
  geom_col(aes(fill = program), position = "dodge") +
  geom_text(aes(label = nb_species_kraken),
            position = position_dodge(0.9),
            vjust = -0.5) +
  facet_grid(ech ~ ., labeller = label_both) +
  labs(x = "Library", y = "number of reads", fill = "Mapped") +
  scale_fill_manual(values=c('#8873B2', '#73B281'), 
                    name = "Programs", labels = c("bowtie + kraken", "kraken")) +
  scale_y_continuous(limits=c(0,8000))
