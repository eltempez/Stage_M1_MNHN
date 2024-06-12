## Révision des graphiques, simplifiés pour la soutenance
# et ré-établis à partir des génomes mcf sans B15 et B16
library(eulerr)
library(ggplot2)

###############################
####### NOMBRE DE READS #######
###############################
data <- read.table("./test_kraken_bowtie/soutenance/results.txt", header = TRUE, sep = "\t")

ggplot(data, aes(x = kefir, y = p_mapped, fill = prg)) +
  geom_col(position = "dodge") +
  labs(x = "", y = "% de lectures alignées", fill = "Programme") +
  theme_minimal()




#########################
####### MISMATCHS #######
#########################
venn_simple <- function(file_mismatch, kefir) {
  # lecture des données
  mis_data <- read.table(file_mismatch, sep = "\t", header = TRUE, skip = 1)[, -1]
  nb_total_reads <- read.table(file_mismatch, nrows = 1, sep = "\t")[, 2]
  
  # distribution
  nb_matched_reads <- nb_total_reads - nrow(mis_data)
  nb_found_by_bowtie <- nrow(mis_data[mis_data$kraken == 0, ])
  nb_found_by_kraken <- nrow(mis_data[mis_data$bowtie == 0, ])
  nb_mismatched_reads <- nrow(mis_data) - nb_found_by_bowtie - nb_found_by_kraken
  
  # graph
  vd1 <- euler(c("Bowtie" = nb_found_by_bowtie, 
                 "Kraken" = nb_found_by_kraken, 
                 "Bowtie&Kraken" = nb_matched_reads + nb_mismatched_reads))
  
  plot(vd1, counts = TRUE,
       fills = c("#80b1d3","#8EDE68"),
       legend = list(space = "right", 
                     columns=1,
                     labels = c("Alignés par Bowtie",
                                "Alignés par Kraken")),
       opacity = .7, 
       main = paste0(kefir, " : alignement des couples de reads"),
       quantities = list(TRUE, type = "percent"))
}


venn_mismatchs <- function(file_mismatch, kefir) {
  # lecture des données
  mis_data <- read.table(file_mismatch, sep = "\t", header = TRUE, skip = 1)[, -1]
  nb_total_reads <- read.table(file_mismatch, nrows = 1, sep = "\t")[, 2]
  
  # distribution
  nb_matched_reads <- nb_total_reads - nrow(mis_data)
  nb_found_by_bowtie <- nrow(mis_data[mis_data$kraken == 0, ])
  nb_found_by_kraken <- nrow(mis_data[mis_data$bowtie == 0, ])
  nb_mismatched_reads <- nrow(mis_data) - nb_found_by_bowtie - nb_found_by_kraken
  
  # graph
  vd2 <- euler(c("Bowtie" = nb_found_by_bowtie, 
                 "Kraken" = nb_found_by_kraken, 
                 "Bowtie&Kraken" = nb_matched_reads + nb_mismatched_reads,
                 "Organismes distincts&Bowtie&Kraken" = nb_mismatched_reads))
  
  plot(vd2, counts = TRUE,
       fills = c("#80b1d3","#8EDE68", "#DED268"),
       legend = list(space= "right", 
                     columns=1,
                     labels = c("Alignés par Bowtie",
                                "Alignés par Kraken",
                                "Lectures assignées à deux\norganismes différents")),
       opacity = .7, 
       main = paste0(kefir, " : alignement des couples de reads"),
       quantities = list(TRUE, type = "percent"))
}


mismatchs_kA <- "./test_kraken_bowtie/soutenance/mismatchs_kA.txt"
mismatchs_mcf <- "./test_kraken_bowtie/soutenance/mismatchs_mcf.txt"

venn_simple(mismatchs_kA, "Kéfir A")
venn_mismatchs(mismatchs_kA, "Kéfir A")

venn_simple(mismatchs_mcf, "Kéfir MCF")
venn_mismatchs(mismatchs_mcf, "Kéfir MCF")





#####################################
####### ALIGNE QUE PAR KRAKEN #######
#####################################
distrib_species_kraken <- function(species_info, file_mismatchs, kefir) {
  mis_data <- read.table(file_mismatchs, header = TRUE, sep = "\t", skip = 1)
  kr_only <- merge(mis_data[mis_data$bowtie == 0, ],
                   species_info,
                   by = "kraken")
  
  # Transformation des données
  prop_kr_only <- table(kr_only$species) / length(kr_only$species)
  df <- as.data.frame(prop_kr_only)
  colnames(df) <- c("species", "proportion")
  df$species <- gsub(" ", "\n", df$species)
  df <- df[order(df$species),]
  
  # plot
  ggplot(df, aes(x = fct_inorder(species), y = proportion, fill = fct_inorder(species))) +
    geom_bar(stat = "identity") +
    labs(title = paste0("Distribution des espèces du ", kefir, " \nalignées uniquement par Kraken2"),
         x = "",
         y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
}


distrib_species_match <- function(species_info, file_all, kefir) {
  all_data <- read.table(file_all, header = TRUE, sep = "\t")
  matched <- merge(all_data[all_data$bowtie == all_data$kraken, ],
                   species_info,
                   by = "kraken")
  
  # Transformation des données
  prop_matched <- table(matched$species) / length(matched$species)
  df <- as.data.frame(prop_matched)
  colnames(df) <- c("species", "proportion")
  df$species <- gsub(" ", "\n", df$species)
  df <- df[order(df$species),]
  
  # plot
  ggplot(df, aes(x = fct_inorder(species), y = proportion, fill = fct_inorder(species))) +
    geom_bar(stat = "identity") +
    labs(title = paste0("Distribution des espèces du ", kefir, " \nalignées par Kraken2 et Bowtie2"),
         x = "",
         y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
}



species_info <- read.table("./distrib_species/metagenome_ncbi_id.txt", sep = "\t")
colnames(species_info) <- c("library", "kraken", "species", "type")

distrib_species_kraken(species_info, mismatchs_kA, "Kéfir A")
file_all_ka <- "./test_kraken_bowtie/soutenance/reads_kA.txt"
distrib_species_match(species_info, file_all_ka, "Kéfir A")

distrib_species_kraken(species_info, mismatchs_mcf, "Kéfir MCF")
file_all_mcf <- "./test_kraken_bowtie/soutenance/reads_mcf.txt"
distrib_species_match(species_info, file_all_mcf, "Kéfir MCF")


####################################
####### DISTRIBUTION GLOBALE #######
####################################
distrib_species <- function(species_info, file_all, kefir, program) {
  all_data <- read.table(file_all, header = TRUE, sep = "\t")
  dat <- all_data[all_data[, program] != 0, ]
  colnames(dat)[which(names(dat) == program)] <- "num"
  sp <- merge(dat, species_info, by = "num")
  
  # Transformation des données
  prop <- table(sp$species) / length(sp$species)
  df <- as.data.frame(prop)
  colnames(df) <- c("species", "proportion")
  df$species <- gsub(" ", "\n", df$species)
  df <- df[order(df$species),]
  
  # plot
  ggplot(df, aes(x = fct_inorder(species), y = proportion, fill = fct_inorder(species))) +
    geom_bar(stat = "identity") +
    labs(title = paste0("Distribution des espèces du ", kefir, " \nalignées par ", program),
         x = "",
         y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
}

species_info <- read.table("./distrib_species/metagenome_ncbi_id.txt", sep = "\t")
colnames(species_info) <- c("library", "num", "species", "type")

distrib_species(species_info, "./test_kraken_bowtie/soutenance/reads_kA.txt", "Kéfir A", "bowtie")
distrib_species(species_info, "./test_kraken_bowtie/soutenance/reads_kA.txt", "Kéfir A", "kraken")

distrib_species(species_info, "./test_kraken_bowtie/soutenance/reads_mcf.txt", "Kéfir MCF", "bowtie")
distrib_species(species_info, "./test_kraken_bowtie/soutenance/reads_mcf.txt", "Kéfir MCF", "kraken")
# Le fichier texte reads_mcf est trop gros, j'ai créé la table de comptage en bash directement
# avec grep -v "^seq id" reads_mcf.txt | awk '$3 != 0 {print $3}' | sort | uniq -c > compt_mcf_kraken.txt


compt_mcf_kraken <- read.table("./test_kraken_bowtie/soutenance/compt_mcf_kraken.txt", header = FALSE)
colnames(compt_mcf_kraken) <- c("compt", "num")
sp <- merge(compt_mcf_kraken, species_info, by = "num")
sp <- sp[, c("compt", "species")]
sp[, "prop"] <- sp$compt / sum(sp$compt)
sp$species <- gsub(" ", "\n", sp$species)
sp <- sp[order(sp$species),]

ggplot(sp, aes(x = fct_inorder(species), y = prop, fill = fct_inorder(species))) +
  geom_bar(stat = "identity") +
  labs(title = "Distribution des espèces du Kéfir MCF\nalignées par kraken",
       x = "",
       y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")


compt_mcf_kraken <- read.table("./test_kraken_bowtie/soutenance/compt_mcf_bowtie.txt", header = FALSE)
colnames(compt_mcf_kraken) <- c("compt", "num")
sp <- merge(compt_mcf_kraken, species_info, by = "num")
sp <- sp[, c("compt", "species")]
sp[, "prop"] <- sp$compt / sum(sp$compt)
sp$species <- gsub(" ", "\n", sp$species)
sp <- sp[order(sp$species),]

ggplot(sp, aes(x = fct_inorder(species), y = prop, fill = fct_inorder(species))) +
  geom_bar(stat = "identity") +
  labs(title = "Distribution des espèces du Kéfir MCF\nalignées par bowtie",
       x = "",
       y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")




