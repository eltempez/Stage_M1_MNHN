library(tidyverse)
library(RColorBrewer)

# Import
data_10h <- read.table("./distrib_species/global_metrics_10h.txt", skip = 1, header = T, sep = "\t")
data_20h <- cbind(library = data_10h[, 1], data_10h[, -1]*2)
  
  


v1 <- c(2, 2, 10)
v2 <- c(3, 2, 5)
v3 <- c(5, 2, 5)
data <- data.frame(v1, v2, v3, row.names = c("a", "b", "c"))


## Transformation des données en format long
data_long <- data %>%
  # Ajout des noms de lignes en colonne
  rownames_to_column(var = "espece") %>%
  # Ajout des colonnes en format long, sauf l'espèce
  pivot_longer(cols = -espece, names_to = "echantillon", values_to = "value")

## Visualisation en cercles
max_val <- max(data_long$value)
min_val <- min(data_long$value)
ggplot(data_long, aes(x = echantillon, y = espece, size = value, col = espece)) +
  geom_count() +
  guides(col = FALSE) +
  scale_size_area(max_size = 30) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Echantillon", y = "Espèce", size = "Nombre de \nreads")

## Visualisation en barplot
sum <- colSums(data)
data_f <- sweep(data, 2, sum, FUN = '/')
data_long_f <- data_f %>%
  # Ajout des noms de lignes en colonne
  rownames_to_column(var = "espece") %>%
  # Ajout des colonnes en format long, sauf l'espèce
  pivot_longer(cols = -espece, names_to = "echantillon", values_to = "freq")
# visualisation
ggplot(data_long_f, aes(x = echantillon, y = freq, fill = espece)) +
  geom_col() +
  scale_fill_brewer(palette = "Dark2")


