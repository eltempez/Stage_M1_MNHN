Cette routine snakemake est à effectuer sur le cluster de calcul, avec comme ligne de commande :

snakemake -s <chemin d'accès au fichier remapping_illu.smk> --core 1 --use-envmodules --use-conda --configfile <chemin d'accès au fichier de configuration config.yaml> --directory /mnt/beegfs/<username>

(Pour faire tourner la routine plusieurs fois en parallèle, changer le combre de cores avec --core 2)

Prérequis : besoin d'activer snakemake à l'aide de la commande :

conda activate snakemake7


Le dossier config_file doit contenir le fichier de configuration config.yaml, le fichier metagenome_ncbi_id.txt, ainsi que 3 dossiers contenant les scripts nécessaires : yml, Rmd et py.

Les seules choses à modifier sont les variables à l'intérieur du fichier de configuration config.yaml, ainsi que le contenu du fichier metagenome_ncbi_id.txt si besoin


Fonctionnement du fichier metagenome_ncbi_id.txt :

Ce fichier répertorie les numéros d'accession des organismes présents dans le métagénome. les colonnes sont séparées par un "tab" et sont : n° d'accession NCBI, numéro interne kéfir, espèce (séparé par espace), type d'organisme.
ATTENTION il est inutile de supprimer des lignes, la routine ne prendra en compte que les organismes présents dans le fichier folder_species, sous condition que le nom de fichier commence par le même nom que la première colonne de metagenome_ncbi_id.txt