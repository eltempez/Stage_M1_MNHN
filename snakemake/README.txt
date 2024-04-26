Cette routine snakemake est à effectuer sur le cluster de calcul, avec comme ligne de commande :
>snakemake -s <chemin d'accès au fichier> --core 1 --use-envmodules --use-conda --configfile <chemin d'accès au fichier de configuration>

Prérequis : besoin d'activer snakemake à l'aide de la commande :
>conda activate snakemake7

Les seules choses à modifier sont les variables du fichier de configuration config.yaml, ainsi que le contenu du fichier metagenome_ncbi_id.txt si besoin

Fonctionnement du fichier metagenome_ncbi_id.txt : 
Ce fichier répertorie les numéros d'accession des organismes présents dans le métagénome. les colonnes sont séparées par un "tab". La 3ème colonne (suivant le #) est facultative.
ATTENTION il est inutile de supprimer des lignes, la routine ne prendra en compte que les organismes présents dans le fichier folder_species