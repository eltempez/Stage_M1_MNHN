Cette routine snakemake est à effectuer sur le cluster da calcul, avec comme ligne de commande snakemake -s <chemin d'accès au fichier> --core 1 --use-envmodules

Les lignes à modifier sont :
- les chemins d'accès aux différents fichiers dans le fichier de configuration config.yaml
- le chemin d'accès au fichier de configuration dans remapping_illu.smk