#!/bin/bash

# Créer une boucle for pour itérer sur la plage de valeurs de 25 à 40 par pas de 5
for ((i=25; i<=40; i+=5))
do
    l=$(($i - 4))
    sp=$(($l / 4 - 1))
    if [$l -gt 31 ]; then
        $l=31
    fi
    if [$sp -gt 7 ]; then
        $sp=7
    fi
    echo $i
    # données taxonomie
    /home/eliotttempez/Documents/tests/kraken_kmers/kraken2-master/kraken2-build --download-taxonomy --db /home/eliotttempez/Documents/tests/kraken_kmers/kefir_library
    # ajout des génomes à la librairie
    for file in /home/eliotttempez/Documents/tests/outputs_mcf10_snk/genome_data/fa/*.fa
    do
        /home/eliotttempez/Documents/tests/kraken_kmers/kraken2-master/kraken2-build --add-to-library $file --db /home/eliotttempez/Documents/tests/kraken_kmers/kefir_library
    done
    # construction de la librairie avec la longueur de k-mers correspondante
    /home/eliotttempez/Documents/tests/kraken_kmers/kraken2-master/kraken2-build --build --db /home/eliotttempez/Documents/tests/kraken_kmers/kefir_library --kmer-len $i --minimizer-len $l --minimizer-spaces $sp

    # alignement
    /home/eliotttempez/Documents/tests/kraken_kmers/kraken2-master/kraken2 --db /home/eliotttempez/Documents/tests/kraken_kmers/kefir_library --threads 8 --paired /home/eliotttempez/Documents/donnees/DATAS_cinétique_MCF/Cinétique_illumina_MCF/NG-32390_MCF_10h_lib666597_10174_2_1.fastq.gz /home/eliotttempez/Documents/donnees/DATAS_cinétique_MCF/Cinétique_illumina_MCF/NG-32390_MCF_10h_lib666597_10174_2_2.fastq.gz > /home/eliotttempez/Documents/tests/kraken_kmers/kraken_output_${i}_kmers.txt

    rm -rf /home/eliotttempez/Documents/tests/kraken_kmers/kefir_library
done
