# DATA TO MODIFY
#----------------
##please comment with sharp the wrong type of data
#TYPE="PE"; #paired-end
TYPE="SE"; #single-end
##complete path to fastq file: assembly input file
# paired-end librairies R1/R2
R1="/home/ajoubert/DATAS/Illumina/MCF/Cinetique/NG-32390_MCF_30h_lib666598_10168_3_1.fastq";
R2="/home/ajoubert/DATAS/Illumina/MCF/Cinetique/NG-32390_MCF_30h_lib666598_10168_3_2.fastq";
# single-end librairy
INPUT="/home/ajoubert/Analyse_Cinetique/MCF_48h_p.fastq";

OUTPUT_NAME="MCF_48h_N"; # no extension - name used in the rest of the outputs data
OUTPUT_FOLDER="/home/ajoubert/KRAKEN/48h/";

DATABASE="KEFIR";
#~ DATABASE="k2protocol_db"; # A CHOISIR - base de ref de Lu et al, 2022, Nature Protocole



### SCRIPT DO NOT TOUCH
#----------------------
KRAKEN_OUT="kraken2_"$OUTPUT_NAME".txt";
EVOLVED="evolved_"$KRAKEN_OUT;
KRONA_OUT="kraken2_"$OUTPUT_NAME".html";
UNATTR="kraken2_unattributed_"$OUTPUT_NAME".txt";
UNATTR_READS="kraken2_unattributed_"$OUTPUT_NAME".fa";
UNATTR_R1="kraken2_unattributed_"$OUTPUT_NAME"_R1.fa";
UNATTR_R2="kraken2_unattributed_"$OUTPUT_NAME"_R2.fa";

# script
export KRAKEN2_DB_PATH="/home/duvernois/DATABASES/KRAKEN/:$KRAKEN2_DB_PATH";
KRAKEN="/home/duvernois/TOOLS/kraken2/kraken2";
KRONA="ktImportTaxonomy";
TAXPATH="/home/duvernois/TOOLS/Krona/KronaTools/taxonomy";
seqtk="/home/duvernois/TOOLS/seqtk/seqtk";

###########################
# *********************** #
###########################
CURRENT_FOLDER=$PWD;
#~ mkdir $OUTPUT_FOLDER;
cd $OUTPUT_FOLDER;

##  launching of kraken
echo "kraken \r"
#k2protocol_db
if [ $TYPE = "PE" ]; then
        echo "KRAKEN IN PAIRED END MODE";
        $KRAKEN --db KEFIR --paired $R1 $R2  --threads 6 > $KRAKEN_OUT;
else
        echo "KRAKEN IN SINGLE END MODE";
        $KRAKEN --db KEFIR $INPUT  --threads 6 > $KRAKEN_OUT;
fi

## reformating
cat $KRAKEN_OUT | cut -f 3,4 > $EVOLVED;

## krona
$KRONA $EVOLVED -o $KRONA_OUT -t 1 -s 2 -tax $TAXPATH;

# get unattributed reads
cat $KRAKEN_OUT | awk -F "\t" '{if ($3 == 0) print $2}' > $UNATTR;
if [ $TYPE = "PE" ]; then
        echo "GET UNATTRIBUTED READS IN PAIRED END MODE";
        $seqtk subseq $R1 $UNATTR > $UNATTR_R1;
        $seqtk subseq $R2 $UNATTR > $UNATTR_R2;
else
        echo "GET UNATTRIBUTED READS IN SINGLE END MODE";
        $seqtk subseq $INPUT $UNATTR > $UNATTR_READS;
fi

