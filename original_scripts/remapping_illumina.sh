# fastq file to analyse
R1="/home/ajoubert/Biodiversite/Illumina/BA/good/NG-24463_kefirba_lib385501_6726_2_1.good.fq";
R2="/home/ajoubert/Biodiversite/Illumina/BA/good/NG-24463_kefirba_lib385501_6726_2_2.good.fq";
READ_LG=151;
# folder with Bx/Yx.fa
# for each species, should contain .fa / .genome files
FOLDER_SPECIES="/home/ajoubert/Analyses_papier/Librairie_MCF/sed"; #/sed/Bowtie/Bowtie_Librairie_MCF_plus";

# INDEX LIBRAIRIE METAGENOME - illumina
# produced independently of that script by :
# bowtie2-build reference.fasta index_name;
BT_BUID="/home/ajoubert/Analyses_papier/Librairie_MCF/sed/Bowtie/Bowtie_Librairie_MCF_plus"; 





#### DO NOT TOUCH
# FILES
#------------
INFILE=`basename $R1`;
BUILD_SHORT=`basename $BT_BUID`;
R1_NAME=`echo $INFILE | cut -d "." -f 1`;
NAME=${R1_NAME/"_R1"/""};
BUILD_NAME=`echo $BUILD_SHORT | cut -d "." -f 1`;

# mapping against metagenome
SAM=$NAME"_"$BUILD_NAME".sam";
SORT=$NAME"_"$BUILD_NAME"_sort.sam";
BH=$NAME"_"$BUILD_NAME"_BH.sam";
UNMAP=$NAME"_"$BUILD_NAME"_unmapped.txt";
UNMAP_FA_R1=$NAME"_"$BUILD_NAME"_unmapped_R1.fasta";
UNMAP_FA_R2=$NAME"_"$BUILD_NAME"_unmapped_R2.fasta";
UNMAP_QUAST="QUAST_"$NAME"_"$BUILD_NAME"_unmapped";


# SCRIPTS
#-----------
BOWTIE="/home/duvernois/TOOLS/bowtie2-2.5.1-linux-x86_64/bowtie2";
SAMTOOLS="samtools";
SEQTK="/home/duvernois/TOOLS/seqtk/seqtk";
QUAST="python3 /home/duvernois/TOOLS/quast-5.2.0/quast.py";
BEDCOV="genomeCoverageBed";
BAMSTAT="/home/duvernois/TOOLS/jvarkit/dist/bamstats04.jar";
PICARD="/home/duvernois/TOOLS/picard/build/libs/picard.jar";
JVARKIT="/home/duvernois/TOOLS/jvarkit/dist/jvarkit.jar";
#----------------#
###  LET'S GO  ###
#----------------#


# mapping against metaG
($BOWTIE -p 10 -x $BT_BUID -1 $R1 -2 $R2 -S $SAM) 2> log_bowtie2;

## rescue of unmapped reads
$SAMTOOLS view -f 4 $SAM | awk -F "\t" '{print $1}' | sort -u > $UNMAP;
$SEQTK subseq $R1  $UNMAP | $SEQTK seq -a - > $UNMAP_FA_R1;
$SEQTK subseq $R2  $UNMAP | $SEQTK seq -a - > $UNMAP_FA_R2;

echo "** Number of unmapped reads (mate - to *2)" >> metrics.txt;
wc -l $UNMAP >> metrics.txt;

## Analysis per species
species=$($SAMTOOLS view -H $SAM | awk -F "\t" '{split($2,a,":"); split(a[2],b,"_");
    if (a[1] == "SN") print b[1]}' | sort -u);

for sp in $species; do
    echo $sp;

    CONCORD=$NAME"_"$BUILD_NAME"_concSH_"$sp".txt";
    CONCORD_R1=$NAME"_"$BUILD_NAME"_concSH_"$sp"_R1.fastq";
    CONCORD_R2=$NAME"_"$BUILD_NAME"_concSH_"$sp"_R2.fastq";

    awk -F "\t" -v s="$sp" '{split($2,b,":"); split(b[2],c,"_"); split($3,a,"_");
        if ((a[1]==s) || ((b[1] == "SN") && (c[1] == s))) print $0}' $SAM > out ;

    # reads concordants + single hit
    $SAMTOOLS view -F 4 out | awk -F "\t" '{print $1}' | sort | uniq -c |
        awk -F " " '{if ($1 == "2") print $2}' > $CONCORD;

    $SEQTK subseq $R1 $CONCORD > $CONCORD_R1;
    $SEQTK subseq $R2 $CONCORD > $CONCORD_R2;

    METRIC=$sp"_metrics.txt";
    echo "Library\t$sp" > $METRIC;
    LG=$(wc -l $CONCORD_R1 | awk -F " " '{print $1}');
    COUNT=$(echo "Number mate\t$LG")
    echo $COUNT >> $METRIC;

    READS=$(expr $LG \*  2);
    COUNT2=$(echo "Number reads\t$READS")
    echo $COUNT2 >> $METRIC;

    NT=$(expr $READS \*  $READ_LG);
    COUNT_NT=$(echo "Number mapped nt\t$NT")
    echo $COUNT_NT >> $METRIC;

    # remapping on reference genome
    GENOME=$FOLDER_SPECIES"/"$sp".fa";
    LG=$FOLDER_SPECIES"/"$sp".genome";

    REMAP_BAM=$NAME"_"$BUILD_NAME"_concSH_"$sp".bam";
    COVFILE=$NAME"_"$BUILD_NAME"_concSH_"$sp".coverage";
    COVIMG=$NAME"_"$BUILD_NAME"_concSH_"$sp".svg";

    java -jar $PICARD FilterSamReads I=$SAM O=int.sam FILTER=includeReadList READ_LIST_FILE=$CONCORD ;
    awk -F "\t" -v s="$sp" '{split($2, a, ":"); split(a[2],b,"_");
        if ((a[1] == "SN") && (b[1] == s)) print $0; else if (a[1] != "SN") print $0}' int.sam > int2.sam;
    $SAMTOOLS view -bS int2.sam | $SAMTOOLS sort - >  $REMAP_BAM;
    rm int2.sam;
    $SAMTOOLS index $REMAP_BAM;

    # creation genome reference data
    $SAMTOOLS faidx $GENOME;
    awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $GENOME.fai > $GENOME.bed;
    java -jar $BAMSTAT -B  $GENOME.bed $REMAP_BAM > $COVFILE;
    MEDCOV=$(cat $COVFILE | sed '1d' | awk -F "\t" '{print $9}'| sort -k 1n,1 | tail -n1| cut -d "." -f 1);
    if  [ $MEDCOV -eq 0 ]
    then
        MEDCOV=5;
    fi
    RATIO=$(expr $MEDCOV \*  30 / 100); # 30%
    MEDGRAPH=$(expr $MEDCOV + $RATIO);

    java -jar $PICARD CreateSequenceDictionary R=$GENOME O=$GENOME.dict;
    java -jar $JVARKIT wgscoverageplotter -C $MEDGRAPH -R $GENOME $REMAP_BAM -o $COVIMG;

done


paste QUAST_*/report.tsv > int;
COL=$(head -1 int | awk '{ FS = "\t" } ; { print NF}' );

CMD="cut -f 1"
for i in $(seq 2 2 $COL); do
    CMD=$CMD","$i;
done

echo "\n** QUAST metrics\n" >> metrics.txt;
$CMD int >> metrics.txt ;
rm int;
