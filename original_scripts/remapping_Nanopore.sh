# fastq file to analyse
INPUT="/home/ajoubert/Analyses_papier/01_03_2023/01_03_2023_demultiplexing_pass/01_03_2023_G10/01_03_2023_G10hp.fastq";
# folder with Bx/Yx.fa
# for each species, should contain .fa / .genome files
FOLDER_SPECIES="/home/ajoubert/Analyses_papier/Librairie_MCF/sed";

# INDEX LIBRAIRIE METAGENOME
# produced independently of that script by :
# minimap2 -x map-ont -d Librairie.mmi Librairie.fa;
MMI="/home/ajoubert/Analyses_papier/Librairie_MCF/sed/Librairie_MCF.mmi";






#### DO NOT TOUCH
# FILES
#------------
INFILE=`basename $INPUT`;
MMI_SHORT=`basename $MMI`;
NAME=`echo $INFILE | cut -d "." -f 1`;
MMI_NAME=`echo $MMI_SHORT | cut -d "." -f 1`;

# mapping against metagenome
SAM=$NAME"_"$MMI_NAME".sam";
SORT=$NAME"_"$MMI_NAME"_sort.sam";
BH=$NAME"_"$MMI_NAME"_BH.sam";
UNMAP=$NAME"_"$MMI_NAME"_unmapped.txt";
UNMAP_FA=$NAME"_"$MMI_NAME"_unmapped.fasta";
UNMAP_QUAST="QUAST_"$NAME"_"$MMI_NAME"_unmapped";


# SCRIPTS
#-----------
MINIMAP="minimap2";
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
$MINIMAP -a -t 15 -x map-ont $MMI $INPUT > $SAM;
$SAMTOOLS view -H $SAM > header.sam;
grep -v '^@' $SAM | sort -k 1,1 > $SORT;

# get BH
awk -F "AS:i:" 'BEGIN{cs = -10; ID =""; line= ""} {split($2,a,"\t"); split($1,b,"\t");
    if (ID != b[1]){
        if (line != ""){
            print(line)
        }
        ID = b[1]
        cs = a[1]
        line=$0;
    } else {
        if (cs <a[1]){
            cs = a[1]
            line = $0;
        }
        else {
            line=line;
        }
    }
} END {print line}' $SORT > int;
cat header.sam int > $BH;
rm int;
rm header.sam;

## rescue of unmapped reads
$SAMTOOLS view -f 4 $SAM | awk -F "\t" '{print $1}' > $UNMAP;
$SEQTK subseq $INPUT  $UNMAP | $SEQTK seq -a - > $UNMAP_FA;
$QUAST --min-contig=1 -o $UNMAP_QUAST $UNMAP_FA;

echo "** Analysis of unmapped reads" >> metrics.txt;
wc -l $UNMAP >> metrics.txt;

## Analysis per species
species=$($SAMTOOLS view -H $SAM | awk -F "\t" '{split($2,a,":"); split(a[2],b,"_");
    if (a[1] == "SN") print b[1]}' | sort -u);

for sp in $species; do
    echo $sp;
    OUT=$NAME"_"$MMI_NAME"_BH_"$sp".txt";
    FQ=$NAME"_"$MMI_NAME"_BH_"$sp".fastq";
    FA=$NAME"_"$MMI_NAME"_BH_"$sp".fasta";
    QUAST_FILE="QUAST_"$NAME"_"$MMI_NAME"_BH_"$sp;

    awk -F "\t" -v s="$sp" '{split($3,a,"_");
        if (a[1]==s) print $1}' $BH > $OUT;
    $SEQTK subseq $INPUT $OUT > $FQ;
    $SEQTK subseq $INPUT $OUT | $SEQTK seq -a - > $FA;
    $QUAST --min-contig=1 -o $QUAST_FILE $FA;

    # remapping on reference genome
    GENOME=$FOLDER_SPECIES"/"$sp".fa";
    LG=$FOLDER_SPECIES"/"$sp".genome";
    REMAP=$NAME"_"$MMI_NAME"_BH_"$sp"_remapping.sam";
    REMAP_BAM=$NAME"_"$MMI_NAME"_BH_"$sp"_remapping.bam";
    COVFILE=$NAME"_"$MMI_NAME"_BH_"$sp"_remapping.coverage";
    COVIMG=$NAME"_"$MMI_NAME"_BH_"$sp"_remapping.svg";

    $MINIMAP -a -x map-ont $GENOME $FQ > $REMAP;
    $SAMTOOLS view -bS $REMAP | $SAMTOOLS sort - > $REMAP_BAM;
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
