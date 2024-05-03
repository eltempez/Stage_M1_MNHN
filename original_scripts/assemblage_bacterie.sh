# DATA TO MODIFY
##complete path to fastq file: assembly input file
INPUT="/home/ajoubert/Assemblage_01082023/FAL18593/FAL18593_demultiplexing_pass/barcode01.fastq";
USE_ILLUMINA=1; # 2 == FALSE
ILLUMINA_R1="/home/ajoubert/Assemblage_01082023/Illumina_good/NG-34158_B1_lib717418_10289_3_1.good.fq";
ILLUMINA_R2="/home/ajoubert/Assemblage_01082023/Illumina_good/NG-34158_B1_lib717418_10289_3_2.good.fq";
OUTPUT_FOLDER="/home/ajoubert/Assemblage_01082023/Assemblage_B1_final";
#GENOME_SIZE="2.9m";
#nanopore kit
KIT="ont-ligation";
# complete path to reference genome(s) - concatenation of genomes if needed
REF_GENOME_INIT="/home/ajoubert/Assemblage_01082023/Genome_Proteome_ref/GCF_000583855.1_ASM58385v1_genomic.fna";
# proteome
REF_PROTEOME_INIT="/home/ajoubert/Assemblage_01082023/Genome_Proteome_ref/GCF_000583855.1_ASM58385v1_protein.faa";
# GC content
GC=61 ; ## check to keep it automatically
DFAST_ORG="B1";

### SCRIPT DO NOT TOUCH
#----------------------
export BUSCO_CONFIG_FILE="/home/duvernois/TOOLS/busco/config/config.ini"
export PATH="/home/duvernois/TOOLS/busco/build/lib/busco/:$PATH";
export PATH="/home/duvernois/TOOLS/metaeuk/bin/:$PATH";
export PATH="/usr/bin/:$PATH";
#~ export PATH="/media/bigdata/duvernois/SOFTS/samtools-1.12/:$PATH";

INFILE=`basename $INPUT`;
NAME=`echo $INFILE | cut -d "." -f 1`;
REFFILE=`basename $REF_GENOME_INIT`;
NAME_REF=`echo $REFFILE | cut -d "." -f 1`;
REFPROT=`basename $REF_PROTEOME_INIT`;
TITLEDB=`echo $REFPROT | cut -d "." -f 1`;
# QC files
QC1=$NAME"_longQC";
QC2=$NAME"_porechop_longQC";
# porechop output file
CLEAN_ASSEMBLY=$NAME"_porechop.fastq";
REF_GENOME="clean_"$REFFILE;
REF_PROTEOME="clean_"$REFPROT;
# FLYE
PREF=$NAME"_porechop"; #prefix
FLYE_FOLDER="FLYE_"$NAME"_porechop"; # output folder for FLYE
FLYE_ASSEMBLY=$FLYE_FOLDER"/"$PREF".contigs.fasta";

# POLISHING
if [ $USE_ILLUMINA -eq 1 ]
then
    # minimap
    #~ REMAP_FILE=$NAME"_consensusRaconMedaka_porechop_minimap.bam";
    PAF_FILE=$NAME"_IlluminaR1_overlap.paf";
    PAF_FILE1=$NAME"_IlluminaR1_overlap1.paf";
    #racon
    RACON_FST=$PREF"_assembly_racon.fasta";
    RACON_FST1=$PREF"_assembly_racon1.fasta";
    #medaka
    MEDAKA_ASSEMBLY=$PREF"_consensusRaconMedaka.fasta";
    #busco
    BUSCO_FILE=$NAME"_assembly_raconMedaka.busco";
    # nucmer
    PREF_NUCMER=$NAME"_consensusRaconMedaka_vs_"$NAME_REF"_nucmer";
    PREFREV_NUCMER=$NAME_REF"_vs_"$NAME"_consensusRaconMedaka_nucmer";
    # glimmer
    ICM=$PREF"_consensusRaconMedaka.contigs.icm";
    longorfs=$PREF"_consensusRaconMedaka.contigs.longorfs";
    train=$PREF"_consensusRaconMedaka.contigs.train";
    ORFS=$PREF"_consensusRaconMedaka.contigs.ORFs";
    PREDICT=$PREF"_consensusRaconMedaka.contigs.ORFs.predict";
    #GFF_AUG=$PREF"_consensusRacon.contigs.augustus.gff";
    GFF=$PREF"_consensusRaconMedaka.contigs.genes.gff";
    ORFs_FA=$PREF"_consensusRaconMedaka.contigs.genes.fa";
    PEP_FA=$PREF"_consensusRaconMedaka.contigs.pep.fa";
    #blast
    BLASTP_ANNOT=$PREF"_consensusRaconMedaka.contigs.pep.blastp";
    BLASTP_GFF=$PREF"_consensusRaconMedaka.contigs.pep.blastp.gff";
    BLASTP_GB=$PREF"_consensusRaconMedaka.contigs.pep.blastp";
    BLASTX_ANNOT=$PREF"_consensusRaconMedaka.contigs.pep.blastx";
    BLASTX_GFF=$PREF"_consensusRaconMedaka.contigs.pep.blastx.gff";
    BLASTX_GB=$PREF"_consensusRaconMedaka.contigs.pep.blastx";
    TBLASTN_ANNOT=$PREF"_consensusRaconMedaka.contigs.pep.tblastn";
    TBLASTN_GFF=$PREF"_consensusRaconMedaka.contigs.pep.tblastn.gff";
    TBLASTN_GB=$PREF"_consensusRaconMedaka.contigs.pep.tblastn";
    ASSEMBLY_GB=$PREF"_consensusRaconMedaka.contigs";
else
    #medaka
    MEDAKA_ASSEMBLY=$PREF"_consensus_noIll.fasta";
    #busco
    BUSCO_FILE=$NAME"_assembly__noIll.busco";
    # nucmer
    PREF_NUCMER=$NAME"_consensus_noIll_vs_"$NAME_REF"_nucmer";
    PREFREV_NUCMER=$NAME_REF"_vs_"$NAME"_consensus_noIll_nucmer";
    # glimmer
    ICM=$PREF"_consensus_noIll.contigs.icm";
    longorfs=$PREF"_consensus_noIll.contigs.longorfs";
    train=$PREF"_consensus_noIll.contigs.train";
    ORFS=$PREF"_consensus_noIll.contigs.ORFs";
    PREDICT=$PREF"_consensus_noIll.contigs.ORFs.predict";
    #GFF_AUG=$PREF"_consensusRacon.contigs.augustus.gff";
    GFF=$PREF"_consensus_noIll.contigs.genes.gff";
    ORFs_FA=$PREF"_consensus_noIll.contigs.genes.fa";
    PEP_FA=$PREF"_consensus_noIll.contigs.pep.fa";
    #blast
    BLASTP_ANNOT=$PREF"_consensus_noIll.contigs.pep.blastp";
    BLASTP_GFF=$PREF"_consensus_noIll.contigs.pep.blastp.gff";
    BLASTP_GB=$PREF"_consensus_noIll.contigs.pep.blastp";
    BLASTX_ANNOT=$PREF"_consensus_noIll.contigs.pep.blastx";
    BLASTX_GFF=$PREF"_consensus_noIll.contigs.pep.blastx.gff";
    BLASTX_GB=$PREF"_consensus_noIll.contigs.pep.blastx";
    TBLASTN_ANNOT=$PREF"_consensus_noIll.contigs.pep.tblastn";
    TBLASTN_GFF=$PREF"_consensus_noIll.contigs.pep.tblastn.gff";
    TBLASTN_GB=$PREF"_consensus_noIll.contigs.pep.tblastn";
    ASSEMBLY_GB=$PREF"_consensus_noIll.contigs";
fi

# QUAST
QUAST_OUT="QUAST_"$MEDAKA_ASSEMBLY;

# DFAST OUTPUT
DFAST_OUT="DFAST_"$MEDAKA_ASSEMBLY;

# REMAPPING
READS_REMAP="remapping_"$MEDAKA_ASSEMBLY".bam";
GENOME_LG=$MEDAKA_ASSEMBLY".genome";
COVERAGE=$MEDAKA_ASSEMBLY".coverage";

# script
#---------
LONGQC="/home/duvernois/TOOLS/LongQC/longQC.py";
DIVRATE="/home/duvernois/TOOLS/UTILS_MAC/getDivRate.sh";
FLYE="/home/duvernois/TOOLS/Flye/bin/flye";
porechop="porechop";
RACON="/home/duvernois/TOOLS/racon/build/bin/racon";
MEDAKA="medaka_consensus";
QUAST="/home/duvernois/TOOLS/quast-5.2.0/quast.py";
nucmer="nucmer";
MUMMERPLOT="/home/duvernois/TOOLS/mummer-4.0.0rc1/mummerplot";
MINIMAP="minimap2";
BUSCO="/home/duvernois/TOOLS/busco/bin/busco";
PLOT_BUSCO="/home/duvernois/TOOLS/busco/scripts/generate_plot.py";

# glimmer prediction
GLIMMER_G3="/home/duvernois/TOOLS/glimmer3.02/sample-run/g3-from-scratch.csh";
GLIMMER_G3_IT="/home/duvernois/TOOLS/glimmer3.02/sample-run/g3-iterated.csh";
GLIMMMER_LONGORF="/home/duvernois/TOOLS/glimmer3.02/bin/long-orfs";
GLIMMMER_EXTRACT="/home/duvernois/TOOLS/glimmer3.02/bin/extract";
GLIMMER_ICM="/home/duvernois/TOOLS/glimmer3.02/bin/build-icm";
GLIMMER="/home/duvernois/TOOLS/glimmer3.02/bin/glimmer3";
GLI2FA="/home/duvernois/TOOLS/UTILS_MAC/glimmer2fasta.py";
# dfast prediction
DFAST="/home/duvernois/TOOLS/dfast_core-1.2.18/dfast ";
    #export PATH=$PATH:/home/duvernois/TOOLS/augustus/bin:/home/duvernois/TOOLS/augustus/scripts;
    #AUGUSTUS="augustus";
    #BEDTOOLS="bedtools";
    #BEDFA=$BEDTOOLS" getfasta";
    #BEDCOV=$BEDTOOLS" genomecov";
    #AUG2PEP="/home/duvernois/TOOLS/UTILS/augustus2fasta.py";
    #AUG2FA="/home/duvernois/TOOLS/UTILS/gff2seq.py";
MBDB="makeblastdb";
BLASTP="blastp";
BLASTX="blastx";
TBLASTN="tblastn";
#~ GETHEADER="/home/duvernois/TOOLS/UTILS/create_headerGFF_for_snapgene.py";
MERGEGB="/home/duvernois/TOOLS/UTILS_MAC/merge_gb_gff.py";
# samtools
SAMTOOLS="samtools";
seqtk="/home/duvernois/TOOLS/seqtk/seqtk";
CREATE_GENOME="/home/duvernois/TOOLS/UTILS/create_genome_file_for_BEDTOOLS.py";

###########################

###########################
# *********************** #
###########################
CURRENT_FOLDER=$PWD;
#mkdir $OUTPUT_FOLDER;
cd $OUTPUT_FOLDER;

## cleaning of files
echo "cleaning \r"
sed 's/\r//' $REF_GENOME_INIT > $REF_GENOME;

## filtering of data
echo "QC and filtering..."
python3 $LONGQC sampleqc -x $KIT -o $QC1 $INPUT;
$porechop -i $INPUT -o $CLEAN_ASSEMBLY;
python3 $LONGQC sampleqc -x $KIT -o $QC2 $CLEAN_ASSEMBLY;

## genome assembly
echo "FLYE assembly..."
### COVERAGE CHANGED car ici on repart des reads du metaG et pas nanopore --> pb de couverture !!!
$FLYE --nano-raw $CLEAN_ASSEMBLY --out-dir $FLYE_FOLDER --threads 8;
# --genome-size $GENOME_SIZE --asm-coverage 200;
cp $FLYE_FOLDER/assembly.fasta $FLYE_ASSEMBLY;

if [ $USE_ILLUMINA -eq 1 ]
then
    #~ ## polishing by racon
    echo "Polishing by Racon..."
    echo "\t round0 Illumina R1 reads"
    #realign by minimap
    $MINIMAP -x sr $FLYE_ASSEMBLY $ILLUMINA_R1 > $PAF_FILE;
    # racon sequence_fastq overlap_sam sequence_to_correct
    $RACON -t 16 $ILLUMINA_R1 $PAF_FILE $FLYE_ASSEMBLY > $RACON_FST;

    echo "\t round1 Illumina R1 reads"
    #realign by minimap
    $MINIMAP -x sr $RACON_FST $ILLUMINA_R1 > $PAF_FILE1;
    # racon sequence_fastq overlap_sam sequence_to_correct
    $RACON -t 16 $ILLUMINA_R1 $PAF_FILE1 $RACON_FST > $RACON_FST1;

    . /home/duvernois/TOOLS/medaka/bin/activate
    $MEDAKA -b 100 -i $ILLUMINA_R1 -d $RACON_FST1 -o medaka-results/ -t 10
    deactivate
    cp medaka-results/consensus.fasta $MEDAKA_ASSEMBLY;
else
    cp $FLYE_ASSEMBLY $MEDAKA_ASSEMBLY;
fi

echo "quast";
python3 $QUAST -o $QUAST_OUT $MEDAKA_ASSEMBLY;

echo "BUSCO"
$BUSCO  -i $MEDAKA_ASSEMBLY -m genome -l bacteria -c 10 --out $BUSCO_FILE;
python3 $PLOT_BUSCO -wd $BUSCO_FILE;


## comparison of all genome
echo "comparison with reference genome..."
$nucmer --prefix=$PREF_NUCMER $REF_GENOME $MEDAKA_ASSEMBLY;
$nucmer --prefix=$PREFREV_NUCMER $MEDAKA_ASSEMBLY $REF_GENOME;
show-coords -rcl $PREF_NUCMER.delta > $PREF_NUCMER.coords
$MUMMERPLOT --postscript --prefix=$PREF_NUCMER $PREF_NUCMER.delta -Q $MEDAKA_ASSEMBLY -R $REF_GENOME --filter --layout

## percentage of divergence
#bash $DIVRATE $PREF_NUCMER $PREFREV_NUCMER > out_divergence.txt;
## to do only if one chromosome and not several on REF genome
nb_chr=$(grep ">" $REF_GENOME | wc -l);
if [ $nb_chr==1 ]
then
        bash $DIVRATE $PREF_NUCMER $PREFREV_NUCMER > out_divergence.txt;
fi

#~ ## gene prediction
echo "gene prediction from scratch..."
$GLIMMER_G3 $MEDAKA_ASSEMBLY $PREF".contigs";
$GLIMMMER_LONGORF -n -t 1.15 $MEDAKA_ASSEMBLY $longorfs;
$GLIMMMER_EXTRACT -t $MEDAKA_ASSEMBLY $longorfs > $train;
$GLIMMER_ICM -r $ICM < $train;
$GLIMMER -o50 -g110 -t30 -C $GC $MEDAKA_ASSEMBLY $ICM $ORFS;
cat $PREDICT | awk '{OFS="\t"} {strand = "+"}{if($4 < 0) strand="-"}{gsub(/[+-]/," ")}{if ($1 ~ /^>/) chr=$1;
        else print chr, "GLIMMER3.02", "gene" , $2 , $3, $5, strand , $4, "ID="$1}' | sed 's/>//' > $GFF;
python3 $GLI2FA $PREDICT $MEDAKA_ASSEMBLY $ORFs_FA;
transeq $ORFs_FA $PEP_FA;

cat $REF_PROTEOME_INIT | sed 's/\r//'  > $REF_PROTEOME;
echo "Annotation running..."
$MBDB -in $REF_PROTEOME -dbtype prot -title $TITLEDB;
echo "\t\tblastx..."
$BLASTX -db $REF_PROTEOME -query $ORFs_FA -outfmt  "7 qacc sacc evalue qstart qend sstart send qseqid stitle" -out $BLASTX_ANNOT -max_target_seqs 5 -evalue 0.01;

# create tmp gff file with only annotated data
cat $BLASTX_ANNOT | grep -v "#" | sed 's/\./,/' | awk -F "\t" '{split($9, a, "["); split(a[1], b ,"."); split(b[2], c, ","); print $1   "\t" c[1]}' | awk '{ gsub(/[ ]+$/,""); print }' | sed 's/1 //' | sort -u | sed 's/ /%%/g' |  awk -F "\t" 'BEGIN{id=""; data = ""}{ if ($1 == id) {data = data "/"$2;} else { print id "\t" data; data = $2 ; id = $1;}} END{print id "\t" data}' |awk -F "\t" '{split($1,a,"_"); print a[1] "_" a[2] "\t" $0}' | sort -k 2,2 > bl.tmp;

# all gene ID
# create empty gff
cat $GFF | awk -F "\t" '{split($9,a, "="); print $1 "_" a[2] "\t" $0 ";Name=unknownFunction"}' | sort -k 1,1 > gff.tmp;
#~ grep ">" $ORFs_FA | sed 's/>//' |    awk -F " " '{split($1, a, "_"); split($2,b, "="); split($3,c,"="); split($4,d,"="); print $1 "\t" a[1] "\tGLIMMER3.02\tgene\t" b[2] "\t" c[2] "\t.\t" substr(d[2],1,1) "\t.\tID=" $1";Name=unknownFunction"}' | sort -k 1,1 > gff_init.tmp
# genes with annotation
awk -F "\t" '{if ($2 != "") print $2}' bl.tmp | sort -k 1,1 > ID_gffannot.tmp;

# cross blast results with previous GFF
join -j1 1 -j2 2 gff.tmp bl.tmp |  awk -F " " '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t ID=" $1 ";Name="$12}' | sed 's/%%/ /g' > bl_annot.tmp

# cross annotation and initial list of genes
join -j1 1 -j2 1 gff.tmp ID_gffannot.tmp | sed 's/ /\t/g' > join.tmp;
diff join.tmp gff.tmp | grep ">" | sed 's/> //' | awk -F "\t" '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' > missing.tmp

cat bl_annot.tmp missing.tmp >$BLASTX_GFF;


# assembly in genbank
echo "create genbank file...";
seqret -sequence $MEDAKA_ASSEMBLY -osformat genbank -osname_outseq $ASSEMBLY_GB -ofdirectory_outseq gbk_file -auto;
# tool to merge assembly in GB with features data
python3 $MERGEGB -f $BLASTX_GFF -g $ASSEMBLY_GB.genbank -n $PREF -o $BLASTX_GB.genbank;
rm *.tmp;


# rescue of unknow genes
grep "unknown" $BLASTX_GFF | awk -F "\t" '{split($9,a,";"); split(a[1], b, "="); print b[2]}' > ID_unknown.txt;
$seqtk subseq $ORFs_FA ID_unknown.txt > ID_unknown.fa;





####
## Gene prediction FOR BACTERIA
## with dfast
python3 $DFAST -g $MEDAKA_ASSEMBLY -o $DFAST_OUT --organism "$DFAST_ORG"

################################
#### GET COVERAGE OF MAPPING ###
################################
#1. remapping reads on assembly
$MINIMAP -a -x map-ont $MEDAKA_ASSEMBLY $CLEAN_ASSEMBLY  | samtools view -bS  - | samtools sort - > $READS_REMAP;

#2. bam and sort data
#~ B1_remapping.sam | samtools sort - > B1_remapping.bam

#3. create genome file
python $CREATE_GENOME -f $MEDAKA_ASSEMBLY  -o $GENOME_LG;

#4. get coverage file
genomeCoverageBed -ibam $READS_REMAP -g $GENOME > $COVERAGE;

# print coverage (dans le terminal - si je fais pas d'erreur de calcul)
awk -F "\t" 'BEGIN{contig = ""; cov = ""; lg = "";}{if ($1 == contig) {cov = cov + ($2*$3); lg = $4; } else { if (lg >0){print contig "\t" cov/lg;} contig = $1; cov = cov + ($2*$3); lg = $4}}' $COVERAGE > out_coverage.txt;


