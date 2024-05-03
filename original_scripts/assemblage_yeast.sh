# DATA TO MODIFY
##complete path to fastq file: assembly input file
INPUT="/home/ajoubert/Assemblage_Yeast_JBB/FAW73023_demultiplexing_pass/NB01.fastq";
ILLUMINA_R1="";
ILLUMINA_R2="";
OUTPUT_FOLDER="/home/ajoubert/Assemblage_Yeast_JBB/Assemblage_NB01/";
#GENOME_SIZE="12.1m";
#nanopore kit
KIT="ont-ligation";
# complete path to reference genome(s) - concatenation of genomes if needed
# ref genome - S. cerevisiae
REF_GENOME_INIT="/home/ajoubert/Assemblage_01082023/Genome_Proteome_ref/GCF_000146045.2_R64_genomic.fna";
# proteome - S. ceverisiae
REF_PROTEOME_INIT="/home/ajoubert/Assemblage_01082023/Genome_Proteome_ref/GCF_000146045.2_R64_protein.faa";
# GC content
GC=38.3; ## check to keep it automat# DATA TO MODIFY
#reference species for augustus
SPECIES="saccharomyces_cerevisiae_S288C";


### SCRIPT DO NOT TOUCH
#----------------------
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
FLYE_ASSEMBLY=$PREF".contigs.fasta";
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
#~ # glimmer
#~ ICM=$PREF"_consensusRaconMedaka.contigs.icm";
#~ longorfs=$PREF"_consensusRaconMedaka.contigs.longorfs";
#~ train=$PREF"_consensusRaconMedaka.contigs.train";
#~ ORFS=$PREF"_consensusRaconMedaka.contigs.ORFs";
#~ PREDICT=$PREF"_consensusRaconMedaka.contigs.ORFs.predict";
# augustus
GFF_AUG=$PREF".contigs.augustus.gff";
GFF=$PREF".contigs.genes.gff";
ORFs_FA=$PREF".contigs.genes.fa";
PEP_FA=$PREF".contigs.pep.fa";
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

#~ export PATH=$PATH:/home/duvernois/TOOLS/augustus/bin:/home/duvernois/TOOLS/augustus/scripts;
AUGUSTUS="augustus";
BEDTOOLS="bedtools";
BEDFA=$BEDTOOLS" getfasta";
BEDCOV=$BEDTOOLS" genomecov";
AUG2PEP="/home/duvernois/TOOLS/UTILS/augustus2fasta.py";
#~ AUG2FA="/home/duvernois/TOOLS/UTILS/gff2seq.py";
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

#~ ## genome assembly
echo "FLYE assembly..."
$FLYE --nano-raw $CLEAN_ASSEMBLY --out-dir $FLYE_FOLDER --threads 8;
# --genome-size $GENOME_SIZE --asm-coverage 50;
mv $FLYE_FOLDER/assembly.fasta $FLYE_ASSEMBLY;

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

$BUSCO  -i $MEDAKA_ASSEMBLY -m genome -l bacteria -c 10 --out $BUSCO_FILE;
python3 $PLOT_BUSCO -wd $BUSCO_FILE;


## comparison of all genome
echo "comparison with reference genome..."
$nucmer --prefix=$PREF_NUCMER $REF_GENOME $MEDAKA_ASSEMBLY;
$nucmer --prefix=$PREFREV_NUCMER $MEDAKA_ASSEMBLY $REF_GENOME;
show-coords -rcl $PREF_NUCMER.delta > $PREF_NUCMER.coords
$MUMMERPLOT --postscript --prefix=$PREF_NUCMER $PREF_NUCMER.delta -Q $MEDAKA_ASSEMBLY -R $REF_GENOME --filter --layout

## to do only if one chromosome and not several on REF genome
#~ nb_chr=$(grep ">" $REF_GENOME | wc -l);
#~ if [ $nb_chr==1 ]
#~ then
bash $DIVRATE $PREF_NUCMER $PREFREV_NUCMER > out_divergence.txt;
#~ fi

## gene prediction
echo "gene prediction with augustus..."
$AUGUSTUS --species=$SPECIES $MEDAKA_ASSEMBLY > $GFF_AUG;
awk -v var="$PREF" -F "\t" '{if ($3 == "gene") print $1 "\t" $2 "\t" $3"\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\tID=" var"_gene_" $9}' $GFF_AUG > $GFF;
$BEDFA -s -fi $FLYE_ASSEMBLY -bed  $GFF_AUG -fo $ORFs_FA;
python3 $AUG2PEP -f $GFF_AUG -o $PEP_FA -n $PREF;


cat $REF_PROTEOME_INIT | sed 's/\r//'  > $REF_PROTEOME;
echo "Annotation running..."
$MBDB -in $REF_PROTEOME -dbtype prot -title $TITLEDB;
echo "\t\tblastp..."
$BLASTP -db $REF_PROTEOME -query $PEP_FA -outfmt  "7 qacc sacc evalue qstart qend sstart send qseqid stitle" -out $BLASTP_ANNOT -max_target_seqs 5 -evalue 0.01;


# create tmp gff file with only annotated data
cat $BLASTP_ANNOT | grep -v "#" | sed 's/\./,/' | awk -F "\t" '{split($9, a, "["); split(a[1], b ,"."); split(b[2], c, ","); print $1   "\t" c[1]}' |
        awk '{ gsub(/[ ]+$/,""); print }' | sed 's/1 //' | sort -u | sed 's/ /%%/g' |
        awk -F "\t" 'BEGIN{id=""; data = ""}{ if ($1 == id) {data = data "/"$2;}
                else { print id "\t" data; data = $2 ; id = $1;}}END{print id "\t" data}' |awk -F "\t" '{split($1,a,"_"); print a[1] "_" a[2] "\t" $0  }' | sort -k 2,2 > bl.tmp;

# all gene ID
# create empty gff
cat $GFF | awk -F "\t" '{split($9,a, "="); print a[2] "\t" $0 ";Name=unknownFunction" }' | sort -k 1,1 > gff.tmp;
#grep ">" ../$PEP_FA | sed 's/>//' | awk -F " " '{split($1, a, "_"); split($2,b, "="); split($3,c,"="); split($4,d,"=");
#       print $1 "\t" a[1] "\taugustus\tgene\t" b[2] "\t" c[2] "\t.\t" substr(d[2],1,1) "\t.\tID=" $1";Name=unknownFunction"}' | sort -k 1,1 > gff_init.tmp
# genes with annotation
awk -F "\t" '{print $2}' bl.tmp | sort -k 1,1 > ID_gffannot.tmp;

# cross blast results with previous GFF
join -j1 1 -j2 2 gff.tmp bl.tmp | awk -F " " '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t ID=" $1 ";Name="$12}' | sed 's/%%/ /g' > bl_annot.tmp


# cross annotation and initial list of genes
join -j1 1 -j2 1 gff.tmp ID_gffannot.tmp | sed 's/ /\t/g' > join.tmp;
diff join.tmp gff.tmp | grep ">" | sed 's/> //' |
    awk -F "\t" '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' > missing.tmp

cat bl_annot.tmp missing.tmp >$BLASTP_GFF;

# assembly in genbank
echo "create genbank file...";
seqret -sequence $OUTPUT_FOLDER/$FLYE_ASSEMBLY -osformat genbank -osname_outseq $ASSEMBLY_GB -ofdirectory_outseq gbk_file -auto;
# tool to merge assembly in GB with features data
python3 $MERGEGB -f $BLASTP_GFF -g $ASSEMBLY_GB.genbank -n barcode02 -o $BLASTP_GB.genbank;
rm *.tmp;

# rescue of unknow genes
grep "unknown" $BLASTP_GFF | awk -F "\t" '{split($9,a,";"); split(a[1], b, "="); print b[2]}' > ID_unknown.txt;
$seqtk subseq ../$PEP_FA ID_unknown.txt > ID_unknown.fa;



################################
#### GET COVERAGE OF MAPPING ###
################################
#1. remapping reads on assembly
$MINIMAP -a -x map-ont $MEDAKA_ASSEMBLY $CLEAN_ASSEMBLY  | samtools view -bS  - | samtools sort - > $READS_REMAP;

#2. create genome file
python $CREATE_GENOME -f $MEDAKA_ASSEMBLY  -o $GENOME_LG;

#3. get coverage file
genomeCoverageBed -ibam $READS_REMAP -g $GENOME > $COVERAGE;

# print coverage (dans le terminal - si je fais pas d'erreur de calcul)
awk -F "\t" 'BEGIN{contig = ""; cov = ""; lg = "";}{if ($1 == contig) {cov = cov + ($2*$3); lg = $4; } else { if (lg >0){print contig "\t" cov/lg;} contig = $1; cov = cov + ($2*$3); lg = $4}}' $COVERAGE > out_coverage.txt;


