#############################
##### files & functions #####
#############################
## Import modules
import os
from os.path import isfile, join


## Functions
def get_prefix(file_path, is_paired_illumina = False):
    # filename without the path
    basename = os.path.basename(file_path)
    # filename without the extension
    prefix = basename.split(".")[0]
    # if paired illumina : delete last _R1 or _1
    if is_paired_illumina:
        last_underscore = prefix.rfind("_")
        prefix = prefix[:last_underscore]
    return prefix

# for all folder paths : if doesn't end with "/", add it at the end
def check_slash(folder_path):
    if folder_path[-1] == "/" or folder_path == "":
        return folder_path
    return folder_path + "/"

# check if the files in a folder are fasta files -> return the list of species for which they are
def check_fasta_files(folder_path):
    species_fasta = []
    # get all files from folder
    onlyfiles = [f for f in os.listdir(folder_path) if isfile(join(folder_path, f))]
    # is they are .fasta files, add the prefix into a list
    for f in onlyfiles:
        if f.split(".")[-1] == "fasta":
            species_fasta.append(get_prefix(f))
    return species_fasta

# get species from species file, in order to iterate on them
def get_species(species_file):
    with open(species_file, "r") as f:
        species = [line.strip() for line in f]
    return species


## Filenames
#inputs
R1 = config["r1"]
R2 = config["r2"]
READ_LG = config["read_lg"]
SPECIES_FOLD = check_slash(config["folder_species"])
BUILD_NAME = config["build_name"]
R1_NAME = get_prefix(R1, True)
OUTPUT_FOLD = check_slash(config["output_folder"])
# outputs
BT_BUILD = f"{OUTPUT_FOLD}genome_data/{BUILD_NAME}"
SPECIES_LST = f"{OUTPUT_FOLD}species.txt"



###########################
##### snakemake rules #####
###########################
rule all:
    input:
        fa_files = expand("{folder}{species}.fa", folder =  SPECIES_FOLD, species = check_fasta_files(SPECIES_FOLD)),
        bowtie_index = expand("{folder}genome_data/{index_name}{extension}.bt2", folder = OUTPUT_FOLD, index_name = BUILD_NAME, extension = [".1", ".2", ".3", ".4", ".rev.1", ".rev.2"]),
        unmap_r1 = f"{OUTPUT_FOLD}{R1_NAME}_{BUILD_NAME}_unmapped_R1.fasta",
        unmap_r2 = f"{OUTPUT_FOLD}{R1_NAME}_{BUILD_NAME}_unmapped_R2.fasta",
        metrics = f"{OUTPUT_FOLD}metrics.txt",
        metrics_species = expand("{fold}mapped/{species}_metrics.txt", fold = OUTPUT_FOLD, species = get_species(SPECIES_LST)),
        covfile = expand("{fold}remapping/{name}_{build}_concSH_{species}.coverage", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME, species = get_species(SPECIES_LST)),
        covimg = expand("{fold}remapping/{name}_{build}_concSH_{species}.svg", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME, species = get_species(SPECIES_LST))      


## setup
# transform every .fasta file in the species folder into .fa files (for easier access)
rule create_fa:
    input:
        fasta_files = expand("{folder}{species}.fasta", folder = SPECIES_FOLD, species = check_fasta_files(SPECIES_FOLD))
    output:
        fa_files = expand("{folder}{species}.fa", folder =  SPECIES_FOLD, species = check_fasta_files(SPECIES_FOLD))
    params: 
        species_folder = SPECIES_FOLD
    shell:
        """
        for file in {input.fasta_files}; do
            mv "$file" "${{file%.fasta}}.fa"
        done
        """

# build bowtie library
rule build_library:
    output:
        bowtie_index = expand("{folder}genome_data/{index_name}{extension}.bt2", folder = OUTPUT_FOLD, index_name = BUILD_NAME, extension = [".1", ".2", ".3", ".4", ".rev.1", ".rev.2"])
    params:
        species_fold = SPECIES_FOLD,
        index_name = BUILD_NAME,
        output_fold = OUTPUT_FOLD + "genome_data/"
    envmodules:
        "biology",
        "bowtie2"
    shell:
        """
        echo "building bowtie library..."
        bowtie2-build -q {params.species_fold}* {params.output_fold}{params.index_name}
        """
    

## mapping against metagenome
rule bowtie_mapping:
    input:
        r1 = R1,
        r2 = R2
    params:
        idx = BT_BUILD
    output:
        sam = f"{OUTPUT_FOLD}{R1_NAME}_{BUILD_NAME}.sam"
    envmodules:
        "biology",
        "bowtie2"
    shell:
        """
        echo "mapping against metagenome..."
        bowtie2 -p 10 -x {params.idx} -1 {input.r1} -2 {input.r2} -S {output.sam}
        """

## rescue of unmapped reads
# get all unmapped reads
rule get_unmapped:
    input:
        sam = f"{OUTPUT_FOLD}{R1_NAME}_{BUILD_NAME}.sam"
    output:
        unmap = f"{OUTPUT_FOLD}{R1_NAME}_{BUILD_NAME}_unmapped.txt"
    envmodules:
        "biology",
        "samtools"
    shell:
        """
        echo "rescue of unmapped reads..."
        samtools view -f 4 {input.sam} |  awk -F "\t" '{{print $1}}' | sort -u > {output.unmap}
        """

# sort depending on R1 and R2
rule sort_unmapped:
    input:
        r1 = R1,
        r2 = R2,
        unmap = f"{OUTPUT_FOLD}{R1_NAME}_{BUILD_NAME}_unmapped.txt"
    output:
        unmap_r1 = f"{OUTPUT_FOLD}{R1_NAME}_{BUILD_NAME}_unmapped_R1.fasta",
        unmap_r2 = f"{OUTPUT_FOLD}{R1_NAME}_{BUILD_NAME}_unmapped_R2.fasta"
    conda:
        "seqtk.yml"
    shell:
        """
        echo "sorting unmapped reads..."
        seqtk subseq {input.r1} {input.unmap} | seqtk seq -a - > {output.unmap_r1}
        seqtk subseq {input.r2} {input.unmap} | seqtk seq -a - > {output.unmap_r2}
        """

# print number of unmapped reads in metrics file
rule print_unmapped:
    input:
        unmap = f"{OUTPUT_FOLD}{R1_NAME}_{BUILD_NAME}_unmapped.txt"
    output:
        metrics = f"{OUTPUT_FOLD}metrics.txt"
    shell:
        """
        echo "** Number of unmapped reads (mate - to *2)" >> {output.metrics}
        wc -l {input.unmap} >> {output.metrics}
        """

## Analysis per species
# extract all species present in sam file and put them in txt file
rule extract_species:
    input:
        sam = f"{OUTPUT_FOLD}{R1_NAME}_{BUILD_NAME}.sam"
    output:
        species_file = SPECIES_LST
    envmodules:
        "biology",
        "samtools"
    shell: 
        """
        echo "extracting species..."
        samtools view -H {input.sam} | awk -F "\t" '{{split($2,a,":"); split(a[2],b,"_"); if (a[1] == "SN") print b[1]}}' | sort -u > {output.species_file}
        """

# for each species : read concordants and put all of their id in a txt file
rule extract_concordants:
    input: 
        species_file = SPECIES_LST,
        sam = f"{OUTPUT_FOLD}{R1_NAME}_{BUILD_NAME}.sam",
    output:
        concord = expand("{fold}mapped/{name}_{build}_concSH_{species}.txt", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME, species = get_species(SPECIES_LST))
    envmodules:
        "biology",
        "samtools"
    params: 
        name = R1_NAME, 
        build = BUILD_NAME,
        folder = OUTPUT_FOLD
    shell:
        """
        species=$(cat {input.species_file})
        echo "looking for concordants..."
        # for each species
        for sp in $species; do
            echo $sp

            # filename
            CONCORD={params.folder}mapped/{params.name}_{params.build}_concSH_${{sp}}.txt

            # if species name present in line, prints line in temporary out file
            awk -F "\t" -v s="${{sp}}" '{{split($2,b,":"); split(b[2],c,"_"); split($3,a,"_");
            if ((a[1]==s) || ((b[1] == "SN") && (c[1] == s))) print $0}}' {input.sam} > {params.folder}out

            # reads concordants
            samtools view -F 4 --verbosity 2 {params.folder}out | awk -F "\t" '{{print $1}}' | sort | uniq -c | awk -F " " '{{if ($1 == "2") print $2}}' > $CONCORD

            # remove temporary file
            rm {params.folder}out
        done
        """

# for each species : sort concordants in fastq files depending on R1 and R2 
rule sort_concordants:
    input:
        species_file = SPECIES_LST,
        r1 = R1,
        r2 = R2,
        concord = expand("{fold}mapped/{name}_{build}_concSH_{species}.txt", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME, species = get_species(SPECIES_LST))
    output:
        concord_r1 = expand("{fold}mapped/{name}_{build}_concSH_{species}_R1.fastq", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME, species = get_species(SPECIES_LST)),
        concord_r2 = expand("{fold}mapped/{name}_{build}_concSH_{species}_R2.fastq", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME, species = get_species(SPECIES_LST))
    conda:
        "seqtk.yml"
    params: 
        name = R1_NAME, 
        build = BUILD_NAME,
        folder = OUTPUT_FOLD
    shell:
        """
        species=$(cat {input.species_file})
        echo "sorting concordants..."
        # for each species
        for sp in $species; do
            echo $sp

            # filenames
            CONCORD={params.folder}mapped/{params.name}_{params.build}_concSH_${{sp}}.txt
            CONCORD_R1={params.folder}mapped/{params.name}_{params.build}_concSH_${{sp}}_R1.fastq
            CONCORD_R2={params.folder}mapped/{params.name}_{params.build}_concSH_${{sp}}_R2.fastq

            # divide concordants into R1 and R2
            seqtk subseq {input.r1} $CONCORD > $CONCORD_R1
            seqtk subseq {input.r2} $CONCORD > $CONCORD_R2
        done
        """

# create metrics file for each species
rule metrics_per_species:
    input:
        species_file = SPECIES_LST,
        concord_r1 = expand("{fold}mapped/{name}_{build}_concSH_{species}_R1.fastq", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME, species = get_species(SPECIES_LST))
    output:
        metrics = expand("{fold}mapped/{species}_metrics.txt", fold = OUTPUT_FOLD, species = get_species(SPECIES_LST))
    params: 
        name = R1_NAME, 
        build = BUILD_NAME,
        folder = OUTPUT_FOLD,
        read_lg = READ_LG
    shell:
        """
        species=$(cat {input.species_file})
        # for each species
        for sp in $species; do
            
            # filenames
            CONCORD_R1={params.folder}mapped/{params.name}_{params.build}_concSH_${{sp}}_R1.fastq
            METRIC={params.folder}mapped/${{sp}}_metrics.txt

            # add infos to metrics file
            echo -e "Library\t$sp" > $METRIC
            # number of lines in fastq file
            LG=$(wc -l $CONCORD_R1 | awk -F " " '{{print $1}}')
            # number of mates = number of lines in R1 fastq file / 4 (since one read takes 4 lines)
            echo -e "Number mates\t$(expr $LG / 4)" >> $METRIC
            # number of reads
            READS=$(expr $LG / 2)
            echo -e "Number reads\t$READS" >> $METRIC
            # number of nt
            echo -e "Number mapped nt\t$(expr $READS \* {params.read_lg})" >> $METRIC
        done
        """

# for each species: remapping on reference genome
# create .bam files
rule setup_remapping:
    input:
        species_file = SPECIES_LST,
        sam = f"{OUTPUT_FOLD}{R1_NAME}_{BUILD_NAME}.sam",
        concord = expand("{fold}mapped/{name}_{build}_concSH_{species}.txt", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME, species = get_species(SPECIES_LST))
    output:
        remap_bam = expand("{fold}remapping/{name}_{build}_concSH_{species}.bam", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME, species = get_species(SPECIES_LST)),
    params: 
        name = R1_NAME, 
        build = BUILD_NAME,
        folder = OUTPUT_FOLD
    envmodules:
        "userspace",
        "biology",
        "java-JDK-OpenJDK/11.0.9",
        "picard",
        "samtools"
    threads: 16
    shell:
        """
        species=$(cat {input.species_file})
        echo "setting up for remapping..."
        # for each species
        for sp in $species; do
            echo $sp

            # filenames
            CONCORD={params.folder}mapped/{params.name}_{params.build}_concSH_${{sp}}.txt
            REMAP_BAM={params.folder}remapping/{params.name}_{params.build}_concSH_${{sp}}.bam

            # filter original sam with concordant reads only for the current species (save in temporaty file)
            java -jar $PICARD FilterSamReads --VERBOSITY ERROR --QUIET true -I {input.sam} -O {params.folder}int.sam -FILTER includeReadList -READ_LIST_FILE $CONCORD
            # keep only headers with contigs regarding the current species (save in temporaty file)
            awk -F "\t" -v s="$sp" '{{split($2, a, ":"); split(a[2],b,"_"); if ((a[1] == "SN") && (b[1] == s)) print $0; else if (a[1] != "SN") print $0}}' {params.folder}int.sam > {params.folder}int2.sam
            # sort and save in bam file
            samtools view -b {params.folder}int2.sam | samtools sort > $REMAP_BAM;
            # delete temporary files
            rm {params.folder}int.sam {params.folder}int2.sam
            # index bam file
            samtools index $REMAP_BAM
        done
        """



# creation of genome reference data
rule remapping:
    input:
        species_file = SPECIES_LST,
        fa_files = expand("{folder}{species}.fa", folder =  SPECIES_FOLD, species = get_species(SPECIES_LST)),
        remap_bam = expand("{fold}remapping/{name}_{build}_concSH_{species}.bam", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME, species = get_species(SPECIES_LST))
    output:
        covfile = expand("{fold}remapping/{name}_{build}_concSH_{species}.coverage", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME, species = get_species(SPECIES_LST)),
        covimg = expand("{fold}remapping/{name}_{build}_concSH_{species}.svg", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME, species = get_species(SPECIES_LST))
    envmodules:
        "userspace",
        "biology",
        "java-JDK-OpenJDK/11.0.9",
        "picard",
        "samtools"
    conda:
        "jvarkit.yml"
    params: 
        name = R1_NAME, 
        build = BUILD_NAME,
        folder = OUTPUT_FOLD,
        species_folder = SPECIES_FOLD
    shell:
        """
        species=$(cat {input.species_file})
        echo "calculating coverage..."
        # for each species
        for sp in $species; do
            echo $sp

            # filenames
            REMAP_BAM={params.folder}remapping/{params.name}_{params.build}_concSH_${{sp}}.bam
            GENOME={params.species_folder}${{sp}}.fa
            COVFILE={params.folder}remapping/{params.name}_{params.build}_concSH_${{sp}}.coverage
            COVIMG={params.folder}remapping/{params.name}_{params.build}_concSH_${{sp}}.svg

            # create reference genome data (.fa.fai file)
            samtools faidx $GENOME
            # create .fa.bed file
            awk -v sp="$sp" 'BEGIN {{FS="\t"}} {{print sp "_" $1 FS "0" FS $2}}' $GENOME.fai > $GENOME.bed
            # calculate coverage statistics
            jvarkit bamstats04 -B $GENOME.bed $REMAP_BAM > $COVFILE
            # calculate median
            MEDCOV=$(cat $COVFILE | sed '1d' | awk -F "\t" '{{print $9}}'| sort -k 1n,1 | tail -n1| cut -d "." -f 1)
            # if median = 0, median = 5
            if  [ $MEDCOV -eq 0 ]
            then
            MEDCOV=5
            fi
            # calculate parameters
            RATIO=$((MEDCOV * 30 / 100))
            MEDGRAPH=$(($MEDCOV + $RATIO))
            # create graph
            java -jar $PICARD CreateSequenceDictionary --VERBOSITY ERROR --QUIET true -R $GENOME -O $GENOME.dict
            jvarkit wgscoverageplotter -C $MEDGRAPH -R $GENOME $REMAP_BAM -o $COVIMG
        done
        """
