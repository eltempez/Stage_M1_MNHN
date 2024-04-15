#############################
##### files & functions #####
#############################
## Import modules
import os

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


## Filenames
# inputs
R1 = config["r1"]
R2 = config["r2"]
READ_LG = config["read_lg"]
SPECIES_FOLD = check_slash(config["folder_species"])
OUTPUT_FOLD = check_slash(config["output_folder"])
BUILD_NAME = config["build_name"]
R1_NAME = get_prefix(R1, True)
# outputs
SPECIES_LST = f"{OUTPUT_FOLD}species.txt"

## wildcards
# .fasta and .fa files in the initial metagenome
FASTA_SPECIES, = glob_wildcards(os.path.join(SPECIES_FOLD, "{fasta_species}.fasta"))
FA_SPECIES, = glob_wildcards(os.path.join(SPECIES_FOLD, "{fa_species}.fa"))
# all species in the metagenome
ALL_SPECIES = FASTA_SPECIES + FA_SPECIES

# function to get output files from rule extract_concordants
def get_concord_files(wildcards):
    checkpoint_output = checkpoints.extract_concordants.get(**wildcards).output[0]
    return expand("{fold}mapping/mapped/{i}.txt", fold = OUTPUT_FOLD, i = glob_wildcards(os.path.join(checkpoint_output, "{i}.txt")).i)

# function to get output files from 

    
###########################
##### snakemake rules #####
###########################
rule all:
    input:
        unmap_r1 = expand("{fold}mapping/unmapped/{name}_{build}.unmapped_R1.fasta", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME),
        unmap_r2 = expand("{fold}mapping/unmapped/{name}_{build}.unmapped_R2.fasta", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME),
        metrics = expand("{fold}mapping/unmapped/metrics.txt", fold = OUTPUT_FOLD),
        concord_r1 = expand("{fold}mapping/mapped/{name}_{build}_concSH_{{species}}_R1.fastq", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME),
        concord_r2 = expand("{fold}mapping/mapped/{name}_{build}_concSH_{{species}}_R2.fastq", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)


## setup
# in case of duplicated contig names in fasta files : add the species name at the beginning of the contig name
rule handle_contigs_fa:
    output:
        fa_files = expand("{fold}genome_data/fa/{species}.fa", fold = OUTPUT_FOLD, species = ALL_SPECIES)
    params:
        species_fold = SPECIES_FOLD,
        output_fold = expand("{out}genome_data/fa/", out = OUTPUT_FOLD)
    run:
        # create output directory if doesn't exist
        shell("mkdir -p {params.output_fold}")
        # for each file in the matagenome
        for file in os.listdir(params.species_fold):
            species_name = get_prefix(file)
            output_name = f"{params.output_fold[0]}{species_name}.fa"
            with open(f"{params.species_fold}{file}", "r") as f_in:
                with open(output_name, "w") as f_out:
                    for line in f_in:
                        if line.startswith(">") and species_name not in line:
                        # change the contig name if needed
                            new_line = f">{species_name}_{line[1:]}"
                            f_out.write(new_line)
                        else:
                            f_out.write(line)


# build bowtie library
rule build_library:
    input:
        fa_files = expand("{fold}genome_data/fa/{species}.fa", fold = OUTPUT_FOLD, species = ALL_SPECIES)
    output:
        bowtie_index = expand("{folder}genome_data/index/{index_name}{extension}.bt2", folder = OUTPUT_FOLD, index_name = BUILD_NAME, extension = [".1", ".2", ".3", ".4", ".rev.1", ".rev.2"])
    params:
        fa_fold = expand("{fold}genome_data/fa/", fold = OUTPUT_FOLD),
        output_fold = expand("{fold}genome_data/index/", fold = OUTPUT_FOLD),
        index_name = BUILD_NAME
    envmodules:
        "biology",
        "bowtie2"
    shell:
        """
        echo "building bowtie library..."
        bowtie2-build -q {params.fa_fold}* {params.output_fold}{params.index_name}
        """
    

## mapping against metagenome
rule bowtie_mapping:
    input:
        r1 = R1,
        r2 = R2,
        bowtie_index = expand("{folder}genome_data/index/{index_name}{extension}.bt2", folder = OUTPUT_FOLD, index_name = BUILD_NAME, extension = [".1", ".2", ".3", ".4", ".rev.1", ".rev.2"])
    params:
        idx = expand("{fold}genome_data/index/{name}", fold = OUTPUT_FOLD, name = BUILD_NAME)
    output:
        sam = expand("{fold}mapping/{name}_{build}.sam", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
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
        sam = expand("{fold}mapping/{name}_{build}.sam", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
    output:
        unmap = expand("{fold}mapping/unmapped/{name}_{build}.unmapped.txt", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
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
        unmap = expand("{fold}mapping/unmapped/{name}_{build}.unmapped.txt", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
    output:
        unmap_r1 = expand("{fold}mapping/unmapped/{name}_{build}.unmapped_R1.fasta", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME),
        unmap_r2 = expand("{fold}mapping/unmapped/{name}_{build}.unmapped_R2.fasta", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
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
        unmap = expand("{fold}mapping/unmapped/{name}_{build}.unmapped.txt", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
    output:
        metrics = expand("{fold}mapping/unmapped/metrics.txt", fold = OUTPUT_FOLD)
    shell:
        """
        echo "** Number of unmapped reads (mate - to *2)" >> {output.metrics}
        wc -l {input.unmap} >> {output.metrics}
        """

## Analysis per species
# extract all species present in sam file and put them in txt file
rule extract_species:
    input:
        sam = expand("{fold}mapping/{name}_{build}.sam", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
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
checkpoint extract_concordants:
    input: 
        species_file = SPECIES_LST,
        sam = expand("{fold}mapping/{name}_{build}.sam", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
    output:
        directory(OUTPUT_FOLD + "mapping/mapped/")
    envmodules:
        "biology",
        "samtools"
    params: 
        name = R1_NAME, 
        build = BUILD_NAME,
        folder = OUTPUT_FOLD
    shell:
        """
        mkdir -p {output[0]}
        species=$(cat {input.species_file})
        echo "looking for concordants..."
        # for each species
        for sp in $species; do
            echo $sp

            # filename
            CONCORD={params.folder}mapping/mapped/{params.name}_{params.build}_concSH_${{sp}}.txt

            # if species name present in line, prints line in temporary out file
            awk -F "\t" -v s="${{sp}}" '{{split($2,b,":"); split(b[2],c,"_"); split($3,a,"_");
            if ((a[1]==s) || ((b[1] == "SN") && (c[1] == s))) print $0}}' {input.sam} > {params.folder}out

            # reads concordants from out file
            samtools view -F 4 --verbosity 2 {params.folder}out | awk -F "\t" '{{print $1}}' | sort | uniq -c | awk -F " " '{{if ($1 == "2") print $2}}' > $CONCORD

            # remove temporary file
            rm {params.folder}out
        done
        """

# for each species : sort concordants in fastq files depending on R1 and R2 
rule sort_concordants:
    input:
        get_concord_files,
        species_file = SPECIES_LST,
        r1 = R1,
        r2 = R2
    output:
        concord_r1 = expand("{fold}mapping/mapped/{name}_{build}_concSH_{{species}}_R1.fastq", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME),
        concord_r2 = expand("{fold}mapping/mapped/{name}_{build}_concSH_{{species}}_R2.fastq", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
    conda:
        "seqtk.yml"
    params: 
        name = R1_NAME, 
        build = BUILD_NAME,
        folder = OUTPUT_FOLD + "mapping/mapped/"
    shell:
        """
        species=$(cat {input.species_file})
        echo "sorting concordants..."
        # for each species
        for sp in $species; do
            echo $sp

            # filenames
            CONCORD={params.folder}{params.name}_{params.build}_concSH_${{sp}}.txt
            CONCORD_R1={params.folder}{params.name}_{params.build}_concSH_${{sp}}_R1.fastq
            CONCORD_R2={params.folder}{params.name}_{params.build}_concSH_${{sp}}_R2.fastq

            # divide concordants into R1 and R2
            seqtk subseq {input.r1} $CONCORD > $CONCORD_R1
            seqtk subseq {input.r2} $CONCORD > $CONCORD_R2
        done
        """

# create metrics file for each species
rule metrics_per_species:
    input:
        species_file = SPECIES_LST,
        concord_r1 = expand("{fold}mapped/{name}_{build}_concSH_{{species}}_R1.fastq", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
    output:
        metrics = expand("{fold}mapped/{{species}}_metrics.txt", fold = OUTPUT_FOLD)
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
        concord = expand("{fold}mapped/{name}_{build}_concSH_{{species}}.txt", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
    output:
        remap_bam = expand("{fold}remapping/{name}_{build}_concSH_{{species}}.bam", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME),
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
        fa_files = expand("{folder}{{species}}.fa", folder =  SPECIES_FOLD),
        remap_bam = expand("{fold}remapping/{name}_{build}_concSH_{{species}}.bam", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
    output:
        covfile = expand("{fold}remapping/{name}_{build}_concSH_{{species}}.coverage", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME),
        covimg = expand("{fold}remapping/{name}_{build}_concSH_{{species}}.svg", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
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
