#---------------------------#
#############################
##### files & functions #####
#############################
#---------------------------#
## Import modules
import os
import re

## Functions
# get the prefix of a file name 
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

# configuration folder
CONFIG_FOLDER = check_slash(config["config_folder"])

# get ncbi accession tax ids from metagenome_ncbi_id.txt file
def get_ncbi_id_metagenome():
    file_path = CONFIG_FOLDER + "metagenome_ncbi_id.txt"
    dict_ncbi = {}
    with open(file_path, "r") as f_id:
        for line in f_id:
                line_split = re.split(r'\t+', line.strip())
                dict_ncbi[line_split[0]] = int(line_split[1])
    return dict_ncbi

# format metagenome dict to use in bash
def bash_dict(dict):
    printed_dict = ""
    for key, val in dict.items():
        printed_dict += f"['{key}']={val} "
    return printed_dict[:-1]


## Filenames
# inputs
KRAKEN_LIB = check_slash(config["kraken_library"])
R1 = config["r1"]
R2 = config["r2"]
READ_LG = config["read_lg"]
SPECIES_FOLD = check_slash(config["folder_species"])
OUTPUT_FOLD = check_slash(config["output_folder"])
BUILD_NAME = config["build_name"]
# check if paired or unpaired
if R2 == "":
    is_paired = False
else:
    is_paired = True
R1_NAME = get_prefix(R1, is_paired)
# outputs
SPECIES_LST = f"{OUTPUT_FOLD}species.txt"

## wildcards
# .fasta and .fa files in the initial metagenome
FASTA_SPECIES, = glob_wildcards(os.path.join(SPECIES_FOLD, "{fasta_species}.fasta"))
FA_SPECIES, = glob_wildcards(os.path.join(SPECIES_FOLD, "{fa_species}.fa"))
# all species in the metagenome
ALL_SPECIES = FASTA_SPECIES + FA_SPECIES

## get output files from checkpoints
# from rule extract_concordants
def get_concord_files(wildcards):
    checkpoint_output = checkpoints.extract_concordants.get(**wildcards).output[0]
    return expand("{fold}mapped/aligned/{i}.txt", fold = OUTPUT_FOLD, i = glob_wildcards(os.path.join(checkpoint_output, "{i}.txt")).i)

# from rule extract_concordants
def get_concordseq_files(wildcards):
    checkpoint_output = checkpoints.extract_concordants.get(**wildcards).output[0]
    return expand("{fold}mapped/aligned/{i}.fastq", fold = OUTPUT_FOLD, i = glob_wildcards(os.path.join(checkpoint_output, "{i}.fastq")).i)

# from rule extract_aligned
def get_aligned_files(wildcards):
    checkpoint_output = checkpoints.extract_aligned.get(**wildcards).output[0]
    return expand("{fold}mapped/aligned/{i}.txt", fold = OUTPUT_FOLD, i = glob_wildcards(os.path.join(checkpoint_output, "{i}.txt")).i)

def get_alignedseq_files(wildcards):
    checkpoint_output = checkpoints.extract_aligned.get(**wildcards).output[0]
    return expand("{fold}mapped/aligned/{i}.fastq", fold = OUTPUT_FOLD, i = glob_wildcards(os.path.join(checkpoint_output, "{i}.fastq")).i)

# chose file based on paired or unpaired
def get_txt_from_mapping():
    if is_paired:
        return get_concord_files
    else:
        return get_aligned_files

def get_fastq_from_mapping():
    if is_paired:
        return get_concordseq_files
    else:
        return get_alignedseq_files


# from rule metrics_per_species
def get_metrics_files(wildcards):
    checkpoint_output = checkpoints.metrics_per_species.get(**wildcards).output[0]
    return expand("{fold}mapped/metrics/{i}.txt", fold = OUTPUT_FOLD, i = glob_wildcards(os.path.join(checkpoint_output, "{i}.txt")).i)

# from rule setup_remapping
def get_bam_files(wildcards):
    checkpoint_output = checkpoints.setup_remapping.get(**wildcards).output[0]
    return expand("{fold}mapped/bam/{i}.bam", fold = OUTPUT_FOLD, i = glob_wildcards(os.path.join(checkpoint_output, "{i}.bam")).i)

# from rule remapping
def get_covfiles(wildcards):
    checkpoint_output = checkpoints.remapping.get(**wildcards).output[0]
    return expand("{fold}mapped/coverage/{i}.coverage", fold = OUTPUT_FOLD, i = glob_wildcards(os.path.join(checkpoint_output, "{i}.coverage")).i)

def get_svgfiles(wildcards):
    checkpoint_output = checkpoints.remapping.get(**wildcards).output[0]
    return expand("{fold}mapped/coverage/{i}.svg", fold = OUTPUT_FOLD, i = glob_wildcards(os.path.join(checkpoint_output, "{i}.svg")).i)

# from rule build_kraken_library
def get_kraken_library(wildcards):
    checkpoint_output = checkpoints.build_kraken_library.get(**wildcards).output[0]
    return expand("{fold}genome_data/kraken_kefir_library/{i}", fold = OUTPUT_FOLD, i = glob_wildcards(os.path.join(checkpoint_output, "{i}")).i)
 


#-------------------------#
###########################
##### snakemake rules #####
###########################
#-------------------------#
if config["use_bowtie"]:
    rule all:
        input:
            report_kraken = expand("{fold}unmapped/kraken_rescue_report.k2report", fold = OUTPUT_FOLD),
            out_kraken = expand("{fold}unmapped/kraken_rescue_output.kraken2", fold = OUTPUT_FOLD),
            covfiles = get_covfiles,
            svgfiles = get_svgfiles,
            glob_metrics = expand("{fold}global_metrics.txt", fold = OUTPUT_FOLD)
elif not config["use_bowtie"]:
    rule all:
        input:
            metrics = expand("{fold}global_metrics.txt", fold = OUTPUT_FOLD)



## setup
# in case of duplicated contig names in fasta files : add the species name at the beginning of the contig name
# + for kraken : add ncbi accession number
rule handle_contigs_fa:
    output:
        fa_files = expand("{fold}genome_data/fa/{species}.fa", fold = OUTPUT_FOLD, species = ALL_SPECIES)
    params:
        species_fold = SPECIES_FOLD,
        output_fold = expand("{out}genome_data/fa/", out = OUTPUT_FOLD),
        dict_ncbi = get_ncbi_id_metagenome()
    run:
        # create output directory if doesn't exist
        shell("mkdir -p {params.output_fold}")
        print("checking fa/fasta files...")
        # for each file in the matagenome
        for file in os.listdir(params.species_fold):
            # don't include hidden files
            if not file.startswith("."):
                # test if is a fa/fasta file
                if not re.search(r"\.(fa|fasta)$", file):
                    print("ERROR : one or several files in the metagenome aren't fa/fasta files")
                    break
                else:
                    species_name = get_prefix(file)
                    output_name = f"{params.output_fold[0]}{species_name}.fa"
                    with open(f"{params.species_fold}{file}", "r") as f_in:
                        with open(output_name, "w") as f_out:
                            for line in f_in:
                                # add the tax id number
                                if line.startswith(">"):
                                    line = line.strip() + f"|kraken:taxid|{params.dict_ncbi[species_name]}\n"
                                # change the contig name if needed
                                if line.startswith(">") and species_name not in line:
                                    new_line = f">{species_name}_{line[1:]}"
                                    f_out.write(new_line)
                                else:
                                    f_out.write(line)

# clean illumina reads
if is_paired:
    rule fastq_cleaning:
        input: 
            r1_brut = R1,
            r2_brut = R2
        output:
            R1_clean = expand("{fold}genome_data/fq_clean/{name}_clean_R1.fq.gz", fold = OUTPUT_FOLD, name = R1_NAME),
            R2_clean = expand("{fold}genome_data/fq_clean/{name}_clean_R2.fq.gz", fold = OUTPUT_FOLD, name = R1_NAME)
        conda:
            CONFIG_FOLDER + "yml/fastp.yml"
        params:
            folder = OUTPUT_FOLD
        shell:
            """
            echo "cleaning illumina reads..."
            fastp -i {input.r1_brut} -I {input.r2_brut} -o {output.R1_clean} -O {output.R2_clean} -j {params.folder}genome_data/fq_clean/fastp.json -h {params.folder}genome_data/fq_clean/fastp.html
            """
else:
    rule fastq_cleaning:
        input: 
            r1_brut = R1
        output:
            R1_clean = expand("{fold}genome_data/fq_clean/{name}_clean_R1.fq.gz", fold = OUTPUT_FOLD, name = R1_NAME)
        conda:
            CONFIG_FOLDER + "yml/fastp.yml"
        params:
            folder = OUTPUT_FOLD
        shell:
            """
            echo "cleaning illumina reads..."
            fastp -i {input.r1_brut} -o {output.R1_clean} -j {params.folder}genome_data/fq_clean/fastp.json -h {params.folder}genome_data/fq_clean/fastp.html
            """
    




################################
########## IF BOWTIE2 ##########
################################
if config["use_bowtie"]:
    # build bowtie library
    rule build_bowtie_library:
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
            bowtie2-build -q {params.fa_fold}*.fa {params.output_fold}{params.index_name}
            """
        

    ## mapping against metagenome
    if is_paired:
        rule bowtie_mapping:
            input:
                r1 = expand("{fold}genome_data/fq_clean/{name}_clean_R1.fq.gz", fold = OUTPUT_FOLD, name = R1_NAME),
                r2 = expand("{fold}genome_data/fq_clean/{name}_clean_R2.fq.gz", fold = OUTPUT_FOLD, name = R1_NAME),
                bowtie_index = expand("{folder}genome_data/index/{index_name}{extension}.bt2", folder = OUTPUT_FOLD, index_name = BUILD_NAME, extension = [".1", ".2", ".3", ".4", ".rev.1", ".rev.2"])
            params:
                idx = expand("{fold}genome_data/index/{name}", fold = OUTPUT_FOLD, name = BUILD_NAME)
            output:
                sam = expand("{fold}{name}_{build}.sam", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
            envmodules:
                "biology",
                "bowtie2"
            shell:
                """
                echo 'mapping against metagenome...'
                bowtie2 -p 10 -x {params.idx} -1 {input.r1} -2 {input.r2} -S {output.sam}
                """
    else:
        rule bowtie_mapping:
            input:
                r1 = expand("{fold}genome_data/fq_clean/{name}_clean_R1.fq.gz", fold = OUTPUT_FOLD, name = R1_NAME),
                bowtie_index = expand("{folder}genome_data/index/{index_name}{extension}.bt2", folder = OUTPUT_FOLD, index_name = BUILD_NAME, extension = [".1", ".2", ".3", ".4", ".rev.1", ".rev.2"])
            params:
                idx = expand("{fold}genome_data/index/{name}", fold = OUTPUT_FOLD, name = BUILD_NAME)
            output:
                sam = expand("{fold}{name}_{build}.sam", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
            envmodules:
                "biology",
                "bowtie2"
            shell:
                """
                echo "mapping against metagenome..."
                bowtie2 -p 10 -x {params.idx} -U {input.r1} -S {output.sam}
                """


    ## Analysis per species
    # extract all species present in sam file and put them in txt file
    rule extract_species:
        input:
            sam = expand("{fold}{name}_{build}.sam", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
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
    
    if is_paired:
        # for each species : read concordants and put all of their id in a txt file, then sort concordants in fastq files depending on R1 and R2 
        checkpoint extract_concordants:
            input: 
                r1 = expand("{fold}genome_data/fq_clean/{name}_clean_R1.fq.gz", fold = OUTPUT_FOLD, name = R1_NAME),
                r2 = expand("{fold}genome_data/fq_clean/{name}_clean_R2.fq.gz", fold = OUTPUT_FOLD, name = R1_NAME),
                species_file = SPECIES_LST,
                sam = expand("{fold}{name}_{build}.sam", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
            output:
                directory(OUTPUT_FOLD + "mapped/aligned/"),
            envmodules:
                "biology",
                "samtools"
            conda:
                CONFIG_FOLDER + "yml/seqtk.yml"
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
                    CONCORD={params.folder}mapped/aligned/{params.name}_{params.build}_concSH_${{sp}}.txt
                    CONCORD_R1={params.folder}mapped/aligned/{params.name}_{params.build}_concSH_${{sp}}_R1.fastq
                    CONCORD_R2={params.folder}mapped/aligned/{params.name}_{params.build}_concSH_${{sp}}_R2.fastq

                    # if species name present in line, prints line in temporary out file
                    awk -F "\t" -v s="${{sp}}" '{{split($2,b,":"); split(b[2],c,"_"); split($3,a,"_");
                    if ((a[1]==s) || ((b[1] == "SN") && (c[1] == s))) print $0}}' {input.sam} > {params.folder}out

                    # reads concordants from out file
                    samtools view -F 4 --verbosity 2 {params.folder}out | awk -F "\t" '{{print $1}}' | sort | uniq -c | awk -F " " '{{if ($1 == "2") print $2}}' > $CONCORD

                    # remove temporary file
                    rm {params.folder}out

                    # divide concordants into R1 and R2
                    seqtk subseq {input.r1} $CONCORD > $CONCORD_R1
                    seqtk subseq {input.r2} $CONCORD > $CONCORD_R2
                done
                """

    # for each species : read aligned and put all of their id in a txt file, then sort aligned reads in fastq file
    else:
        checkpoint extract_aligned:
            input:
                species_file = SPECIES_LST,
                sam = expand("{fold}{name}_{build}.sam", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME),
                r1 = expand("{fold}genome_data/fq_clean/{name}_clean_R1.fq.gz", fold = OUTPUT_FOLD, name = R1_NAME)
            output:
                directory(OUTPUT_FOLD + "mapped/aligned/")
            params: 
                name = R1_NAME, 
                build = BUILD_NAME,
                folder = OUTPUT_FOLD
            envmodules:
                "biology",
                "samtools"
            conda:
                CONFIG_FOLDER + "yml/seqtk.yml"
            shell:
                """
                mkdir -p {output[0]}
                species=$(cat {input.species_file})
                echo "sorting aligned reads..."
                # for each species
                for sp in $species; do
                    echo $sp

                    # filenames
                    ALIGNED={params.folder}mapped/aligned/{params.name}_{params.build}_concSH_${{sp}}.txt
                    ALIGNED_SEQ={params.folder}mapped/aligned/{params.name}_{params.build}_concSH_${{sp}}_R1.fastq

                    # if species name present in line, prints line in temporary out file
                    awk -F "\t" -v s="${{sp}}" '{{split($2,b,":"); split(b[2],c,"_"); split($3,a,"_");
                    if ((a[1]==s) || ((b[1] == "SN") && (c[1] == s))) print $0}}' {input.sam} > {params.folder}out

                    # reads all reads from temporary out file
                    samtools view -F 4 --verbosity 2 {params.folder}out | awk -F "\t" '{{print $1}}' | sort > $ALIGNED
                    # remove temporary file
                    rm {params.folder}out

                    # extract sequences
                    seqtk subseq {input.r1} $ALIGNED > $ALIGNED_SEQ
                done
                """

    # create metrics file for each species
    checkpoint metrics_per_species:
        input:
            species_file = SPECIES_LST,
            fastq = get_fastq_from_mapping(),
            r1_brut = R1
        output:
            directory(OUTPUT_FOLD + "mapped/metrics/")
        params: 
            name = R1_NAME, 
            build = BUILD_NAME,
            folder = OUTPUT_FOLD,
            read_lg = READ_LG,
            is_paired = int(is_paired)
        shell:
            """
            echo "calculating metrics per species..."
            mkdir -p {output[0]}
            species=$(cat {input.species_file})

            # total number of reads in whole sample
            NB_TOTAL_READS_R1=$(zcat {input.r1_brut} | echo $((`wc -l`/4)))
            if [ {params.is_paired} -eq 1 ]; then
                NB_TOTAL_READS=$((NB_TOTAL_READS_R1 * 2))
            else
                NB_TOTAL_READS=$NB_TOTAL_READS_R1
            fi

            # for each species
            for sp in $species; do
                    
                # filenames
                CONCORD_R1={params.folder}mapped/aligned/{params.name}_{params.build}_concSH_${{sp}}_R1.fastq
                METRIC={params.folder}mapped/metrics/${{sp}}_metrics.txt

                # number of reads for the species
                NB_READS_R1=$(($(wc -l < $CONCORD_R1) / 4))
                if [ {params.is_paired} -eq 1 ]; then
                    NB_READS=$((NB_READS_R1 * 2))
                else
                    NB_READS=$NB_READS_R1
                fi

                # print in file
                echo -e "library\tNumber reads\tNumber mapped nt\tPercents reads" > $METRIC
                echo -e "$sp\t$NB_READS\t$(($NB_READS * {params.read_lg}))\t$(echo "scale=2; $NB_READS * 100 / $NB_TOTAL_READS" | bc)" >> $METRIC

            done
                """

    # for each species: remapping on reference genome
    # create .bam files
    checkpoint setup_remapping:
        input:
            species_file = SPECIES_LST,
            sam = expand("{fold}{name}_{build}.sam", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME),
            concord = get_txt_from_mapping()
        output:
            directory(OUTPUT_FOLD + "mapped/bam/")
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
        shell:
            """
            mkdir -p {output[0]}
            species=$(cat {input.species_file})
            echo "setting up for remapping..."
            # for each species
            for sp in $species; do
                echo $sp

                # filenames
                CONCORD={params.folder}mapped/aligned/{params.name}_{params.build}_concSH_${{sp}}.txt
                REMAP_BAM={params.folder}mapped/bam/{params.name}_{params.build}_concSH_${{sp}}.bam

                # filter original sam with concordant reads only for the current species (save in temporaty file)
                java -jar $PICARD FilterSamReads --VERBOSITY ERROR --QUIET true -I {input.sam} -O {params.folder}int.sam -FILTER includeReadList -READ_LIST_FILE $CONCORD 2> {params.folder}log.txt
                
                # keep only headers with contigs regarding the current species (save in temporaty file)
                awk -F "\t" -v s="$sp" '{{split($2, a, ":"); split(a[2],b,"_"); if ((a[1] == "SN") && (b[1] == s)) print $0; else if (a[1] != "SN") print $0}}' {params.folder}int.sam > {params.folder}int2.sam
                # sort and save in bam file
                samtools view -b {params.folder}int2.sam | samtools sort > $REMAP_BAM;
                # delete temporary files
                rm {params.folder}int.sam {params.folder}int2.sam
                # index bam file
                samtools index $REMAP_BAM

                rm {params.folder}log.txt
            done
            """

    # creation of genome reference data + coverage data
    checkpoint remapping:
        input:
            species_file = SPECIES_LST,
            fa_files = expand("{fold}genome_data/fa/{species}.fa", fold = OUTPUT_FOLD, species = ALL_SPECIES),
            remap_bam = get_bam_files
        output:
            directory(OUTPUT_FOLD + "mapped/coverage/")
        envmodules:
            "userspace",
            "biology",
            "java-JDK-OpenJDK/11.0.9",
            "picard",
            "samtools"
        conda:
            CONFIG_FOLDER + "yml/jvarkit.yml"
        params: 
            name = R1_NAME, 
            build = BUILD_NAME,
            folder = OUTPUT_FOLD,
            species_folder = OUTPUT_FOLD + "genome_data/fa/"
        shell:
            """
            mkdir -p {output[0]}
            species=$(cat {input.species_file})
            echo "calculating coverage..."
            # for each species
            for sp in $species; do
                echo $sp

                # filenames
                REMAP_BAM={params.folder}mapped/bam/{params.name}_{params.build}_concSH_${{sp}}.bam
                GENOME={params.species_folder}${{sp}}.fa
                COVFILE={params.folder}mapped/coverage/{params.name}_{params.build}_concSH_${{sp}}.coverage
                COVIMG={params.folder}mapped/coverage/{params.name}_{params.build}_concSH_${{sp}}.svg

                # create reference genome data (.fa.fai file)
                samtools faidx $GENOME
                # create .fa.bed file
                awk 'BEGIN {{FS="\t"}} {{print $1 FS "0" FS $2}}' $GENOME.fai > $GENOME.bed
                # calculate coverage statistics
                jvarkit bamstats04 -B $GENOME.bed $REMAP_BAM > $COVFILE 2> {params.folder}log.txt
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
                java -jar $PICARD CreateSequenceDictionary -VERBOSITY ERROR -QUIET true -R $GENOME -O $GENOME.dict 2> {params.folder}log.txt
                jvarkit wgscoverageplotter -C $MEDGRAPH -R $GENOME $REMAP_BAM -o $COVIMG 2> {params.folder}log.txt

                rm {params.folder}log.txt
            done
            """


    ## rescue of unmapped reads
    # get all unmapped reads
    rule get_unmapped:
        input:
            sam = expand("{fold}{name}_{build}.sam", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
        output:
            unmap = expand("{fold}unmapped/{name}_{build}.unmapped.txt", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
        envmodules:
            "biology",
            "samtools"
        shell:
            """
            echo "rescue of unmapped reads..."
            samtools view -f 4 {input.sam} |  awk -F "\t" '{{print $1}}' | sort -u > {output.unmap}
            """

    # sort depending on R1 and R2
    if is_paired:
        rule sort_unmapped:
            input:
                r1 = expand("{fold}genome_data/fq_clean/{name}_clean_R1.fq.gz", fold = OUTPUT_FOLD, name = R1_NAME),
                r2 = expand("{fold}genome_data/fq_clean/{name}_clean_R2.fq.gz", fold = OUTPUT_FOLD, name = R1_NAME),
                unmap = expand("{fold}unmapped/{name}_{build}.unmapped.txt", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
            output:
                unmap_r1 = expand("{fold}unmapped/{name}_{build}.unmapped_R1.fasta", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME),
                unmap_r2 = expand("{fold}unmapped/{name}_{build}.unmapped_R2.fasta", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
            conda:
                CONFIG_FOLDER + "yml/seqtk.yml"
            shell:
                """
                echo "sorting unmapped reads..."
                seqtk subseq {input.r1} {input.unmap} | seqtk seq -a - > {output.unmap_r1}
                seqtk subseq {input.r2} {input.unmap} | seqtk seq -a - > {output.unmap_r2}
                """
    else:
        rule sort_unmapped:
            input:
                r1 = expand("{fold}genome_data/fq_clean/{name}_clean_R1.fq.gz", fold = OUTPUT_FOLD, name = R1_NAME),
                unmap = expand("{fold}unmapped/{name}_{build}.unmapped.txt", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
            output:
                unmap_r1 = expand("{fold}unmapped/{name}_{build}.unmapped_R1.fasta", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
            conda:
                CONFIG_FOLDER + "yml/seqtk.yml"
            shell:
                """
                echo "sorting unmapped reads..."
                seqtk subseq {input.r1} {input.unmap} | seqtk seq -a - > {output.unmap_r1}
                """

    # run unmapped reads through Kraken2
    if is_paired:
        rule kraken_unmapped:
            input:
                unmap_r1 = expand("{fold}unmapped/{name}_{build}.unmapped_R1.fasta", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME),
                unmap_r2 = expand("{fold}unmapped/{name}_{build}.unmapped_R2.fasta", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
            output:
                report = expand("{fold}unmapped/kraken_rescue_report.k2report", fold = OUTPUT_FOLD),
                out = expand("{fold}unmapped/kraken_rescue_output.kraken2", fold = OUTPUT_FOLD)
            params:
                kraken_library = KRAKEN_LIB
            conda:
                CONFIG_FOLDER + "yml/kraken2.yml"
            shell:
                """
                echo "aligning unmapped reads..."
                kraken2 --db {params.kraken_library} --threads 8 --report {output.report} --report-minimizer-data --paired {input.unmap_r1} {input.unmap_r2} > {output.out}
                """
    else:
        rule kraken_unmapped:
            input:
                unmap_r1 = expand("{fold}unmapped/{name}_{build}.unmapped_R1.fasta", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
            output:
                report = expand("{fold}unmapped/kraken_rescue_report.k2report", fold = OUTPUT_FOLD),
                out = expand("{fold}unmapped/kraken_rescue_output.kraken2", fold = OUTPUT_FOLD)
            params:
                kraken_library = KRAKEN_LIB
            conda:
                CONFIG_FOLDER + "yml/kraken2.yml"
            shell:
                """
                echo "aligning unmapped reads..."
                kraken2 --db {params.kraken_library} --threads 8 --report {output.report} --report-minimizer-data {input.unmap_r1} > {output.out}
                """
                

    ## global metrics
    rule global_metrics:
        input:
            r1_brut = R1,
            r1_clean = expand("{fold}genome_data/fq_clean/{name}_clean_R1.fq.gz", fold = OUTPUT_FOLD, name = R1_NAME),
            species_file = SPECIES_LST,
            metrics_mapped = get_metrics_files
        output:
            glob_metrics = expand("{fold}global_metrics.txt", fold = OUTPUT_FOLD)
        params: 
            name = R1_NAME, 
            build = BUILD_NAME,
            folder = OUTPUT_FOLD,
            read_lg = READ_LG,
            is_paired = int(is_paired)
        shell:
            """
            # total number of reads in whole sample
            NB_TOTAL_READS_R1=$(zcat {input.r1_brut} | echo $((`wc -l`/4)))
            if [ {params.is_paired} -eq 1 ]; then
                NB_TOTAL_READS=$((NB_TOTAL_READS_R1 * 2))
            else
                NB_TOTAL_READS=$NB_TOTAL_READS_R1
            fi

            NB_TOTAL_READS_R1_CLEAN=$(zcat {input.r1_clean} | echo $((`wc -l`/4)))
            if [ {params.is_paired} -eq 1 ]; then
                NB_TOTAL_READS_CLEAN=$((NB_TOTAL_READS_R1_CLEAN * 2))
            else
                NB_TOTAL_READS_CLEAN=$NB_TOTAL_READS_R1_CLEAN
            fi

            # print number of reads at the beginning of file
            echo -e "total reads\t$NB_TOTAL_READS" > {output.glob_metrics}
            echo -e "clean reads\t$NB_TOTAL_READS_CLEAN" >> {output.glob_metrics}
            
            echo -e "library\tNumber reads\tNumber mapped nt\tPercents reads" >> {output.glob_metrics}

            # for each species
            species=$(cat {input.species_file})
            for sp in $species; do
                # filename
                METRIC={params.folder}mapped/metrics/${{sp}}_metrics.txt
                # get the second line
                data=$(sed -n '2p' "$METRIC")
                # add it to the global file
                echo -e "$data" >> {output.glob_metrics}
            done
            """
    
    

################################
########## IF KRAKEN2 ##########
################################
elif not config["use_bowtie"]:
    # build kraken2 library
    checkpoint build_kraken_library:
        input:
            fa_files = expand("{fold}genome_data/fa/{species}.fa", fold = OUTPUT_FOLD, species = ALL_SPECIES)
        output:
            directory(OUTPUT_FOLD + "genome_data/kraken_kefir_library")
        conda:
            CONFIG_FOLDER + "yml/kraken2.yml"
        shell:
            """
            echo "building kraken2 library..."
            # add taxonomy data
            kraken2-build --download-taxonomy --db {output}
            # add fasta files to library
            for file in {input.fa_files}; do kraken2-build --add-to-library $file --db {output};  done
            # create library
            kraken2-build --build --db {output}
            """

    # align reads with kraken library + rescue unmapped reads
    if is_paired:
        rule kraken_alignment:
            input:
                r1 = expand("{fold}genome_data/fq_clean/{name}_clean_R1.fq.gz", fold = OUTPUT_FOLD, name = R1_NAME),
                r2 = expand("{fold}genome_data/fq_clean/{name}_clean_R2.fq.gz", fold = OUTPUT_FOLD, name = R1_NAME),
                kefir_library = get_kraken_library
            conda:
                CONFIG_FOLDER + "yml/kraken2.yml"
            output:
                report = expand("{fold}mapped/kraken_report.k2report", fold = OUTPUT_FOLD),
                out = expand("{fold}mapped/kraken_output.kraken2", fold = OUTPUT_FOLD),
                unmap_r1 = expand("{fold}unmapped/{name}_{build}_unmapped_1.fastq", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME),
                unmap_r2 = expand("{fold}unmapped/{name}_{build}_unmapped_2.fastq", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
            params:
                kefir_library = directory(OUTPUT_FOLD + "genome_data/kraken_kefir_library"),
                name = R1_NAME, 
                build = BUILD_NAME,
                folder = OUTPUT_FOLD
            shell:
                """
                UNMAP={params.folder}unmapped/{params.name}_{params.build}_unmapped#.fastq
                echo "aligning with kraken2..."
                kraken2 --db {params.kefir_library} --threads 8 --report {output.report} --report-minimizer-data --paired --unclassified-out $UNMAP {input.r1} {input.r2} > {output.out}
                """
    
    else:
        rule kraken_alignment:
            input:
                r1 = expand("{fold}genome_data/fq_clean/{name}_clean_R1.fq.gz", fold = OUTPUT_FOLD, name = R1_NAME),
                kefir_library = get_kraken_library
            params:
                kefir_library = directory(OUTPUT_FOLD + "genome_data/kraken_kefir_library"),
                name = R1_NAME, 
                build = BUILD_NAME,
                folder = OUTPUT_FOLD
            conda:
                CONFIG_FOLDER + "yml/kraken2.yml"
            output:
                report = expand("{fold}mapped/kraken_report.k2report", fold = OUTPUT_FOLD),
                out = expand("{fold}mapped/kraken_output.kraken2", fold = OUTPUT_FOLD),
                unmap_r1 = expand("{fold}unmapped/{name}_{build}_unmapped.fastq", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
            shell:
                """
                echo "aligning with kraken2..."
                kraken2 --db {params.kefir_library} --threads 8 --report {output.report} --report-minimizer-data --unclassified-out {params.folder}unmapped/{params.name}_{params.build}_unmapped.fastq {input.r1} > {output.out}
                """

    ## Analysis per species
    # global metrics per species
    rule global_metrics:
        input:
            r1_brut = R1,
            r1_clean = expand("{fold}genome_data/fq_clean/{name}_clean_R1.fq.gz", fold = OUTPUT_FOLD, name = R1_NAME),
            report = expand("{fold}mapped/kraken_report.k2report", fold = OUTPUT_FOLD)
        output:
            glob_metrics = expand("{fold}global_metrics.txt", fold = OUTPUT_FOLD)
        params:
            read_lg = READ_LG,
            tax_dict = bash_dict(get_ncbi_id_metagenome()),
            is_paired = int(is_paired)
        shell:
            """
            declare -A dict=({params.tax_dict})
            echo "calculating metrics..."

            # total number of reads in whole sample
            NB_TOTAL_READS_R1=$(zcat {input.r1_brut} | echo $((`wc -l`/4)))
            if [ {params.is_paired} -eq 1 ]; then
                NB_TOTAL_READS=$((NB_TOTAL_READS_R1 * 2))
            else
                NB_TOTAL_READS=$NB_TOTAL_READS_R1
            fi

            NB_TOTAL_READS_R1_CLEAN=$(zcat {input.r1_clean} | echo $((`wc -l`/4)))
            if [ {params.is_paired} -eq 1 ]; then
                NB_TOTAL_READS_CLEAN=$((NB_TOTAL_READS_R1_CLEAN * 2))
            else
                NB_TOTAL_READS_CLEAN=$NB_TOTAL_READS_R1_CLEAN
            fi

            # print number of reads at the beginning of file
            echo -e "total reads\t$NB_TOTAL_READS" > {output.glob_metrics}
            echo -e "clean reads\t$NB_TOTAL_READS_CLEAN" >> {output.glob_metrics}
            
            echo -e "library\tNumber reads\tNumber mapped nt\tPercents reads" >> {output.glob_metrics}

            # for each species
            for sp in "${{!dict[@]}}"; do
                id="${{dict[$sp]}}"

                # get line corresponding to the species
                ligne=$(awk -F '\t' -v id="$id" '{{if ($7 == id) print}}' {input.report})
                if [ ! -z "$ligne" ]; then
                    # extract right column
                    percent=$(echo "$ligne" | cut -f1)
                    nb_reads=$(echo "$ligne" | cut -f2)
                    echo -e "$sp\t$nb_reads\t$(($nb_reads * {params.read_lg}))\t$$(echo "scale=2; $nb_reads * 100 / $NB_TOTAL_READS" | bc)" >> {output.glob_metrics}
                fi
            done    

            """
            







