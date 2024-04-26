#---------------------------#
#############################
##### files & functions #####
#############################
#---------------------------#
## Import modules
import os
import re
import importlib

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


## Filenames
# inputs
KRAKEN_LIB = check_slash(config["kraken_library"])
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

## get output files from checkpoints
# from rule extract_concordants
def get_concord_files(wildcards):
    checkpoint_output = checkpoints.extract_concordants.get(**wildcards).output[0]
    return expand("{fold}mapped/concord/{i}.txt", fold = OUTPUT_FOLD, i = glob_wildcards(os.path.join(checkpoint_output, "{i}.txt")).i)

# from rule sort_concordants
def get_concordR_files(wildcards):
    checkpoint_output = checkpoints.sort_concordants.get(**wildcards).output[0]
    return expand("{fold}mapped/concord_r/{i}.fastq", fold = OUTPUT_FOLD, i = glob_wildcards(os.path.join(checkpoint_output, "{i}.fastq")).i)

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

# from rule remapping
def get_svgfiles(wildcards):
    checkpoint_output = checkpoints.remapping.get(**wildcards).output[0]
    return expand("{fold}mapped/coverage/{i}.svg", fold = OUTPUT_FOLD, i = glob_wildcards(os.path.join(checkpoint_output, "{i}.svg")).i)
 


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
            fa_files = expand("{fold}genome_data/fa/{species}.fa", fold = OUTPUT_FOLD, species = ALL_SPECIES),
            dir = directory(OUTPUT_FOLD + "genome_data/kraken_kefir_library")



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
    rule bowtie_mapping:
        input:
            r1 = R1,
            r2 = R2,
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
            bowtie2 -p 10 -x {params.idx} -1 {input.r1} -2 {input.r2} -S {output.sam}
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
    rule sort_unmapped:
        input:
            r1 = R1,
            r2 = R2,
            unmap = expand("{fold}unmapped/{name}_{build}.unmapped.txt", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
        output:
            unmap_r1 = expand("{fold}unmapped/{name}_{build}.unmapped_R1.fasta", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME),
            unmap_r2 = expand("{fold}unmapped/{name}_{build}.unmapped_R2.fasta", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
        conda:
            CONFIG_FOLDER + "seqtk.yml"
        shell:
            """
            echo "sorting unmapped reads..."
            seqtk subseq {input.r1} {input.unmap} | seqtk seq -a - > {output.unmap_r1}
            seqtk subseq {input.r2} {input.unmap} | seqtk seq -a - > {output.unmap_r2}
            """

    # print number of unmapped reads in metrics file
    rule print_unmapped:
        input:
            unmap = expand("{fold}unmapped/{name}_{build}.unmapped.txt", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
        output:
            metrics = expand("{fold}unmapped/metrics.txt", fold = OUTPUT_FOLD)
        params:
            read_lg = READ_LG
        shell:
            """
            echo "calculating metrics of unmapped reads..."
            echo -e "library\tNumber mates\tNumber reads\tNumber mapped nt" > {output.metrics}
            # number of lines
            LG=$(wc -l {input.unmap} | awk -F " " '{{print $1}}')
            # write data
            echo -e "unmapped\t$(expr $LG / 2)\t$LG\t$(expr $LG \* {params.read_lg})" >> {output.metrics}
            """

    # run unmapped reads through Kraken2
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
            CONFIG_FOLDER + "kraken2.yml"
        shell:
            """
            echo "aligning unmapped reads..."
            kraken2 --db {params.kraken_library} --threads 8 --report {output.report} --report-minimizer-data --paired {input.unmap_r1} {input.unmap_r2} > {output.out}
            """


    ## Analysis per species
    # extract all species present in sam file and put them in txt file
    rule extract_species:
        input:
            sam = expand("{fold}{name}_{build}.sam", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
        output:
            species_file = temp(SPECIES_LST)
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
            sam = expand("{fold}{name}_{build}.sam", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME)
        output:
            directory(OUTPUT_FOLD + "mapped/concord/"),
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
                CONCORD={params.folder}mapped/concord/{params.name}_{params.build}_concSH_${{sp}}.txt

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
    checkpoint sort_concordants:
        input:
            get_concord_files,
            species_file = SPECIES_LST,
            r1 = R1,
            r2 = R2
        output:
            directory(OUTPUT_FOLD + "mapped/concord_r/")
        conda:
            CONFIG_FOLDER + "seqtk.yml"
        params: 
            name = R1_NAME, 
            build = BUILD_NAME,
            folder = OUTPUT_FOLD
        shell:
            """
            mkdir -p {output[0]}
            species=$(cat {input.species_file})
            echo "sorting concordants..."
            # for each species
            for sp in $species; do
                echo $sp

                # filenames
                CONCORD={params.folder}mapped/concord/{params.name}_{params.build}_concSH_${{sp}}.txt
                CONCORD_R1={params.folder}mapped/concord_r/{params.name}_{params.build}_concSH_${{sp}}_R1.fastq
                CONCORD_R2={params.folder}mapped/concord_r/{params.name}_{params.build}_concSH_${{sp}}_R2.fastq

                # divide concordants into R1 and R2
                seqtk subseq {input.r1} $CONCORD > $CONCORD_R1
                seqtk subseq {input.r2} $CONCORD > $CONCORD_R2
            done
            """

    # create metrics file for each species
    checkpoint metrics_per_species:
        input:
            species_file = SPECIES_LST,
            concord_r = get_concordR_files
        output:
            directory(OUTPUT_FOLD + "mapped/metrics/")
        params: 
            name = R1_NAME, 
            build = BUILD_NAME,
            folder = OUTPUT_FOLD,
            read_lg = READ_LG
        shell:
            """
            echo "calculating metrics per species..."
            mkdir -p {output[0]}
            species=$(cat {input.species_file})
            # for each species
            for sp in $species; do
                
                # filenames
                CONCORD_R1={params.folder}mapped/concord_r/{params.name}_{params.build}_concSH_${{sp}}_R1.fastq
                METRIC={params.folder}mapped/metrics/${{sp}}_metrics.txt

                echo -e "library\tNumber mates\tNumber reads\tNumber mapped nt" > $METRIC
                # number of lines
                LG=$(wc -l $CONCORD_R1 | awk -F " " '{{print $1}}')
                # number of reads
                READS=$(expr $LG / 2)
                # write data
                echo -e "$sp\t$(expr $LG / 2)\t$READS\t$(expr $READS \* {params.read_lg})" >> $METRIC
            done
            """

    # for each species: remapping on reference genome
    # create .bam files
    checkpoint setup_remapping:
        input:
            species_file = SPECIES_LST,
            sam = expand("{fold}{name}_{build}.sam", fold = OUTPUT_FOLD, name = R1_NAME, build = BUILD_NAME),
            concord = get_concord_files
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
        threads: 16
        shell:
            """
            mkdir -p {output[0]}
            species=$(cat {input.species_file})
            echo "setting up for remapping..."
            # for each species
            for sp in $species; do
                echo $sp

                # filenames
                CONCORD={params.folder}mapped/concord/{params.name}_{params.build}_concSH_${{sp}}.txt
                REMAP_BAM={params.folder}mapped/bam/{params.name}_{params.build}_concSH_${{sp}}.bam

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
            CONFIG_FOLDER + "jvarkit.yml"
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
                java -jar $PICARD CreateSequenceDictionary -VERBOSITY ERROR -QUIET true -R $GENOME -O $GENOME.dict 
                jvarkit wgscoverageplotter -C $MEDGRAPH -R $GENOME $REMAP_BAM -o $COVIMG
            done
            """

    rule global_metrics:
        input:
            r1 = R1,
            species_file = SPECIES_LST,
            metrics_unmapped = expand("{fold}unmapped/metrics.txt", fold = OUTPUT_FOLD),
            metrics_mapped = get_metrics_files
        output:
            glob_metrics = expand("{fold}global_metrics.txt", fold = OUTPUT_FOLD)
        params: 
            name = R1_NAME, 
            build = BUILD_NAME,
            folder = OUTPUT_FOLD,
            read_lg = READ_LG
        shell:
            """
            # print number of reads at the beginning of file
            MATES=$(zcat {input.r1} | echo $((`wc -l`/4)))
            echo $(expr $MATES \* 2) > {output.glob_metrics}
            # copy unmapped metrics file
            cat {input.metrics_unmapped} >> {output.glob_metrics}
            species=$(cat {input.species_file})
            # for each species
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
    rule build_kraken_library:
        input:
            fa_files = expand("{fold}genome_data/fa/{species}.fa", fold = OUTPUT_FOLD, species = ALL_SPECIES)
        output:
            directory(OUTPUT_FOLD + "genome_data/kraken_kefir_library")
        conda:
            CONFIG_FOLDER + "kraken2.yml"
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
            

