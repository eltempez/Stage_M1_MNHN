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

## Filenames
#input
R1 = config["r1"]
R2 = config["r2"]
BT_BUILD = config["bt_build"]
# Use functions to create output filenames
R1_NAME = get_prefix(R1, True)
BUILD_NAME = get_prefix(BT_BUILD)
OUTPUT_FOLD = config["output_folder"]
SAM = f"{OUTPUT_FOLD}{R1_NAME}_{BUILD_NAME}.sam"
SORT = f"{OUTPUT_FOLD}{R1_NAME}_{BUILD_NAME}_sort.sam"
BH = f"{OUTPUT_FOLD}{R1_NAME}_{BUILD_NAME}_BH.sam"
UNMAP = f"{OUTPUT_FOLD}{R1_NAME}_{BUILD_NAME}_unmapped.txt"
UNMAP_FA_R1 = f"{OUTPUT_FOLD}{R1_NAME}_{BUILD_NAME}_unmapped_R1.fasta"
UNMAP_FA_R2 = f"{OUTPUT_FOLD}{R1_NAME}_{BUILD_NAME}_unmapped_R2.fasta"
UNMAP_QUAST = f"{OUTPUT_FOLD}QUAST_{R1_NAME}_{BUILD_NAME}_unmapped"
METRICS = f"{OUTPUT_FOLD}metrics.txt"


###########################
##### snakemake rules #####
###########################
rule all:
    input:
        METRICS
   
## mapping against metagenome
rule bowtie_mapping:
    input:
        r1 = R1,
        r2 = R2,
    params:
        idx = BT_BUILD
    output:
        sam = SAM,
        log_bowtie = OUTPUT_FOLD + "log_bowtie2"
    envmodules:
        "biology",
        "bowtie2"
    shell:
        "bowtie2 -p 10 -x {params.idx} -1 {input.r1} -2 {input.r2} -S {output.sam} 2> {output.log_bowtie}"

## rescue of unmapped reads
# get all unmapped reads
rule get_unmapped:
    input:
        sam = SAM
    output:
        unmap = UNMAP
    envmodules:
        "biology",
        "samtools"
    shell:
        "samtools view -f 4 {input.sam} |  awk -F \"\t\" '{{print $1}}' | sort -u > {output.unmap}"

# sort depending on R1 and R2
rule sort_unmapped_R1:
    input:
        r1 = R1,
        unmap = UNMAP
    output:
        unmap_r1 = UNMAP_FA_R1
    conda:
        "seqtk.yml"
    shell:
        "seqtk subseq {input.r1} {input.unmap} | seqtk seq -a - > {output.unmap_r1}"

rule sort_unmapped_R2:
    input:
        r2 = R2,
        unmap = UNMAP
    output:
        unmap_r2 = UNMAP_FA_R2
    conda:
        "seqtk.yml"
    shell:
        "seqtk subseq {input.r2} {input.unmap} | seqtk seq -a - > {output.unmap_r2}"

# print number of unmapped reads in metrics file
rule print_unmapped:
    input:
        unmap = UNMAP
    output:
        metrics = METRICS
    shell:
        """
        echo "** Number of unmapped reads (mate - to *2)" >> {output.metrics}
        wc -l {input.unmap} >> {output.metrics}
        """

## Analysis per species
    