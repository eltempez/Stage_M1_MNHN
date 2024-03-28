## Import modules
import os

## Configuration file
configfile: "/trinity/home/etempez/test_remapping/smk_files/config.yaml"

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
SAM = f"{R1_NAME}_{BUILD_NAME}.sam"
SORT = f"{R1_NAME}_{BUILD_NAME}_sort.sam"
BH = f"{R1_NAME}_{BUILD_NAME}_BH.sam"
UNMAP = f"{R1_NAME}_{BUILD_NAME}_unmapped.txt"
UNMAP_FA_R1 = f"{R1_NAME}_{BUILD_NAME}_unmapped_R1.fasta"
UNMAP_FA_R2 = f"{R1_NAME}_{BUILD_NAME}_unmapped_R2.fasta"
UNMAP_QUAST = f"QUAST_{R1_NAME}_{BUILD_NAME}_unmapped"

## snamake rules
rule all:
    input:
        SAM
   
# Mapping against metagenome
rule bowtie_mapping:
    input:
        r1 = R1,
        r2 = R2,
    params:
        idx = BT_BUILD
    output:
        sam = SAM
    envmodules:
        "userspace/tr17.10",
        "biology",
        "bowtie2"
    shell:
        "bowtie2 -p 10 -x {params.idx} -1 {input.r1} -2 {input.r2} -S {output.sam} 2> log_bowtie2"


