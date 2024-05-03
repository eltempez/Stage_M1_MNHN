import os

# get prefix from file name 
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

## Variables
# list of taxonomy ids
dict_ncbi = {"B3": 483199, "B4": 318683, "B6": 38307, "B1": 442, "B5": 28448, "B2": 436, "B11": 1597, "B18": 1588, "B10": 152331, "B7": 1245, "B14": 399370, "B13": 82688, "B17": 259059, "B21": 2203724, "B12": 304207, "BB": 5007, "Y3": 4926, "Y1": 4932, "Y2": 48255, "B16": 542, "B15": 336988}
# paths
bowtie_folder = "/home/eliotttempez/Documents/tests/outputs_mcf10_snk/mapped/aligned/"
kraken_file = "/home/eliotttempez/Documents/tests/kraken_kmers/kraken_output_40_kmers.txt"


# build dicts with read name as key, and tax id as value
global_dict = {}
# fill with kraken values
with open(kraken_file, "r") as kfile:
    for line in kfile:
        line_split = line.strip().split("\t")
        if line_split[0] == "C":
            global_dict[line_split[1]] = {}
            global_dict[line_split[1]]["kraken"] = int(line_split[2])
# fill with bowtie values
for file_name in os.listdir(bowtie_folder):
    species_name = get_prefix(file_name).split("_")[-1]
    file_path = os.path.join(bowtie_folder, file_name)
    if file_path.split(".")[-1] == "txt":
        with open(file_path, "r") as bfile:
            for line in bfile:
                if line.strip() not in global_dict:
                    global_dict[line.strip()] = {}
                global_dict[line.strip()]["bowtie"] = dict_ncbi[species_name]




# get mismatchs
# list of mismatched reads 
mismatched_reads = []
nb_matched = 0
for read in global_dict:
    if ("bowtie" not in global_dict[read]) or ("kraken" not in global_dict[read]) or (global_dict[read]["bowtie"] != global_dict[read]["kraken"]):
        mismatched_reads.append(read)
    else:
        nb_matched +=1

# output in file
with open("./mismatchs_bt_kr.txt", "w") as f_out:
    f_out.write(f"total concordant reads\t{nb_matched}\n")
    f_out.write(f"seq id\tbowtie\tkraken\n")
    for mis in mismatched_reads:
        if "bowtie" not in global_dict[mis]:
            f_out.write(f"{mis}\t0\t{global_dict[mis]['kraken']}\n")
        elif "kraken" not in global_dict[mis]:
            f_out.write(f"{mis}\t{global_dict[mis]['bowtie']}\t0\n")
        else:
            f_out.write(f"{mis}\t{global_dict[mis]['bowtie']}\t{global_dict[mis]['kraken']}\n")
