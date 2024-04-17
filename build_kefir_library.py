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
dict_ncbi = {"B3": 483199, "B4": 318683, "B6": 38307, "B1": 442, "B5": 28448, "B2": 436, "B11": 1597, "B18": 1588, "B10": 152331, "B7": 1245, "B14": 399370, "B13": 82688, "B17": 259059, "B21": 2203724, "B12": 304207, "BB": 5007, "Y3": 4926, "Y1": 4932, "Y2": 48255}
# directory of .fa/.fasta files
dir = "/home/eliotttempez/Documents/tests/outputs_snk1/genome_data/fa/"


## Main
# for each file in the matagenome
for file in os.listdir(dir):
    species_name = get_prefix(file)
    with open(dir + file) as fa_file:
        for line in fa_file:
            if line.startswith(">"):
                line += f"|kraken:taxid|{dict_ncbi[species_name]}"


