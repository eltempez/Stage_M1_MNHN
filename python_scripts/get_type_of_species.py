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
# list of microorganisms
dict_micro = {"aceto": ["B3", "B4", "B6", "B1", "B5", "B2"], "lacto": ["B11", "B18", "B10", "B7", "B14", "B13", "B17", "B21", "B12", "B15"], "yeasts": ["BB", "Y1", "Y2", "Y3"]}
# directory of metrics files
dir = "/home/eliotttempez/Documents/tests/outputs_snk1/mapping/mapped/metrics_kA/"
# number of total reads for kA
nb_reads = 8612900*2
# number of total reads for MCF
#nb_reads = 4611205*2


## Main
nb_aceto, nb_lacto, nb_yeasts = 0, 0, 0
for file_name in os.listdir(dir):
    species_name = file_name.split("_")[0]
    file_path = os.path.join(dir, file_name)
    with open(file_path, "r") as metrics_file:
        ligne = metrics_file.readlines()[2].strip()
        nb_mapped = int(ligne.split("\t")[-1])
    if species_name in dict_micro["aceto"]:
        nb_aceto += nb_mapped
    elif species_name in dict_micro["lacto"]:
        nb_lacto += nb_mapped
    elif species_name in dict_micro["yeasts"]:
        nb_yeasts += nb_mapped

cov_aceto = round((nb_aceto/nb_reads)*100)
cov_lacto = round((nb_lacto/nb_reads)*100)
cov_yeasts = round((nb_yeasts/nb_reads)*100)
print("lacto aceto yeasts")
print(f"{cov_lacto} {cov_aceto} {cov_yeasts}")


