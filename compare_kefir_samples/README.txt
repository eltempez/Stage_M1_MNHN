The program compare_kefir_samples.Rmd takes 2 types of inputs :
- global metrics files from the snakemake illumina remapping pipeline, following the flag -i (minimum one file, no maximum)
- the metagenome_ncbi_id.txt file, following the flag -m

The conda environment must first be activated (and installed if needed with the command line "conda env create --name visu_rmd --file=visu_rmd.yml") with "conda activate visu_rmd"

The command line is : R -e "rmarkdown::render('<rmarkdown script>', output_file='<output file name>')" -i <metrics file 1> -i <metrics file 2> -m <metagenome file> 1> /dev/null

Every argument which doesn't follow a -i or -m flag will be ignored. If no output_file is specified, the default will be compare_kefir_samples.html.
/!\ the -i files must be named <sample>_global_metrics.txt, with <sample> being 20 characters at most, so that the graphs don't become cropped