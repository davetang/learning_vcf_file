## README

The `run.sh` script first creates the example files by running the parameterised R Markdown file `variant.Rmd`. Next, the script uses the conda environment created from `../environment.yml` (create the environment before running the script) to map the example reads and call variants using `BCFtools`, `FreeBayes`, and `GATK`. You will need to download the GATK zip file and unzip the files into this directory (`analysis`). If everything works accordingly, the final VCF files are generated and moved to the directory `../result`.

