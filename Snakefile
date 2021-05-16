#Snakefile for the Internal Imputation pipeline
#
#Author: Massimiliano Cocca
#
#
# configfile: "config.yaml"
   

input_prefix=config['input_file_prefix']
output_folder=config['output_folder']
cohort_name=config['cohort_name']

# define a scatter gather rule to work by chromosome
CHR_COUNT=23

scattergather:
    split=CHR_COUNT

#define parameter useful to cluster job submission
localrules: all

#define rules
rule all:
    input:
        # config["output_folder"]+"/"+config["pop"]+"/" + config["pop"] + ".pipe.done"
        # lambda wildcards: config["chr"][wildcards.chrom],
        # expand(config["output_folder"]+"/"+config["pop"]+"/{chrom}.pipe.done", chrom=config["chr"])
        # config["output_folder"]+"/"+config["pop"]+"/" + config["chr"] + ".pipe.done"
        scatter.split("{{output_folder}}/00.splitted_input/{scatteritem}_{{cohort_name}}.bed"),
        scatter.split("{{output_folder}}/00.splitted_input/{scatteritem}_{{cohort_name}}.bim"),
        scatter.split("{{output_folder}}/00.splitted_input/{scatteritem}_{{cohort_name}}.fam")

# MODULES
include_prefix="rules"
include:
    include_prefix + "/functions.py"
include:
    include_prefix + "/preproc.smk"
# include:
#     include_prefix + "/phasing.smk"
# include:
#     include_prefix + "/impute.smk"


onsuccess:
    print("The workflow finished without errors!")

onerror:
    print("An error occurred in the current workflow execution!!")