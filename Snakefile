#Snakefile for the Internal Imputation pipeline
#
#Author: Massimiliano Cocca
#
#
# configfile: "config.yaml"

# recover some fixed variables from config file
output_folder = config['paths']["output_folder"]
log_folder = config['paths']["log_folder"]
cohort_name = config["cohort_name"]
input_prefix = config['paths']["input_file_prefix"]
chrs = config['chromosomes']
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
        expand(output_folder+"/00.splitted_input/"+cohort_name+"_{chr}.{ext}", ext=['bed','bim','fam'],chr=chrs)
        # scatter.split(output_folder+"/00.splitted_input/{scatteritem}_"+cohort_name+".bed"),
        # scatter.split(output_folder+"/00.splitted_input/{scatteritem}_"+cohort_name+".bim"),
        # scatter.split(output_folder+"/00.splitted_input/{scatteritem}_"+cohort_name+".fam")

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