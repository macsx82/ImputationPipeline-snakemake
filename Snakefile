#Snakefile for the Internal Imputation pipeline
#
#Author: Massimiliano Cocca
#
#
# configfile: "config.yaml"
import pathlib

# recover some fixed variables from config file
output_folder = config['paths']["output_folder"]
log_folder = config['paths']["log_folder"]
cohort_name = config["cohort_name"]
input_prefix = config['paths']["input_file_prefix"]
chrs = config['chromosomes']
ref_panel=config['ref_panel']
# define a scatter gather rule to work by chromosome
CHR_COUNT=23

scattergather:
    split=CHR_COUNT

#define parameter useful to cluster job submission
localrules: all

#define rules
rule all:
    # wildcard_constraints:
    #     g_chunk='\d+'
    input:
        # config["output_folder"]+"/"+config["pop"]+"/" + config["pop"] + ".pipe.done"
        # lambda wildcards: config["chr"][wildcards.chrom],
        # expand(config["output_folder"]+"/"+config["pop"]+"/{chrom}.pipe.done", chrom=config["chr"])
        # config["output_folder"]+"/"+config["pop"]+"/" + config["chr"] + ".pipe.done"
        # expand(output_folder+"/00.splitted_input/"+cohort_name+"_{chr}.{ext}", ext=['bed','bim','fam'],chr=chrs)
        # expand(output_folder+"/01.refAlign/"+ref_panel+"/{chr}_shapeit_"+ref_panel+".alignments.snp.{ext}", ext=['strand','strand.exclude'],chr=chrs)
        # expand(output_folder+"/01.refAlign/"+ref_panel+"/{chr}_shapeit_rsids.to_flip",chr=chrs)
        # expand(output_folder + "/02.flipped_input/" + ref_panel + "/"+ cohort_name+"_{chr}_flipped.{ext}",ext=['bed','bim','fam'],chr=chrs)
        # expand(output_folder+ "/03.phased_data/" + ref_panel + "/chr{chr}.{ext}" , ext=["haps.gz","sample"], chr=chrs)
        # expand(output_folder+"/04.impute_intervals/{chr}/{chr}.{g_chunk}.int",chr=chrs,g_chunk=list(range(1,11)))
        # expand(output_folder+"/04.impute_intervals/{chr}/{chr}.{{g_chunk}}.pippo",chr=chrs)
        output_folder+"/04.impute_intervals/{{chr}}/{{chr}}.{{g_chunk}}.pippo"
        # directory(expand(output_folder+"/04.impute_intervals/{chr}/",chr=chrs))
        # expand(output_folder+"/04.impute_intervals/{chr}/{chr}.{{g_chunk}}.pippo",chr=chrs)

# MODULES
include_prefix="rules"
include:
    include_prefix + "/functions.py"
include:
    include_prefix + "/preproc.smk"
include:
    include_prefix + "/phasing.smk"
include:
    include_prefix + "/impute.smk"


onsuccess:
    print("The workflow finished without errors!")

onerror:
    print("An error occurred in the current workflow execution!!")