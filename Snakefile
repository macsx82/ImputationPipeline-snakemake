#Snakefile for the Integrated Prioritization Workflow pipeline
#
#Author: Massimiliano Cocca
#
#
# configfile: "config.yaml"

print (config['chr_to_phase'].keys())


def generate_shapeit_out_files(key):
    chr_phased= "%s/%s/%s/chr%s.haps.gz" % (config["output_folder"],config["pop"],key,key)
    samples= "%s/%s/%s/chr%s.samples" % (config["output_folder"],config["pop"],key,key)

    return chr_phased,samples

def generate_end_of_pipeline_files(key):
    return "%s/%s/chr%s.pipe.done" % (config["output_folder"],config["pop"],key)

#define parameter useful to cluster job submission
localrules: all

#define rules
rule all:
    input:
        config["output_folder"]+"/"+config["pop"]+"/"+config["chr"] +".pipe.done"

#We assume all our data has been already strand oriented with plink
#We will orient the data to match the reference panel orientation, using shapeit
# Input files will be the plink genotypes. We will get them from a config file
rule snp_flip:
    input:
        ug_bed=config["input_folder"] + "/" + config["chr"]+ ".bed",
        ug_bim=config["input_folder"] + "/" + config["chr"]+ ".bim",
        ug_fam=config["input_folder"] + "/" + config["chr"]+ ".fam",
        rp_hap=config["ref_panel_base_folder"]+ "/" + config["ref_panel"]+ "/" + config["chr"]+ "/" + config["chr"] + "." + config["ref_panel"] + ".hap.gz",
        rp_legend=config["ref_panel_base_folder"]+ "/" + config["ref_panel"]+ "/" + config["chr"]+ "/" + config["chr"] + "." + config["ref_panel"] + ".legend.gz",
        rp_samples=config["ref_panel_base_folder"]+ "/" + config["ref_panel"]+ "/" + config["chr"]+ "/" + config["chr"] + "." + config["ref_panel"] + ".samples"
    params:
        g_map="/netapp/nfs/resources/1000GP_phase3/impute/genetic_map_chr"+config["chr"]+"_combined_b37.txt",
        output_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_shapeit_refpanel.alignments"
    output:
        # expand(config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]+ "/chr"+config["chr"]+"_relate_pos_sel{ext}", ext=[".freq",".lin",".sele"])
        expand(config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_shapeit_refpanel.alignments" , ext=[".strand",".strand.exclude"])
    shell:
        "{config[shapeit_path]} --check --input-bed {input.ug_bed} {input.ug_bim} {input.ug_fam} -M {params.g_map} --input-ref {input.rp_hap} {input.rp_legend} {input.rp_samples} --output-log {params.output_prefix}"

# rule phase:
#     input:
#         # g_map=config["genetic_map_chr"],
#         # lambda wildcards: config["chr_to_phase"][wildcards.chr]
#         # chr_to_phase="{input_folder}/{chr}.vcf.gz"
#         # "{params.input_f}/{chr}.vcf.gz"
#         # config["input_folder"] + "/{chr}.vcf.gz"
#         config["input_folder"] + "/" + config["chr"]+ ".vcf.gz"
#     params:
#         input_f=config["input_folder"],
#         g_map="/netapp/nfs/resources/1000GP_phase3/impute/genetic_map_chr"+config["chr"]+"_combined_b37.txt",
#         # base_out=config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]
#     output:
#         # generate_shapeit_out_files("{input.chr}")
#         generate_shapeit_out_files(config["chr"])
#         # chr_phased=config["output_folder"] + "/" + config["pop"] + "/" + config["chr"] + "/chr"+config["chr"]+".haps.gz",
#         # samples= params[base_out] + "/chr"+config["chr"]+".samples"
#     threads: 4
#     shell:
#         # "shapeit -V {input_f}/{input} -M {g_map} -O {output.chr_phased} {output.samples} -T {threads}"
#         # "{config[shapeit_path]} -V {input} -M {params.g_map} -O {output.chr_phased} {output.samples} -T {threads}"
#         # shapeit --input-bed gwas.bed gwas.bim gwas.fam \
#         # -M genetic_map.txt \
#         # -O gwas.phased
#         "{config[shapeit_path]} -V {input} -M {params.g_map} -O {output[0]} {output[1]} -T {threads}"
#         # "touch {output.chr_phased} {output.samples}"

rule pipe_finish:
    input:
        expand(config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_shapeit_refpanel.alignments" , ext=[".strand",".strand.exclude"])
    output:
        config["output_folder"]+"/"+config["pop"]+"/"+config["chr"] +".pipe.done"
    shell:
        "touch {output}"

onsuccess:
    print("The workflow finished without errors!")

onerror:
    print("An error occurred in the current workflow execution!!")