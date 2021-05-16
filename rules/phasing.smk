#Snakefile for the Internal Imputation pipeline
#
#Author: Massimiliano Cocca

# Phasing rule
rule phase:
    input:
        # expand(config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_flipped" , ext=[".bim",".bed",".fam"]),
        rules.snp_flip.output[0],
        rules.snp_flip.output[1],
        rules.snp_flip.output[2]
    params:
        g_map="/netapp/nfs/resources/1000GP_phase3/impute/genetic_map_chr"+config["chr"]+"_combined_b37.txt",
        input_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_" + config["chr"] + "_flipped"
    output:
        # generate_shapeit_out_files("{input.chr}")
        touch(config["output_folder"]+"/"+config["pop"]+"/"+config["chr"] +".pipe.done"),
        generate_shapeit_out_files(config["chr"])
    threads: 16
    benchmark:
        config["output_folder"]+"/"+config["pop"]+"/"+config["chr"] +".phase_rule.tsv"
    shell:
        # shapeit --input-bed gwas.bed gwas.bim gwas.fam \
        # -M genetic_map.txt \
        # -O gwas.phased
        "{config[shapeit_path]} -B {params.input_prefix} -M {params.g_map} -O {output[1]} {output[2]} -T {threads}"

# rule pipe_finish:
#     input:
#         expand(config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_flipped" , ext=[".bim",".bed",".fam"]),
#         rules.snp_flip.output.strand_rsid
#         # expand(config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_flipped" , ext=[".bim",".bed",".fam"])
#         # expand(config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_shapeit_refpanel.alignments" , ext=[".strand",".strand.exclude"])
#     output:
#         config["output_folder"]+"/"+config["pop"]+"/"+config["chr"] +".pipe.done"
#     shell:
#         "touch {output}"

onsuccess:
    print("The workflow finished without errors!")

onerror:
    print("An error occurred in the current workflow execution!!")