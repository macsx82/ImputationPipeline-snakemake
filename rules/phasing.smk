#Snakefile for the Internal Imputation pipeline
#
#Author: Massimiliano Cocca

# Phasing rule
rule phase:
    output:
        # generate_shapeit_out_files("{input.chr}")
        touch(output_folder+"/"+config["pop"]+"/{chr}.pipe.done"),
        expand(output_folder+ "03.phased_data/" + ref_panel + "/chr{chr}.{ext}" , ext=["haps.gz","sample"])
        # generate_shapeit_out_files("{chr}")
    input:
        rules.snpFlip.output[0],
        rules.snpFlip.output[1],
        rules.snpFlip.output[2]
    params:
        g_map="/netapp/nfs/resources/1000GP_phase3/impute/genetic_map_chr{chr}_combined_b37.txt",
        input_prefix=output_folder + "/02.flipped_input/" + ref_panel + "/"+ cohort_name+"_{chr}_flipped"
        shapeit=config['tools']['shapeit']
    threads: 16
    benchmark:
        output_folder+"/benchmarks/{chr}.phase_rule.tsv"
    shell:
        # shapeit --input-bed gwas.bed gwas.bim gwas.fam \
        # -M genetic_map.txt \
        # -O gwas.phased
        "{params.shapeit} -B {params.input_prefix} -M {params.g_map} -O {output[1]} {output[2]} -T {threads}"

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