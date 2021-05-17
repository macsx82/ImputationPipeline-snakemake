#Snakefile for the Internal Imputation pipeline
#
#Author: Massimiliano Cocca

# Phasing rule
rule phase:
    output:
        # generate_shapeit_out_files("{input.chr}")
        output_folder+ "/03.phased_data/" + ref_panel + "/chr{chr}.haps.gz",
        output_folder+ "/03.phased_data/" + ref_panel + "/chr{chr}.sample"
        # generate_shapeit_out_files("{chr}")
    input:
        rules.snpFlip.output[0],
        rules.snpFlip.output[1],
        rules.snpFlip.output[2]
    params:
        g_map=config['paths']['genetic_map_path']+"/genetic_map_chr{chr}_combined_b37.txt",
        input_prefix=output_folder + "/02.flipped_input/" + ref_panel + "/"+ cohort_name+"_{chr}_flipped",
        shapeit=config['tools']['shapeit']
    threads: 16
    benchmark:
        output_folder+"/benchmarks/{chr}.phase_rule.tsv"
    shell:
        "{params.shapeit} -B {params.input_prefix} -M {params.g_map} -O {output[0]} {output[1]} -T {threads}"

