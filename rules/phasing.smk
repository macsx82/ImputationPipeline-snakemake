#Snakefile for the Internal Imputation pipeline
#
#Author: Massimiliano Cocca

# Phasing rule tailored to use SHAPEIT2
# rule phase:
#     output:
#         # generate_shapeit_out_files("{input.chr}")
#         output_folder+ "/03.phased_data/" + ref_panel + "/chr{chr}.haps.gz",
#         output_folder+ "/03.phased_data/" + ref_panel + "/chr{chr}.sample"
#         # generate_shapeit_out_files("{chr}")
#     input:
#         rules.snpFlip.output[0],
#         rules.snpFlip.output[1],
#         rules.snpFlip.output[2]
#     params:
#         g_map=config['paths']['genetic_map_path']+"/genetic_map_chr{chr}_combined_b37.txt",
#         input_prefix=output_folder + "/02.flipped_input/" + ref_panel + "/"+ cohort_name+"_{chr}_flipped",
#         shapeit=config['tools']['shapeit']
#     threads: 16
#     benchmark:
#         output_folder+"/benchmarks/{chr}.phase_rule.tsv"
#     shell:
#         "{params.shapeit} -B {params.input_prefix} -M {params.g_map} -O {output[0]} {output[1]} -T {threads}"

# Phasing rule tailored to use SHAPEIT4
rule phase:
    output:
        output_folder+ "/04.phased_data/" + ref_panel + "/"+ cohort_name +"_{chr}_phased.vcf.gz",
        output_folder+ "/04.phased_data/" + ref_panel + "/"+ cohort_name +"_{chr}_phased.vcf.gz.tbi"
    input:
        # rules.plink2vcf.output[0]
        rules.vcfAnnotate.output[0]
    params:
        g_map=config['paths']['genetic_map_path']+"/chr{chr}.b37.gmap.gz",
        mcmc_iterations=config['rules']['phase']['mcmc_iterations'],
        pbwt_depth=config['rules']['phase']['pbwt_depth'],
        phasing_tool=config['tools']['phasing_tool']
    threads: 16
    log:
        output_folder+ "/04.phased_data/" + ref_panel + "/"+ cohort_name +"_{chr}_phase.log"
    benchmark:
        output_folder+"/benchmarks/{chr}.phase_rule.tsv"
    shell:
        """
        {params.phasing_tool} --input {input[0]} --map {params.g_map} --region {wildcards.chr} --output {output[0]} --thread {threads} --log {log[0]} --mcmc-iterations {params.mcmc_iterations} --pbwt-depth {params.pbwt_depth}
        tabix -p vcf {output[0]}
        """

