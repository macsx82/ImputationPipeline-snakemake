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
        output_folder + "/03.flipped_input/" + ref_panel + "/VCF/"+ cohort_name+"_{chr}_fixRef_sorted_rsID.vcf.gz"
        # rules.vcfAnnotate.output[0]
    params:
        g_map=config['paths']['genetic_map_path']+"/chr{chr}.b37.gmap.gz",
        mcmc_iterations=config['rules']['phase']['mcmc_iterations'],
        pbwt_depth=config['rules']['phase']['pbwt_depth'],
        additional_args=config['rules']['phase']['additional_args'],
        phasing_tool=config['tools']['phasing_tool'],
        region_chr=getChrForPhasing
    threads: 16
    log:
        output_folder+ "/04.phased_data/" + ref_panel + "/"+ cohort_name +"_{chr}_phase.log",
        stdout=log_folder+"/phase_{chr}.o",
        stderr=log_folder+"/phase_{chr}.e"
    benchmark:
        output_folder+"/benchmarks/{chr}.phase_rule.tsv"
    shell:
        """
        {params.phasing_tool} --input {input[0]} --map {params.g_map} --region {params.region_chr} --output {output[0]} --thread {threads} --log {log[0]} {params.mcmc_iterations} {params.pbwt_depth} {params.additional_args} 2> {log.stderr}
        tabix -p vcf {output[0]}
        """

