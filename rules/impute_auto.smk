# rule to run imputation for each chunk
rule impute:
	wildcard_constraints:
		g_chunk='\d+',
		chr='\d+'
	output:
		expand(output_folder+"/06.imputed/{{chr}}/{{chr}}.{{g_chunk}}.{ext}", ext=["vcf.gz","log"])
	input:
		interval_file=output_folder+"/05.impute_intervals/{chr}/splitted/{chr}.{g_chunk}.int",
		ref_panel=config["paths"]["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".vcf.gz",
		study_geno=rules.phase.output[0]
		# interval_file=rules.chunkIntervalFileGenerator.output
	threads:
		config["rules"]["impute"]["threads"]
	resources:
		mem_mb=lambda wildcards, attempt: get_mem_mb(config["rules"]["impute"]["mem"], attempt)
	benchmark:
		output_folder+"/benchmarks/{chr}.{g_chunk}.impute_rule.tsv"
	params:
		impute=config['tools']['impute'],
		ne=config['rules']['impute']['ne'],
		pbwt_depth=config['rules']['impute']['pbwt_depth'],
		buffer_size=config['rules']['impute']['buffer_size'],
		interval= lambda wildcards: get_imputation_interval("{output_folder}/05.impute_intervals/{chr}/splitted/{chr}.{g_chunk}.int".format(chr=wildcards.chr, g_chunk=wildcards.g_chunk, output_folder=output_folder)),
		# interval=get_imputation_interval('{input.interval_file}'),
		impute_options=config['rules']['impute']['options'],
		g_map=config['paths']['genetic_map_path']+"/chr{chr}.b37.gmap.gz",
		out_prefix=output_folder+"/06.imputed/{chr}/{chr}.{g_chunk}",
		chrx_str=''
	log:
		stdout=log_folder+"/impute_{chr}_{g_chunk}.o",
		stderr=log_folder+"/impute_{chr}_{g_chunk}.e"
	shell:
		"""
		{params.impute} {params.impute_options} --m {params.g_map} --h {input.ref_panel} --g {input.study_geno} --r {params.interval} --o {params.out_prefix}.vcf.gz --l {params.out_prefix}.log --b {params.buffer_size} --threads {threads} --ne {params.ne} --pbwt-depth {params.pbwt_depth} {params.chrx_str} > {log.stdout} 2> {log.stderr}
		"""


# # rule to concat back data imputed by chromosome
rule concatImputed:
	wildcard_constraints:
		g_chunk='\d+',
		chr='\d+'
	output:
		expand(output_folder+"/06.imputed/MERGED/{{chr}}/{{chr}}.{ext}", ext=["vcf.gz","vcf.gz.tbi","stats"])
	input:
		collect_imputed_chunks
	params:
		bcftools_bin=config['tools']['bcftools'],
		temp=define_tmp(config['rules']['concatImputed']['temp'])
	resources:
		mem_mb=10000
	log:
		stdout=log_folder+"/concatImputed_{chr}.o",
		stderr=log_folder+"/concatImputed_{chr}.e"
	shell:
		"""
		temp=$(mktemp -u -d -p {params.temp})
		{params.bcftools_bin} concat {input}| {params.bcftools_bin} sort -T ${{temp}} -O z -o {output[0]} > {log.stdout} 2> {log.stderr}
		{params.bcftools_bin} index -t {output[0]} >> {log.stdout} 2>> {log.stderr}
		{params.bcftools_bin} stats {output[0]} > {output[2]} 2>> {log.stderr}
		"""
