
# first thing we need is to generate chunks for each chromosome
# this rule will generate a file that will contain the interval string to be used in the imputation. we are using a method similar to the scattergather implementation
# of snakemake, since we want to be able to run multiple chunks togethter in the next rules
rule chunkGenerator:
	wildcard_constraints:
		g_chunk='\d+',
		chr='\d+'
	output:
		# output_folder+"/05.impute_intervals/{chr}/{chr}.{g_chunk}.int"
		coord_by_chunker=output_folder+"/05.impute_intervals/{chr}/{chr}.coordinates.txt"
		# directory(output_folder+"/04.impute_intervals/{chr}")
	input:
		# regardless of the format of the panel, here we use the reference panel itself
		ref_panel=config["paths"]["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".vcf.gz",
		study_geno=rules.phase.output[0]
	params:
		chunker_tool=config['tools']['chunker_tool'],
		win_size=config['rules']['chunkGenerator']['window-size'],
		win_count=config['rules']['chunkGenerator']['window-count']
		# coord_by_chunker=lambda wildcards: output_folder+"/05.impute_intervals/{chr}/{chr}.coordinates.txt".format(chr=wildcards.chr)
	log:
		stdout=log_folder+"/chunkGenerator_{chr}.o",
		stderr=log_folder+"/chunkGenerator_{chr}.e"
	run:
		chunk_cmd="%s --h %s --r %s --g %s --window-size %s --window-count %s --o %s --l %s 2> %s" %(params.chunker_tool,input.ref_panel, wildcards.chr, input.study_geno,params.win_size,params.win_count,output.coord_by_chunker,log.stdout,log.stderr)
		shell(chunk_cmd)
		# read the generated file and proceed as we did before
		# with open(params.coord_by_chunker) as chunk_file:
		# 	for line in chunk_file:
		# 		chunk=int(line.strip().split("\t")[0])+1
		# 		chrom=line.strip().split("\t")[1]
		# 		interval=line.strip().split("\t")[3]
		# 		out_file=output_folder+"/05.impute_intervals/"+chrom+"/"+chrom+"."+"{:02d}".format(chunk) +".int"
		# 		open(out_file,"w").write(interval)

checkpoint chunkIntervalFileGenerator:
	wildcard_constraints:
		g_chunk='\d+',
		chr='\d+'
	output:
		intervals=directory(output_folder+"/05.impute_intervals/{chr}/splitted")
		# expand(output_folder+"/05.impute_intervals/{{chr}}/{{chr}}.{{g_chunk}}.int")
		# directory(output_folder+"/04.impute_intervals/{chr}")
	input:
		# regardless of the format of the panel, here we use the reference panel itself
		rules.chunkGenerator.output[0]
		# ref_panel=config["paths"]["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".vcf.gz",
		# study_geno=rules.phase.output[0]
	params:
		# coord_by_chunker=lambda wildcards: output_folder+"/05.impute_intervals/{chr}/{chr}.coordinates.txt".format(chr=wildcards.chr)
	log:
		stdout=log_folder+"/chunkIntervalFileGenerator_{chr}.o",
		stderr=log_folder+"/chunkIntervalFileGenerator_{chr}.e"
	run:
		# read the generated file and proceed as we did before
		with open(input[0]) as chunk_file:
			for line in chunk_file:
				chunk=int(line.strip().split("\t")[0])+1
				chrom=line.strip().split("\t")[1]
				interval=line.strip().split("\t")[3]
				out_file=output_folder+"/05.impute_intervals/"+chrom+"/splitted/"+chrom+"."+"{:02d}".format(chunk) +".int"
				createAndOpen(out_file,"w").write(interval)

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
	log:
		stdout=log_folder+"/concatImputed_{chr}.o",
		stderr=log_folder+"/concatImputed_{chr}.e"
	shell:
		"""
		temp=$(mktemp -u -d -p {params.temp})
		{params.bcftools_bin} concat {input}| {params.bcftools_bin} sort -T ${{temp}} -O z -o {output[0]} > {log.stdout} 2> {log.stderr}
		tabix -p vcf {output[0]} >> {log.stdout} 2>> {log.stderr}
		{params.bcftools_bin} stats {output[0]} > {output[2]} 2>> {log.stderr}
		"""
