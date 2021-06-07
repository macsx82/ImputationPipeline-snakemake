
# first thing we need is to generate chunks for each chromosome
# this rule will generate a file that will contain the interval string to be used in the imputation. we are using a method similar to the scattergather implementation
# of snakemake, since we want to be able to run multiple chunks togethter in the next rules
rule chunkGenerator:
	wildcard_constraints:
		g_chunk='\d+',
		chr='\d+'
	output:
		output_folder+"/05.impute_intervals/{chr}/{chr}.{g_chunk}.int"
		# directory(output_folder+"/04.impute_intervals/{chr}")
	input:
		# regardless of the format of the panel, here we use the reference panel itself
		ref_panel=config["paths"]["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".vcf.gz",
		study_geno=rules.phase.output[0]
	params:
		chunker_tool=config['tools']['chunker_tool'],
		coord_by_chunker=output_folder+"/05.impute_intervals/{chr}/{chr}.coordinates.txt".format(chr=wildcards.chr)
	run:
		chunk_cmd="%s --h %s --r %s --g %s --o %s" %(params.chunker_tool,input.ref_panel, wildcards.chr, input.study_geno, params.coord_by_chunker)
		shell(chunk_cmd)
		# read the generated file and proceed as we did before
		with open(params.coord_by_chunker) as chunk_file:
			for line in chunk_file:
				chunk=int(line.strip().split("\t")[0])+1
				chrom=line.strip().split("\t")[1]
				interval=line.strip().split("\t")[3]
				out_file=output_folder+"/05.impute_intervals/"+chrom+"/"+chrom+"."+"{:02d}".format(chunk) +".int"
				open(out_file,"w").write(interval)
	# run:
	# 	# here we will generate the interval string
	# 	# get chr start and end and how many chunks we need for the current chr
	# 	chrom,start,end,chunk_num= get_chunk_num(input.ref_legend,params.chunk_size)
	# 	for chunk in list(range(1,chunk_num+1)):
	# 		out_file=output_folder+"/05.impute_intervals/"+chrom+"/"+chrom+"."+"{:02d}".format(chunk) +".int"
	# 		interval=create_chunks(input.ref_legend,params.chunk_size,chunk)
	# 		open(out_file,"w").write(interval)

# rule to run imputation for each chunk
# rule impute:
# 	wildcard_constraints:
# 		g_chunk='\d+',
# 		chr='\d+'
# 	output:
# 		expand(output_folder+"/05.imputed/{{chr}}/{{chr}}.{{g_chunk}}.{ext}", ext=["gen.gz","gen_info","gen_info_by_sample","gen_samples","gen_summary","gen_warnings"])
# 	input:
# 		ref_hap=config["paths"]["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".hap.gz",
# 		ref_legend=config["paths"]["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".legend.gz",
# 		study_geno=rules.phase.output[0],
# 		study_samples=rules.phase.output[1],
# 		# interval_file=rules.chunkGenerator.output
# 	threads:
# 		config["rules"]["impute"]["threads"]
# 	resources:
# 		mem_mb=config["rules"]["impute"]["mem"]
# 	benchmark:
# 		output_folder+"/benchmarks/{chr}.{g_chunk}.impute_rule.tsv"
# 	params:
# 		impute=config['tools']['impute'],
# 		ne=config['rules']['impute']['ne'],
# 		iterations=config['rules']['impute']['iter'],
# 		burnin=config['rules']['impute']['burnin'],
# 		k_hap=config['rules']['impute']['k_hap'],
# 		buffer_size=config['rules']['impute']['buffer_size'],
# 		interval= lambda wildcards: get_imputation_interval("{output_folder}/04.impute_intervals/{chr}/{chr}.{g_chunk}.int".format(chr=wildcards.chr, g_chunk=wildcards.g_chunk, output_folder=output_folder)),
# 		# interval=get_imputation_interval('{input.interval_file}'),
# 		impute_options=config['rules']['impute']['options'],
# 		gen_map=config['paths']['genetic_map_path']+"/genetic_map_chr{chr}_combined_b37.txt",
# 		out_prefix=output_folder+"/05.imputed/{chr}/{chr}.{g_chunk}.gen",
# 		chrx_str=''
# 	shell:
# 		"""
# 		{params.impute} {params.impute_options} -m {params.gen_map} -h {input.ref_hap} -l {input.ref_legend} -known_haps_g {input.study_geno} -sample_g {input.study_samples} -iter {params.iterations} -burnin {params.burnin} -k_hap {params.k_hap} -int {params.interval} -Ne {params.ne} -buffer {params.buffer_size} -o {params.out_prefix} {params.chrx_str}
# 		"""
# Rules preserved here and custom made for impute2 and impute2 reference panels
# 
# rule chunkGenerator:
# 	wildcard_constraints:
# 		g_chunk='\d+',
# 		chr='\d+'
# 	output:
# 		output_folder+"/04.impute_intervals/{chr}/{chr}.{g_chunk}.int"
# 		# directory(output_folder+"/04.impute_intervals/{chr}")
# 	input:
# 		ref_hap=config["paths"]["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".hap.gz",
# 		ref_legend=config["paths"]["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".legend.gz"
# 	params:
# 		chunk_size=config['rules']['impute']['chunk_size']
# 	run:
# 		# here we will generate the interval string
# 		# get chr start and end and how many chunks we need for the current chr
# 		chrom,start,end,chunk_num= get_chunk_num(input.ref_legend,params.chunk_size)
# 		for chunk in list(range(1,chunk_num+1)):
# 			out_file=output_folder+"/04.impute_intervals/"+chrom+"/"+chrom+"."+"{:02d}".format(chunk) +".int"
# 			interval=create_chunks(input.ref_legend,params.chunk_size,chunk)
# 			open(out_file,"w").write(interval)

# # rule to run imputation for each chunk
# rule impute:
# 	wildcard_constraints:
# 		g_chunk='\d+',
# 		chr='\d+'
# 	output:
# 		expand(output_folder+"/05.imputed/{{chr}}/{{chr}}.{{g_chunk}}.{ext}", ext=["gen.gz","gen_info","gen_info_by_sample","gen_samples","gen_summary","gen_warnings"])
# 	input:
# 		ref_hap=config["paths"]["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".hap.gz",
# 		ref_legend=config["paths"]["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".legend.gz",
# 		study_geno=rules.phase.output[0],
# 		study_samples=rules.phase.output[1],
# 		# interval_file=rules.chunkGenerator.output
# 	threads:
# 		config["rules"]["impute"]["threads"]
# 	resources:
# 		mem_mb=config["rules"]["impute"]["mem"]
# 	benchmark:
# 		output_folder+"/benchmarks/{chr}.{g_chunk}.impute_rule.tsv"
# 	params:
# 		impute=config['tools']['impute'],
# 		ne=config['rules']['impute']['ne'],
# 		iterations=config['rules']['impute']['iter'],
# 		burnin=config['rules']['impute']['burnin'],
# 		k_hap=config['rules']['impute']['k_hap'],
# 		buffer_size=config['rules']['impute']['buffer_size'],
# 		interval= lambda wildcards: get_imputation_interval("{output_folder}/04.impute_intervals/{chr}/{chr}.{g_chunk}.int".format(chr=wildcards.chr, g_chunk=wildcards.g_chunk, output_folder=output_folder)),
# 		# interval=get_imputation_interval('{input.interval_file}'),
# 		impute_options=config['rules']['impute']['options'],
# 		gen_map=config['paths']['genetic_map_path']+"/genetic_map_chr{chr}_combined_b37.txt",
# 		out_prefix=output_folder+"/05.imputed/{chr}/{chr}.{g_chunk}.gen",
# 		chrx_str=''
# 	shell:
# 		"""
# 		{params.impute} {params.impute_options} -m {params.gen_map} -h {input.ref_hap} -l {input.ref_legend} -known_haps_g {input.study_geno} -sample_g {input.study_samples} -iter {params.iterations} -burnin {params.burnin} -k_hap {params.k_hap} -int {params.interval} -Ne {params.ne} -buffer {params.buffer_size} -o {params.out_prefix} {params.chrx_str}
# 		"""