
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
