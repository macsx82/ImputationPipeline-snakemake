
# first thing we need is to generate chunks for each chromosome
# this rule will generate a file that will contain the interval string to be used in the imputation. we are using a method similar to the scattergather implementation
# of snakemake, since we want to be able to run multiple chunks togethter in the next rules
rule chunkGenerator:
	wildcard_constraints:
		g_chunk='\d+',
		chr='\d+'
	output:
		output_folder+"/04.impute_intervals/{chr}/{chr}.{g_chunk}.int"
		# directory(output_folder+"/04.impute_intervals/{chr}")
	input:
		ref_hap=config["paths"]["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".hap.gz",
		ref_legend=config["paths"]["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".legend.gz"
	params:
		chunk_size=config['rules']['impute']['chunk_size']
	run:
		# here we will generate the interval string
		# get chr start and end and how many chunks we need for the current chr
		chrom,start,end,chunk_num= get_chunk_num(input.ref_legend,params.chunk_size)
		for chunk in list(range(1,chunk_num+1)):
			out_file=output_folder+"/04.impute_intervals/"+chrom+"/"+chrom+"."+"{:02d}".format(chunk) +".int"
			interval=create_chunks(input.ref_legend,params.chunk_size,chunk)
			open(out_file,"w").write(interval)

# rule to run imputation for each chunk
rule impute:
	wildcard_constraints:
		g_chunk='\d+',
		chr='\d+'
	output:
		expand(output_folder+"/05.imputed/{{chr}}/{{chr}}.{{g_chunk}}.{ext}", ext=["gen.gz","gen_info","gen_info_by_sample","gen_samples","gen_summary","gen_warnings"])
	input:
		ref_hap=config["paths"]["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".hap.gz",
		ref_legend=config["paths"]["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".legend.gz",
		study_geno=rules.phase.output[0],
		study_samples=rules.phase.output[1],
		# interval_file=rules.chunkGenerator.output
	threads:
		config["rules"]["impute"]["threads"]
	resources:
		mem_mb=config["rules"]["impute"]["mem"]
	benchmark:
		output_folder+"/benchmarks/{chr}.{g_chunk}.impute_rule.tsv"
	params:
		impute=config['tools']['impute'],
		ne=config['rules']['impute']['ne'],
		iterations=config['rules']['impute']['iter'],
		burnin=config['rules']['impute']['burnin'],
		k_hap=config['rules']['impute']['k_hap'],
		buffer_size=config['rules']['impute']['buffer_size'],
		interval= lambda wildcards: get_imputation_interval("{output_folder}/04.impute_intervals/{chr}/{chr}.{g_chunk}.int".format(chr=wildcards.chr, g_chunk=wildcards.g_chunk, output_folder=output_folder)),
		# interval=get_imputation_interval('{input.interval_file}'),
		impute_options=config['rules']['impute']['options'],
		gen_map=config['paths']['genetic_map_path']+"/genetic_map_chr{chr}_combined_b37.txt",
		chrx_str=''
	shell:
		"""
		{params.impute} {params.impute_options} -m {params.gen_map} -h {input.ref_hap} -l {input.ref_legend} -known_haps_g {input.study_geno} -sample_g {input.study_samples} $extra_str -iter {params.iterations} -burnin {params.burnin} -k_hap {params.k_hap} -int {params.interval} -Ne {params.ne} -buffer {params.buffer_size} -o {output[0]} {params.chrx_str}
		"""

# this rule is used to gzip the resulting gen files. we decided to not include this step in the previous rule to maximize parallelisation
# but we need to be sure we will have enough disk space
# rule genGzip:
# 	output:
# 	input:
# 	params:
# 	shell:

# let "chunk_num=($chr_end - $chr_begin)/$chunk_size" # bash rounds automatically
# if [[ $chunk_num <1 ]]; then
# 	chunk_num=1
# fi

# #check if command list file exists and remove it
# if [[ -s $imputedir/chr${chr}_command.list ]];then
# 	rm $imputedir/chr${chr}_command.list
# fi

# for chunk in `seq 1 $chunk_num`; do
# 	chunkStr=`printf "%02d" $chunk`
# 	if [[ -e $imputedir/chr$chr.$chunkStr.log ]]; then
# 		continue
# 	fi
# 	if [[ $chunk -gt 1 ]]; then
# 		chunk_begin=`echo "$chr_begin+($chunk-1)*$chunk_size+1" | bc`
# 	else
# 		chunk_begin=$chr_begin
# 	fi
# 	if [[ $chunk -eq $chunk_num ]]; then
# 		mem=${m} #12000
# 		queue=${q}
# 		chunk_end=$chr_end
# 	else
# 		mem=${m} #12000
# 		queue=${q}
# 		chunk_end=`echo "$chr_begin+($chunk*$chunk_size)" | bc`
# 	fi
	
# 	if [[ $by_chunk == "Y" ]]; then
# 		refhap=$refdir/$refname/chr$chr.${chunkStr}$postfix.hap.gz
# 		reflegend=$refdir/$refname/chr$chr.${chunkStr}$postfix.legend.gz
# 	fi
# 	gen_map=${genmap_dir}/genetic_map_chr${chr}_combined_b37.txt

# 	if [[ -s $imputedir/chr$chr.$chunkStr.gen.gz ]];then
# 		echo "Chunk $imputedir/chr$chr.$chunkStr.gen.gz already imputed!!!Skip!"
# 	else
# 	echo -e "#!/usr/bin/env bash
# 	\nmkdir -p ${imputedir}
# 	#$impute2 -allow_large_regions -m ${gen_map} -h $refhap -l $reflegend -known_haps_g $phasedir/chr$chr.haps.gz -sample_g $phasedir/chr$chr.sample $extra_str -use_prephased_g -k_hap $k_hap -int $chunk_begin $chunk_end -Ne 20000 -buffer $buffer_size -o $imputedir/chr$chr.$chunkStr.gen $chrX_impute_str
# 	\n$impute2 -allow_large_regions -m ${gen_map} -h $refhap -l $reflegend -known_haps_g $phasedir/chr$chr.haps.gz -sample_g $phasedir/chr$chr.sample $extra_str -use_prephased_g -iter ${iter} -burnin ${burnin} -k_hap $k_hap -int $chunk_begin $chunk_end -Ne 20000 -buffer $buffer_size -o $imputedir/chr$chr.$chunkStr.gen $chrX_impute_str
# 	\ngzip -f $imputedir/chr$chr.$chunkStr.gen
# 	\nif [[ -e $imputedir/chr$chr.$chunkStr.gen_allele_probs ]]; then
# 	\ngzip $imputedir/chr$chr.$chunkStr.gen_allele_probs $imputedir/chr$chr.$chunkStr.gen_haps
# 	\nfi
# 	\nN_info=\`awk 'NR>1' $imputedir/chr$chr.$chunkStr.gen_info | wc -l | awk '{printf \$1}'\`
# 	\nN_gen=\`zcat $imputedir/chr$chr.$chunkStr.gen.gz | wc -l | awk '{printf \$1}'\`
#             if [[ \$N_info != \$N_gen ]]; then
#                     echo \"chr$chr $chunkStr: \$N_info for info, \$N_gen for gen\" > $imputedir/chr$chr.$chunkStr.ERR
#             fi
#             " > $imputedir/chr$chr.$chunkStr.cmd
#             chmod ug+x $imputedir/chr$chr.$chunkStr.cmd
# 	# cd $imputedir
	
# 	ls $imputedir/chr$chr.$chunkStr.cmd >> $imputedir/chr${chr}_command.list.tmp
# 	fi
# done


