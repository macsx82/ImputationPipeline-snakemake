
# first thing we need is to generate chunks for each chromosome
# this rule will generate a file that will contain the interval string to be used in the imputation. we are using a method similar to the scattergather implementation
# of snakemake, since we want to be able to run multiple chunks togethter in the next rules
rule chunkGenerator:
	wildcard_constraints:
    	chunk='\d+'
	output:
		output_folder+"/04.impute_intervals/{chr}.{chunk}.int"
	input:
		ref_hap=config["paths"]["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".hap.gz",
		ref_legend=config["paths"]["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".legend.gz",
		chunk_size=config['rules']['impute']['chunk_size']
	run:
		# here we will generate the interval string
		# for chr in chrs:
		# get chr start and end and how many chunks we need for the current chr
		start,end,chunk_num= get_chunk_num(params.ref_legend,params.chunk_size)
		for chunk in list(range(1,chunk_num+1)):
			out_file="%s/04.impute_intervals/{chr}.%s.int" % (output_folder,"{:02d}".format(chunk))
			open(out_file,"w").write(create_chunks(params.ref_legend,params.chunk_size,chunk))
			open(out_file, 'a').close()
			# create_chunks(params.ref_legend,params.chunk_size,chunk) > output_folder+"/04.impute_intervals/{chr}..int"

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


