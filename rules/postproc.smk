# rules to generate a report to have a general feel of the imputation quality. By chromosome and by chunk
rule infoStatsChunks:
	wildcard_constraints:
		g_chunk='\d+',
		chr='\d+'
	output:
		output_folder+"/07.stats/{chr}/CHUNKS/{chr}_{g_chunk}_impute_summary_by_maf_by_info.csv",
		output_folder+"/07.stats/{chr}/CHUNKS/{chr}_{g_chunk}_impute_summary_by_maf.csv",
		output_folder+"/07.stats/{chr}/CHUNKS/{chr}_{g_chunk}_impute_summary.pdf",
		output_folder+"/07.stats/{chr}/CHUNKS/{chr}_{g_chunk}_impute_manhattan.png"
	input:
		output_folder+"/06.imputed/{chr}/{chr}.{g_chunk}.vcf.gz"
	params:
		bcftools_bin=config['tools']['bcftools'],
		scripts_folder=config['paths']['scripts'],
		tab_prefix=output_folder+"/07.stats/{chr}/CHUNKS/{chr}_{g_chunk}_impute_summary"
	log:
		stdout=log_folder+"/infoStatsChunks_{chr}_{g_chunk}.o",
		stderr=log_folder+"/infoStatsChunks_{chr}_{g_chunk}.e"
	shell:
		"""
		{params.bcftools_bin} query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%INFO/INFO\n" {input} | {params.scripts_folder}/imputationStats.py --tab {params.tab_prefix} --fig {output[2]} > {log.stdout} 2> {log.stderr}
		{params.bcftools_bin} query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%INFO/INFO\n" {input} | {params.scripts_folder}/plot_manhattan.py --no-log --cols 0,1,6 --title 'Impute INFO score chr{wildcards.chr} chunk {wildcards.g_chunk}' --image {output[3]} --ymax=1.2 - >> {log.stdout} 2>> {log.stderr}
		"""

rule infoStatsChrom:
	wildcard_constraints:
		g_chunk='\d+',
		chr='\d+'
	output:
		output_folder+"/07.stats/{chr}/{chr}_impute_summary_by_maf_by_info.csv",
		output_folder+"/07.stats/{chr}/{chr}_impute_summary_by_maf.csv",
		output_folder+"/07.stats/{chr}/{chr}_impute_summary.pdf",
		output_folder+"/07.stats/{chr}/{chr}_impute_manhattan.png"
	input:
		output_folder+"/06.imputed/MERGED/{chr}/{chr}.vcf.gz",
		
	params:
		bcftools_bin=config['tools']['bcftools'],
		scripts_folder=config['paths']['scripts'],
		tab_prefix=output_folder+"/07.stats/{chr}/{chr}_impute_summary"
	log:
		stdout=log_folder+"/infoStatsChrom_{chr}.o",
		stderr=log_folder+"/infoStatsChrom_{chr}.e"
	shell:
		"""
		{params.bcftools_bin} query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%INFO/INFO\n" {input[0]} | {params.scripts_folder}/imputationStats.py --tab {params.tab_prefix} --fig {output[2]} > {log.stdout} 2> {log.stderr}
		{params.bcftools_bin} query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%INFO/INFO\n" {input[0]} | {params.scripts_folder}/plot_manhattan.py --no-log --cols 0,1,6 --title 'Impute INFO score chr{wildcards.chr}' --image {output[3]} --ymax=1.2 - >> {log.stdout} 2>> {log.stderr}
		"""



# convert VCF to BIMBAM format, to be used with GEMMA and other softwares for analyses
rule convertBimbam:
	output:
		output_folder+"/06.imputed/BIMBAM/{chr}/{chr}.bimbam.gz",
		output_folder+"/06.imputed/BIMBAM/{chr}/{chr}.pos"
	input:
		output_folder+"/06.imputed/MERGED/{chr}/{chr}.vcf.gz"
	params:
		bcftools_bin=config['tools']['bcftools']
	log:
		stdout=log_folder+"/convertBimbam_{chr}.o",
		stderr=log_folder+"/convertBimbam_{chr}.e"
	shell:
		"""
		#We just need to format the vcf file and extract the DS field already present in the VCF, for the main bimbam
		{params.bcftools_bin} query -f'%ID,%REF,%ALT[,%DS]\\n' {input} | gzip --best -c > {output[0]} 2> {log.stderr}
		
		#We just need to format the vcf file and extract the position
		{params.bcftools_bin} query -f'%ID,%POS,%CHROM\\n' {input} -o {output[1]} >> {log.stdout} 2>> {log.stderr}
		"""

