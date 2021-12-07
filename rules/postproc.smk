# rules to generate a report to have a general feel of the imputation quality. By chromosome and by chunk
rule infoStatsChunks:
	wildcard_constraints:
		g_chunk='\d+',
		chr='\d+'
	output:
		output_folder+"/07.stats/{chr}/CHUNKS/{chr}_{g_chunk}_impute_summary_by_maf_by_info.csv",
		output_folder+"/07.stats/{chr}/CHUNKS/{chr}_{g_chunk}_impute_summary_by_maf.csv",
		output_folder+"/07.stats/{chr}/CHUNKS/{chr}_{g_chunk}_impute_summary.png",
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
	priority: 50
	shell:
		"""
		{params.bcftools_bin} query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%INFO/INFO\n" {input} | {params.scripts_folder}/imputationStats.py --tab {params.tab_prefix} --fig {output[2]} > {log.stdout} 2> {log.stderr}
		{params.bcftools_bin} query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%INFO/INFO\n" {input} | {params.scripts_folder}/plot_manhattan.py --chunk --no-log --cols 0,1,6 --title 'Impute INFO score chr{wildcards.chr} chunk {wildcards.g_chunk}' --image {output[3]} --ymax=1.2 - >> {log.stdout} 2>> {log.stderr}
		"""

rule infoStatsChrom:
	wildcard_constraints:
		g_chunk='\d+',
		chr='\d+'
	output:
		output_folder+"/07.stats/{chr}/{chr}_impute_summary_by_maf_by_info.csv",
		output_folder+"/07.stats/{chr}/{chr}_impute_summary_by_maf.csv",
		output_folder+"/07.stats/{chr}/{chr}_impute_summary.png",
		output_folder+"/07.stats/{chr}/{chr}_impute_manhattan.png"
	input:
		output_folder+"/06.imputed/MERGED/{chr}/{chr}.vcf.gz"
	params:
		bcftools_bin=config['tools']['bcftools'],
		scripts_folder=config['paths']['scripts'],
		tab_prefix=output_folder+"/07.stats/{chr}/{chr}_impute_summary"
	log:
		stdout=log_folder+"/infoStatsChrom_{chr}.o",
		stderr=log_folder+"/infoStatsChrom_{chr}.e"
	priority: 49
	shell:
		"""
		{params.bcftools_bin} query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%INFO/INFO\n" {input[0]} | {params.scripts_folder}/imputationStats.py --tab {params.tab_prefix} --fig {output[2]} > {log.stdout} 2> {log.stderr}
		{params.bcftools_bin} query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%INFO/INFO\n" {input[0]} | {params.scripts_folder}/plot_manhattan.py --no-log --cols 0,1,6 --title 'Impute INFO score chr{wildcards.chr}' --image {output[3]} --ymax=1.2 - >> {log.stdout} 2>> {log.stderr}
		"""


# #aggregator rule to get all data from all chunks and generate a single pdf file
# rule pdfReportCunks:
# 	wildcard_constraints:
# 		g_chunk='\d+',
# 		chr='\d+'
# 	output:
# 		output_folder+"/07.stats/{chr}/{chr}_impute_summary_report_by_chunk.pdf"
# 	input:
# 		# chunk_stats_by_maf_by_info=expand(output_folder+"/07.stats/{{chr}}/CHUNKS/{{chr}}_{g_chunk}_impute_summary_by_maf_by_info.csv",zip,**glob_wildcards(os.path.join(output_folder+"/06.imputed/{chr}/", "{chr}.{g_chunk}.vcf.gz"))._asdict()),
# 		# chunk_stats_by_maf=expand(output_folder+"/07.stats/{{chr}}/CHUNKS/{{chr}}_{g_chunk}_impute_summary_by_maf.csv",zip,**glob_wildcards(os.path.join(output_folder+"/06.imputed/{chr}/", "{chr}.{g_chunk}.vcf.gz"))._asdict()),
# 		# info_af=expand(output_folder+"/07.stats/{{chr}}/CHUNKS/{{chr}}_{g_chunk}_impute_summary.png",zip,**glob_wildcards(os.path.join(output_folder+"/06.imputed/{chr}/", "{chr}.{g_chunk}.vcf.gz"))._asdict()),
# 		# manhattan=expand(output_folder+"/07.stats/{{chr}}/CHUNKS/{{chr}}_{g_chunk}_impute_manhattan.png",zip,**glob_wildcards(os.path.join(output_folder+"/06.imputed/{chr}/", "{chr}.{g_chunk}.vcf.gz"))._asdict())
# 		# chunk_stats_by_maf_by_info=lambda wildcards: expand(output_folder+"/07.stats/{chr}/CHUNKS/{chr}_{g_chunk}_impute_summary_by_maf_by_info.csv",chr=wildcards.chr,g_chunk=glob_wildcards(os.path.join(output_folder+"/06.imputed/"+wildcards.chr+"/", wildcards.chr+".{g_chunk}.vcf.gz")).g_chunk),
# 		# chunk_stats_by_maf=lambda wildcards: expand(output_folder+"/07.stats/{chr}/CHUNKS/{chr}_{g_chunk}_impute_summary_by_maf.csv",chr=wildcards.chr,g_chunk=glob_wildcards(os.path.join(output_folder+"/06.imputed/"+wildcards.chr+"/", wildcards.chr+".{g_chunk}.vcf.gz")).g_chunk),
# 		# info_af=lambda wildcards: expand(output_folder+"/07.stats/{chr}/CHUNKS/{chr}_{g_chunk}_impute_summary.png",chr=wildcards.chr,g_chunk=glob_wildcards(os.path.join(output_folder+"/06.imputed/"+wildcards.chr+"/", wildcards.chr+".{g_chunk}.vcf.gz")).g_chunk),
# 		# manhattan=lambda wildcards: expand(output_folder+"/07.stats/{chr}/CHUNKS/{chr}_{g_chunk}_impute_manhattan.png",chr=wildcards.chr,g_chunk=glob_wildcards(os.path.join(output_folder+"/06.imputed/"+wildcards.chr+"/", wildcards.chr+".{g_chunk}.vcf.gz")).g_chunk)
# 		chunk_stats_by_maf_by_info=lambda wildcards: expand(rules.infoStatsChunks.output[0],chr=wildcards.chr,g_chunk=glob_wildcards(os.path.join(output_folder+"/06.imputed/"+wildcards.chr+"/", wildcards.chr+".{g_chunk}.vcf.gz")).g_chunk),
# 		chunk_stats_by_maf=lambda wildcards: expand(rules.infoStatsChunks.output[1],chr=wildcards.chr,g_chunk=glob_wildcards(os.path.join(output_folder+"/06.imputed/"+wildcards.chr+"/", wildcards.chr+".{g_chunk}.vcf.gz")).g_chunk),
# 		info_af=lambda wildcards: expand(rules.infoStatsChunks.output[2],chr=wildcards.chr,g_chunk=glob_wildcards(os.path.join(output_folder+"/06.imputed/"+wildcards.chr+"/", wildcards.chr+".{g_chunk}.vcf.gz")).g_chunk),
# 		manhattan=lambda wildcards: expand(rules.infoStatsChunks.output[3],chr=wildcards.chr,g_chunk=glob_wildcards(os.path.join(output_folder+"/06.imputed/"+wildcards.chr+"/", wildcards.chr+".{g_chunk}.vcf.gz")).g_chunk)
	
# 	params:
# 		stat_base_folder=output_folder+"/07.stats/{chr}/CHUNKS"
# 	log:
# 		stdout=log_folder+"/pdfReportCunks_{chr}.o",
# 		stderr=log_folder+"/pdfReportCunks_{chr}.e"	
# 	priority: 1
# 	run:
# 		chunk_number=len(input.chunk_stats_by_maf_by_info)
# 		pdf_report(wildcards.chr,params.stat_base_folder,chunk_number,output[0])



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

#Add the correct annotation name for the impute info score. Impute5 creates a INFO/INFO tag, in the VCf, but this tag is not read correctly by SAIGE and probably other softwares
# we will add a R2 INFO tag, so we can use that information.
rule convertInfoToR2:
	output:
		output_folder+"/06.imputed/R2/{chr}/{chr}.vcf.gz",
		output_folder+"/06.imputed/R2/{chr}/{chr}.vcf.gz.csi",
		output_folder+"/06.imputed/R2/{chr}/{chr}.vcf.gz.tbi"
	input:
		output_folder+"/06.imputed/MERGED/{chr}/{chr}.vcf.gz"
	params:
		bcftools_bin=config['tools']['bcftools']
	log:
		stdout=log_folder+"/convertInfoToR2_{chr}.o",
		stderr=log_folder+"/convertInfoToR2_{chr}.e"
	shell:
		"""
		#We need to add the correct annotation name for the impute info score. We will call it R2, so other softwares can use this field (e.g. SAIGE and others)
		{params.bcftools_bin} annotate -c INFO/R2:=INFO/INFO -a {input} {input}| {params.bcftools_bin} view -O z -o {output[0]} 1> {log.stdout} 2> {log.stderr}
		{params.bcftools_bin} index -t {output[0]}
		{params.bcftools_bin} index -c {output[0]}
		"""

#add a rule to generate also a table similar to the info_score table from impute, with stats from the imputed vcf
rule imputeTabStat:
	output:
		output_folder+"/07.stats/{chr}/{chr}.info_stats.gz"
	input:
		output_folder+"/06.imputed/MERGED/{chr}/{chr}.vcf.gz"
	params:
		bcftools_bin=config['tools']['bcftools'],
	log:
		stdout=log_folder+ "imputeTabStat_{chr}.o"
		stderr=log_folder+ "imputeTabStat_{chr}.e"
	shell:
		"""
		{params.bcftools_bin} query -H -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%INFO/INFO\n" {input} | gzip -c > {output}
		"""