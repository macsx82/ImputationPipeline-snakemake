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
	resources:
		mem_mb=10000
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
	resources:
		mem_mb=10000
	log:
		stdout=log_folder+"/infoStatsChrom_{chr}.o",
		stderr=log_folder+"/infoStatsChrom_{chr}.e"
	priority: 49
	shell:
		"""
		{params.bcftools_bin} query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%INFO/INFO\n" {input[0]} | {params.scripts_folder}/imputationStats.py --tab {params.tab_prefix} --fig {output[2]} > {log.stdout} 2> {log.stderr}
		{params.bcftools_bin} query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%INFO/INFO\n" {input[0]} | {params.scripts_folder}/plot_manhattan.py --no-log --cols 0,1,6 --title 'Impute INFO score chr{wildcards.chr}' --image {output[3]} --ymax=1.2 - >> {log.stdout} 2>> {log.stderr}
		"""


#aggregator rule to get all data from all chunks and generate a single pdf file
rule pdfReportChunks:
	wildcard_constraints:
		g_chunk='\d+',
		chr='\d+'
	output:
		output_folder+"/07.stats/{chr}/{chr}_impute_summary_report_by_chunk.pdf"
	input:
		# chunk_stats_by_maf_by_info=expand(output_folder+"/07.stats/{{chr}}/CHUNKS/{{chr}}_{g_chunk}_impute_summary_by_maf_by_info.csv",zip,**glob_wildcards(os.path.join(output_folder+"/06.imputed/{chr}/", "{chr}.{g_chunk}.vcf.gz"))._asdict()),
		# chunk_stats_by_maf=expand(output_folder+"/07.stats/{{chr}}/CHUNKS/{{chr}}_{g_chunk}_impute_summary_by_maf.csv",zip,**glob_wildcards(os.path.join(output_folder+"/06.imputed/{chr}/", "{chr}.{g_chunk}.vcf.gz"))._asdict()),
		# info_af=expand(output_folder+"/07.stats/{{chr}}/CHUNKS/{{chr}}_{g_chunk}_impute_summary.png",zip,**glob_wildcards(os.path.join(output_folder+"/06.imputed/{chr}/", "{chr}.{g_chunk}.vcf.gz"))._asdict()),
		# manhattan=expand(output_folder+"/07.stats/{{chr}}/CHUNKS/{{chr}}_{g_chunk}_impute_manhattan.png",zip,**glob_wildcards(os.path.join(output_folder+"/06.imputed/{chr}/", "{chr}.{g_chunk}.vcf.gz"))._asdict())
		# chunk_stats_by_maf_by_info=lambda wildcards: expand(output_folder+"/07.stats/{chr}/CHUNKS/{chr}_{g_chunk}_impute_summary_by_maf_by_info.csv",chr=wildcards.chr,g_chunk=glob_wildcards(os.path.join(output_folder+"/06.imputed/"+wildcards.chr+"/", wildcards.chr+".{g_chunk}.vcf.gz")).g_chunk),
		# chunk_stats_by_maf=lambda wildcards: expand(output_folder+"/07.stats/{chr}/CHUNKS/{chr}_{g_chunk}_impute_summary_by_maf.csv",chr=wildcards.chr,g_chunk=glob_wildcards(os.path.join(output_folder+"/06.imputed/"+wildcards.chr+"/", wildcards.chr+".{g_chunk}.vcf.gz")).g_chunk),
		# info_af=lambda wildcards: expand(output_folder+"/07.stats/{chr}/CHUNKS/{chr}_{g_chunk}_impute_summary.png",chr=wildcards.chr,g_chunk=glob_wildcards(os.path.join(output_folder+"/06.imputed/"+wildcards.chr+"/", wildcards.chr+".{g_chunk}.vcf.gz")).g_chunk),
		# manhattan=lambda wildcards: expand(output_folder+"/07.stats/{chr}/CHUNKS/{chr}_{g_chunk}_impute_manhattan.png",chr=wildcards.chr,g_chunk=glob_wildcards(os.path.join(output_folder+"/06.imputed/"+wildcards.chr+"/", wildcards.chr+".{g_chunk}.vcf.gz")).g_chunk)
		chunk_stats_by_maf_by_info=lambda wildcards: expand(rules.infoStatsChunks.output[0],chr=wildcards.chr,g_chunk=glob_wildcards(os.path.join(output_folder+"/06.imputed/"+wildcards.chr+"/", wildcards.chr+".{g_chunk}.vcf.gz")).g_chunk),
		chunk_stats_by_maf=lambda wildcards: expand(rules.infoStatsChunks.output[1],chr=wildcards.chr,g_chunk=glob_wildcards(os.path.join(output_folder+"/06.imputed/"+wildcards.chr+"/", wildcards.chr+".{g_chunk}.vcf.gz")).g_chunk),
		info_af=lambda wildcards: expand(rules.infoStatsChunks.output[2],chr=wildcards.chr,g_chunk=glob_wildcards(os.path.join(output_folder+"/06.imputed/"+wildcards.chr+"/", wildcards.chr+".{g_chunk}.vcf.gz")).g_chunk),
		manhattan=lambda wildcards: expand(rules.infoStatsChunks.output[3],chr=wildcards.chr,g_chunk=glob_wildcards(os.path.join(output_folder+"/06.imputed/"+wildcards.chr+"/", wildcards.chr+".{g_chunk}.vcf.gz")).g_chunk)
	
	params:
		stat_base_folder=output_folder+"/07.stats/{chr}/CHUNKS"
	log:
		stdout=log_folder+"/pdfReportChunks_{chr}.o",
		stderr=log_folder+"/pdfReportChunks_{chr}.e"	
	priority: 1
	run:
		chunk_number=len(input.chunk_stats_by_maf_by_info)
		pdf_report(wildcards.chr,params.stat_base_folder,chunk_number,output[0])
