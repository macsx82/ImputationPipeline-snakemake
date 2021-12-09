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
	resources:
		mem_mb=10000
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
	resources:
		mem_mb=10000
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
		output_folder+"/06.imputed/R2/{chr}/{chr}.vcf.gz"
	params:
		bcftools_bin=config['tools']['bcftools']
	resources:
		mem_mb=5000
	log:
		stdout=log_folder+ "/imputeTabStat_{chr}.o",
		stderr=log_folder+ "/imputeTabStat_{chr}.e"
	shell:
		"""
		{params.bcftools_bin} query -H -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%INFO/R2\t%INFO/IMP\n" {input} | awk 'BEGIN{FS="\t";OFS="\t"}{if($8==".") print $1,$2,$3,$4,$5,$6,$7,0;else print $0}' | gzip -c > {output}
		"""

#add a release rule
rule release:
	output:
		directory(release_folder+"/08.release/00.CLEANED_INPUT"),
		directory(release_folder+"/08.release/01.FLIPPED_INPUT"),
		directory(release_folder+"/08.release/02.PHASED"),
		directory(release_folder+"/08.release/03.IMPUTED/VCF"),
		directory(release_folder+"/08.release/03.IMPUTED/BIMBAM"),
		directory(release_folder+"/08.release/04.STATS")
	input:
		cleaned_input_1=rules.indelsRemove.output,
		cleaned_input_2=rules.plinkSplit.output,
		flipped_input=expand(output_folder + "/03.flipped_input/" + ref_panel + "/VCF/"+ cohort_name+"_{chr}_fixRef_sorted_rsID.{ext}", chr=chrs, ext=["vcf.gz","vcf.gz.tbi"]),
		phased_input=expand(output_folder+ "/04.phased_data/" + ref_panel + "/"+ cohort_name +"_{chr}_phased.{ext}",chr=chrs, ext=["vcf.gz","vcf.gz.tbi"]),
		imputed_1=expand(output_folder+"/06.imputed/R2/{chr}/{chr}.{ext}", ext=["vcf.gz","vcf.gz.tbi","vcf.gz.csi"],chr=chrs),
		imputed_2=expand(output_folder+"/06.imputed/BIMBAM/{chr}/{chr}.{ext}",ext=["bimbam.gz","pos"],chr=chrs),
		stats=expand(output_folder+"/07.stats/{chr}/{chr}.info_stats.gz", chr=chrs)
	resources:
		mem_mb=3000
	log:
		stdout=log_folder+ "/data_release.o",
		stderr=log_folder+ "/data_release.e"
	shell:
		"""

		#1) cp cleaned input files, only the first and the last
		mkdir -p {output[0]}
		for o_file in {input.cleaned_input_1}
		do
			rsync -avP ${{o_file}} {output[0]}/.
		done 1> {log.stdout} 2> {log.stderr}

		for o_file in {input.cleaned_input_2}
		do
			rsync -avP ${{o_file}} {output[0]}/.
		done 1>> {log.stdout} 2>> {log.stderr}
		
		#2) cp flipped input files, only the rs annotated vcf files
		mkdir -p {output[1]}
		for o_file in {input.flipped_input}
		do
			rsync -avP ${{o_file}} {output[1]}/.
		done 1>> {log.stdout} 2>> {log.stderr}
		
		#3) cp phased data
		mkdir -p {output[2]}
		for o_file in {input.phased_input}
		do
			rsync -avP ${{o_file}} {output[2]}/.
		done 1>> {log.stdout} 2>> {log.stderr}
		
		#4) cp Imputed results: IMPUTE/R2 folder (that will be renamed as VCF) and BIMBAM formatted files
		mkdir -p {output[3]}
		for o_file in {input.imputed_1}
		do
			rsync -avP ${{o_file}} {output[3]}/.
		done 1>> {log.stdout} 2>> {log.stderr}
		
		mkdir -p {output[4]}
		for o_file in {input.imputed_2}
		do
			rsync -avP ${{o_file}} {output[4]}/.
		done 1>> {log.stdout} 2>> {log.stderr}
		
		#5) cp info stats in tab format
		mkdir -p {output[5]}
		for o_file in {input.stats}
		do
			rsync -avP ${{o_file}} {output[5]}/.
		done 1>> {log.stdout} 2>> {log.stderr}
		"""