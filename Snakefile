#Snakefile for the Internal Imputation pipeline
#
#Author: Massimiliano Cocca
#
#
# configfile: "config.yaml"
import pathlib
import logging


# recover some fixed variables from config file
output_folder = config['paths']["output_folder"]
log_folder = config['paths']["log_folder"]
cohort_name = config["cohort_name"]
input_prefix = config['paths']["input_file_prefix"]
chrs = config['chromosomes']
ref_panel=config['ref_panel']
ref_panel_base_folder=config["paths"]["ref_panel_base_folder"]
# chunk_size=config['rules']['impute']['chunk_size']
# define a scatter gather rule to work by chromosome
# CHR_COUNT=23

# scattergather:
#     split=CHR_COUNT
# MODULES
include_prefix="rules"
include:
    include_prefix + "/functions.py"

# define a dictionary containing chromosomes and relative chuk size
# chunked={}
# for chrom in chrs:
#     legend_file="%s/%s/%s/%s.%s.legend.gz" % (ref_panel_base_folder,ref_panel,chrom,chrom,ref_panel)
#     chunked[chrom]=get_chunk_by_chr(chrom,legend_file,chunk_size)
# print(chunked)
#define parameter useful to cluster job submission
localrules: all

#define rules
rule all:
    wildcard_constraints:
        g_chunk='\d+',
        chr='\d+'
    input:
        # config["output_folder"]+"/"+config["pop"]+"/" + config["pop"] + ".pipe.done"
        # lambda wildcards: config["chr"][wildcards.chrom],
        # expand(config["output_folder"]+"/"+config["pop"]+"/{chrom}.pipe.done", chrom=config["chr"])
        # config["output_folder"]+"/"+config["pop"]+"/" + config["chr"] + ".pipe.done"
        # expand(output_folder+"/00.splitted_input/"+cohort_name+"_{chr}.{ext}", ext=['bed','bim','fam'],chr=chrs)
        # expand(output_folder+"/01.refAlign/"+ref_panel+"/{chr}_shapeit_"+ref_panel+".alignments.snp.{ext}", ext=['strand','strand.exclude'],chr=chrs)
        # expand(output_folder+"/01.refAlign/"+ref_panel+"/{chr}_shapeit_rsids.to_flip",chr=chrs)
        # expand(output_folder + "/02.flipped_input/" + ref_panel + "/"+ cohort_name+"_{chr}_flipped.{ext}",ext=['bed','bim','fam'],chr=chrs)
        # expand(output_folder+ "/03.phased_data/" + ref_panel + "/chr{chr}.{ext}" , ext=["haps.gz","sample"], chr=chrs)
        # expand(output_folder+"/04.impute_intervals/{chr}/{chr}.{g_chunk}.int",chr=chrs,g_chunk=list(range(1,11)))
        # expand(output_folder+"/04.impute_intervals/{chrom}/{chrom}.{{g_chunk}}.pippo",chrom=chrs)
        # [ "%s/04.impute_intervals/%s/%s.%s.int" % (output_folder,key,key,value) for key, value in chunked.items()]
        # expand(output_folder + "/03.flipped_input/" + ref_panel + "/VCF/"+ cohort_name+"_{chr}_allFix_flipped.vcf.gz", chr=chrs),
        # expand(output_folder + "/03.flipped_input/" + ref_panel + "/VCF/"+ cohort_name+"_{chr}_allFix_flipped.vcf.gz.tbi", chr=chrs)
        # expand(output_folder + "/04.phased_data/" + ref_panel + "/"+ cohort_name +"_{chr}_phased.vcf.gz", chr=chrs),
        # expand(output_folder + "/04.phased_data/" + ref_panel + "/"+ cohort_name +"_{chr}_phased.vcf.gz.tbi", chr=chrs)
        # expand(output_folder+"/05.impute_intervals/{chr}/{chr}.{g_chunk}.int",chr=chrs,g_chunk=list(range(1,11)))
        # expand(output_folder+"/05.impute_intervals/{chr}/{chr}.{g_chunk}.int",chr=chrs,g_chunk=["{:02d}".format(chunk) for chunk in list(range(1,11))])
        # expand(output_folder+"/05.impute_intervals/{chr}/{chr}.{{g_chunk}}.int",chr=chrs)
        # output_folder+"/05.impute_intervals/{chr}/{chr}.{g_chunk}.int"
        # [ expand(output_folder+"/05.imputed/{chr}/{chr}.{g_chunk}.{ext}", ext=["gen.gz","gen_info","gen_info_by_sample","gen_samples","gen_summary","gen_warnings"], chr=key, g_chunk=["{:02d}".format(chunk) for chunk in list(range(1,value+1))]) for key, value in chunked.items()] 
        # [ expand(output_folder+"/05.imputed/{chr}/{chr}.{g_chunk}.{ext}", ext=["vcf.gz","log"], chr=key, g_chunk=["{:02d}".format(chunk) for chunk in list(range(1,value+1))]) for key, value in chunked.items()] 
        # [ expand(output_folder+"/06.imputed/{chr}/{chr}.{g_chunk}.{ext}", ext=["vcf.gz","log"], chr=chrs, g_chunk=["{:02d}".format(chunk) for chunk in list(range(1,5))])] 
        # expand(output_folder + "/03.flipped_input/" + ref_panel + "/VCF/"+ cohort_name+"_{chr}_fixRef_sorted_rsID.vcf.gz", chr=chrs),
        # expand(output_folder + "/03.flipped_input/" + ref_panel + "/VCF/"+ cohort_name+"_{chr}_fixRef_sorted_rsID.vcf.gz.tbi", chr=chrs)
        expand(output_folder+"/06.imputed/MERGED/{chr}/{chr}.{ext}", ext=["vcf.gz","vcf.gz.tbi"],chr=chrs),
        expand(output_folder+"/06.imputed/BIMBAM/{chr}/{chr}.{ext}", ext=["bimbam.gz","pos"],chr=chrs),
        expand(output_folder+"/07.stats/{chr}/{chr}_{ext}",chr=chrs,ext=['impute_summary_by_maf_by_info.csv','impute_summary_by_maf.csv','impute_summary.pdf','impute_manhattan.png']),
        # expand(output_folder+"/07.stats/{chr}/CHUNKS/{chr}_{g_chunk}_{ext}",chr=chrs,ext=['impute_summary_by_maf_by_info.csv','impute_summary_by_maf.csv','impute_summary.pdf','impute_manhattan.png'],g_chunk=["{:02d}".format(chunk) for chunk in list(range(1,lambda wildcards : getChunkNumByChr(wildcards)))])
        # expand(output_folder+"/07.stats/{chr}/CHUNKS/{chr}_{g_chunk}_{ext}",chr=chrs,ext=['impute_summary_by_maf_by_info.csv','impute_summary_by_maf.csv','impute_summary.pdf','impute_manhattan.png'],g_chunk=glob_wildcards(os.path.join(output_folder+"/06.imputed/{chr}/", "{chr}.{g_chunk}.vcf.gz")).g_chunk)
        expand(output_folder+"/07.stats/{chr}/CHUNKS/{chr}_{g_chunk}_impute_summary_by_maf_by_info.csv",zip,**glob_wildcards(os.path.join(output_folder+"/06.imputed/{chr1}/", "{chr}.{g_chunk}.vcf.gz"))._asdict()),
        expand(output_folder+"/07.stats/{chr}/CHUNKS/{chr}_{g_chunk}_impute_summary_by_maf.csv",zip,**glob_wildcards(os.path.join(output_folder+"/06.imputed/{chr1}/", "{chr}.{g_chunk}.vcf.gz"))._asdict()),
        expand(output_folder+"/07.stats/{chr}/CHUNKS/{chr}_{g_chunk}_impute_summary.pdf",zip,**glob_wildcards(os.path.join(output_folder+"/06.imputed/{chr1}/", "{chr}.{g_chunk}.vcf.gz"))._asdict()),
        expand(output_folder+"/07.stats/{chr}/CHUNKS/{chr}_{g_chunk}_impute_manhattan.png",zip,**glob_wildcards(os.path.join(output_folder+"/06.imputed/{chr1}/", "{chr}.{g_chunk}.vcf.gz"))._asdict())

        # expand(output_folder+"/07.stats/{chr}/CHUNKS/{chr}_{{g_chunk}}_{ext}",chr=chrs,ext=['impute_summary.csv','impute_summary.pdf','impute_manhattan.pdf'])

        # [ output_folder+"/04.impute_intervals/{key}/{key}.{value}.int" for key, value in chunked.items()]
        # directory(expand(output_folder+"/04.impute_intervals/{chr}/",chr=chrs))
        # expand(output_folder+"/04.impute_intervals/{chr}/{chr}.{{g_chunk}}.pippo",chr=chrs)

# MODULES
# include_prefix="rules"
# include:
#     include_prefix + "/functions.py"
# include:
    # include_prefix + "/preproc.smk"
include:
    include_prefix + "/phasing.smk"
include:
    include_prefix + "/impute.smk"
include:
    include_prefix + "/postproc.smk"


onsuccess:
    print("The workflow finished without errors!")

onerror:
    print("An error occurred in the current workflow execution!!")