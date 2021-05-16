#Snakefile for the Internal Imputation pipeline
#
#Author: Massimiliano Cocca
#

#We assume all our data has been already strand oriented with plink
#We will orient the data to match the reference panel orientation, using shapeit
# Input files will be the plink genotypes. We will get them from a config file
###############################################################################################

# Split plink formatted input files by chromosome
rule plinkSplit:
    output:
        expand("{{output_folder}}/00.splitted_input/{{cohort_name}}_{chr}.{ext}", ext=['bed','bim','fam'],chr=chrs)
        # o_bed=scatter.split(output_folder+"/00.splitted_input/{scatteritem}_"+cohort_name+".bed"),
        # o_bim=scatter.split(output_folder+"/00.splitted_input/{scatteritem}_"+cohort_name+".bim"),
        # o_fam=scatter.split(output_folder+"/00.splitted_input/{scatteritem}_"+cohort_name+".fam")
    input:
        expand(input_prefix+".{ext}", ext=['map','ped'])
    params:
        # scatter_chr= lambda w, output : re.search('(\d+-of-\d+)',output[0]).group(1).split('-of-')[0] ,
        output_prefix=output_folder+"/00.splitted_input/"+cohort_name,
        i_prefix=input_prefix,
        plink=config['tools']['plink']
    # log:
    #     stdout=log_folder+"/plinkSplit_{scatteritem}.stdout",
    #     stderr=log_folder+"/plinkSplit_{scatteritem}.stderr"
    run:
        for chr in chrs:
            cmd="%s --file %s --chr %s --make-bed --out %s_%s" % (params.plink,params.i_prefix,chr,params.output_prefix,chr)
            shell(cmd)
        

rule snpCheck:
    output:
        # expand("{{output_folder}}/01.refAlign/{{ref_panel}}/{chr}_shapeit_refpanel.alignments.snp.{ext}", ext=['strand','strand.exclude'])
        output_folder+"/01.refAlign/"+ref_panel+"/{chr}_shapeit_"+ref_panel+".alignments.snp.strand",
        output_folder+"/01.refAlign/"+ref_panel+"/{chr}_shapeit_"+ref_panel+".alignments.snp.strand.exclude"
        # config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_shapeit_refpanel.alignments.snp.strand",
        # config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_shapeit_refpanel.alignments.snp.strand.exclude"
    input:
        ug_bed=output_folder + "/00.splitted_input/"+cohort_name+"_{chr}.bed",
        ug_bim=output_folder + "/00.splitted_input/"+cohort_name+"_{chr}.bim",
        ug_fam=output_folder + "/00.splitted_input/"+cohort_name+"_{chr}.fam",
        rp_hap=config["paths"]["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".hap.gz",
        rp_legend=config["paths"]["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".legend.gz",
        rp_samples=config["paths"]["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".samples",
        g_map="/netapp/nfs/resources/1000GP_phase3/impute/genetic_map_chr{chr}_combined_b37.txt",
    params:
        output_prefix=output_folder+"/01.refAlign/"+ref_panel+"/{chr}_shapeit_"+ref_panel+".alignments",
        shapeit=config['tools']['shapeit']
        # output_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_shapeit_refpanel.alignments"
    shell:
        # {config[shapeit_path]} -check --input-bed {input.ug_bed} {input.ug_bim} {input.ug_fam} \
        """
        set +e
        {params.shapeit} -check --input-bed {input.ug_bed} {input.ug_bim} {input.ug_fam} \
        -M {input.g_map} \
        --input-ref {input.rp_hap} {input.rp_legend} {input.rp_samples} \
        --output-log {params.output_prefix}
        exitcode=$?
        if [ $exitcode -eq 0 ]
        then
            echo "No error found..exiting correctly"
            exit 0
        else
            echo "WARNING....The software raised some errors or warning, be careful and check the results. (EXIT CODE ${{exitcode}})"
            exit 0
        fi
        """

# rule snp_flip_file:
#     input:
#         rules.snp_check.output[0],
#     output:
#         strand_rsid=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_rsids.to_flip"
#     params:
#         bfiles_prefix=config["input_folder"] + "/" + config["chr"] + "/" + config["chr"],
#         bfiles_flipped_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_" + config["chr"] +  "_flipped"
#     run:
#         get_flippable(input[0],output.strand_rsid)
    
# rule snp_flip:
#     input:
#         rules.snp_flip_file.output[0],
#         ug_bed=config["input_folder"] + "/" + config["chr"] + "/" + config["chr"]+ ".bed",
#         ug_bim=config["input_folder"] + "/" + config["chr"] + "/" + config["chr"]+ ".bim",
#         ug_fam=config["input_folder"] + "/" + config["chr"] + "/" + config["chr"]+ ".fam"
#     output:
#         config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_" + config["chr"]+ "_flipped.bim",
#         config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_" + config["chr"]+ "_flipped.bed",
#         config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_" + config["chr"]+ "_flipped.fam",
#         # strand_rsid=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_rsids.to_flip"
#     params:
#         bfiles_prefix=config["input_folder"] + "/" + config["chr"] + "/" + config["chr"],
#         bfiles_flipped_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_" + config["chr"] +  "_flipped"
#     shell:
#         """
#         set +e
#         #we need this file and it could be empty, so we will touch it!
#         # touch {input[0]}
#         plink --bfile {params.bfiles_prefix} --flip {input[0]} --make-bed --out {params.bfiles_flipped_prefix}
#         exitcode=$?
#         if [ $exitcode -eq 0 ]
#         then
#             echo "No error found..exiting correctly"
#             exit 0
#         else
#             echo "WARNING....The software raised some errors or warning, be careful and check the results. (EXIT CODE ${{exitcode}})"
#             exit 0
#         fi
#         """
