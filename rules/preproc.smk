#Snakefile for the Internal Imputation pipeline
#
#Author: Massimiliano Cocca
#

#We assume all our data has been already strand oriented with plink
#We will orient the data to match the reference panel orientation, using shapeit
# Input files will be the plink genotypes. We will get them from a config file
###############################################################################################

# We want to remove indels from plink files, since we cannot update alleles in a consistent way, at the mment (will be a feature for next release)
rule indelsRemove:
    output:
        o_ped=temp(output_folder+"/00.splitted_input/"+cohort_name+"_snps_only.ped"),
        o_map=temp(output_folder+"/00.splitted_input/"+cohort_name+"_snps_only.map")
        # o_fam=scatter.split(output_folder+"/00.splitted_input/{scatteritem}_"+cohort_name+".fam")
    input:
        expand(input_prefix+".{ext}", ext=['map','ped'])
    params:
        # scatter_chr= lambda w, output : re.search('(\d+-of-\d+)',output[0]).group(1).split('-of-')[0] ,
        output_prefix=output_folder+"/00.splitted_input/"+cohort_name+"_snps_only",
        i_prefix=input_prefix,
        plink=config['tools']['plink']
    run:
        cmd="%s --file %s --snps-only 'just-acgt' --recode --out %s" % (params.plink,params.i_prefix,params.output_prefix)
        shell(cmd)

rule plinkSplit:
    output:
        expand(output_folder+"/00.splitted_input/"+cohort_name+"_{chr}.{ext}", ext=['bed','bim','fam'],chr=chrs)
        # o_bed=scatter.split(output_folder+"/00.splitted_input/{scatteritem}_"+cohort_name+".bed"),
        # o_bim=scatter.split(output_folder+"/00.splitted_input/{scatteritem}_"+cohort_name+".bim"),
        # o_fam=scatter.split(output_folder+"/00.splitted_input/{scatteritem}_"+cohort_name+".fam")
    input:
        rules.indelsRemove.output[0],
        rules.indelsRemove.output[1]
        # expand(input_prefix+".{ext}", ext=['map','ped'])
    params:
        # scatter_chr= lambda w, output : re.search('(\d+-of-\d+)',output[0]).group(1).split('-of-')[0] ,
        # i_prefix=input_prefix,
        output_prefix=output_folder+"/00.splitted_input/"+cohort_name,
        i_prefix=output_folder+"/00.splitted_input/"+cohort_name+"_snps_only",
        plink=config['tools']['plink']
    # log:
    #     stdout=log_folder+"/plinkSplit_{scatteritem}.stdout",
    #     stderr=log_folder+"/plinkSplit_{scatteritem}.stderr"
    run:
        for chr in chrs:
            cmd="%s --file %s --chr %s --make-bed --out %s_%s" % (params.plink,params.i_prefix,chr,params.output_prefix,chr)
            shell(cmd)


# align aleles to 1000G data and retrieve alleles names removed because monomorphic
rule allFix:
    output:
        output_folder + "/00.splitted_input/" + ref_panel + "/"+ cohort_name+"_{chr}_allFix.bim",
        output_folder + "/00.splitted_input/" + ref_panel + "/"+ cohort_name+"_{chr}_allFix.bed",
        output_folder + "/00.splitted_input/" + ref_panel + "/"+ cohort_name+"_{chr}_allFix.fam",
        temp(output_folder + "/00.splitted_input/"+cohort_name+"_{chr}_a1.bim"),
        temp(output_folder + "/00.splitted_input/"+cohort_name+"_{chr}_a1.bed"),
        temp(output_folder + "/00.splitted_input/"+cohort_name+"_{chr}_a1.fam")
    input:
        # rules.snpFlip.output[0],
        # rules.snpFlip.output[1],
        # rules.snpFlip.output[2]
        ug_bed=output_folder + "/00.splitted_input/"+cohort_name+"_{chr}.bed",
        ug_bim=output_folder + "/00.splitted_input/"+cohort_name+"_{chr}.bim",
        ug_fam=output_folder + "/00.splitted_input/"+cohort_name+"_{chr}.fam"
    params:
        bfiles_prefix=output_folder + "/00.splitted_input/"+cohort_name+"_{chr}",
        bfiles_prefix_a1=output_folder + "/00.splitted_input/"+cohort_name+"_{chr}_a1",
        bfiles_allFix_prefix=output_folder+"/00.splitted_input/"+ ref_panel + "/" + cohort_name+"_{chr}_allFix",
        plink=config['tools']['plink'],
        update_a1_str=config['paths']['allele_recode_file']+" 5 3 '#'",
        update_a2_str=config['paths']['allele_recode_file']+" 4 3 '#'"
    shell:
        """
        {params.plink} --bfile {params.bfiles_prefix} --a2-allele {params.update_a2_str} --make-bed --out {params.bfiles_prefix_a1}
        {params.plink} --bfile {params.bfiles_prefix_a1} --keep-allele-order --a1-allele {params.update_a1_str} --make-bed --out {params.bfiles_allFix_prefix}
        """

rule snpCheck:
    output:
        # expand("{{output_folder}}/01.refAlign/{{ref_panel}}/{chr}_shapeit_refpanel.alignments.snp.{ext}", ext=['strand','strand.exclude'])
        output_folder+"/01.refAlign/"+ref_panel+"/{chr}_shapeit_"+ref_panel+".alignments.snp.strand",
        output_folder+"/01.refAlign/"+ref_panel+"/{chr}_shapeit_"+ref_panel+".alignments.snp.strand.exclude"
        # config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_shapeit_refpanel.alignments.snp.strand",
        # config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_shapeit_refpanel.alignments.snp.strand.exclude"
    input:
        ug_bed=output_folder + "/00.splitted_input/" + ref_panel + "/"+ cohort_name+"_{chr}_allFix.bim",
        ug_bim=output_folder + "/00.splitted_input/" + ref_panel + "/"+ cohort_name+"_{chr}_allFix.bed",
        ug_fam=output_folder + "/00.splitted_input/" + ref_panel + "/"+ cohort_name+"_{chr}_allFix.fam",
        rp_hap=config["paths"]["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".hap.gz",
        rp_legend=config["paths"]["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".legend.gz",
        rp_samples=config["paths"]["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".samples",
        g_map=config['paths']['genetic_map_path']+"/genetic_map_chr{chr}_combined_b37.txt"
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

rule snpFlipFile:
    output:
        strand_rsid=output_folder+"/01.refAlign/"+ref_panel+"/{chr}_shapeit_rsids.to_flip"
    input:
        rules.snpCheck.output[0]
    run:
        get_flippable(input[0],output.strand_rsid)

# flip snps recovered from the previous file
rule snpFlip:
    output:
        output_folder + "/02.flipped_input/" + ref_panel + "/"+ cohort_name+"_{chr}_allFix_flipped.bim",
        output_folder + "/02.flipped_input/" + ref_panel + "/"+ cohort_name+"_{chr}_allFix_flipped.bed",
        output_folder + "/02.flipped_input/" + ref_panel + "/"+ cohort_name+"_{chr}_allFix_flipped.fam"
        # strand_rsid=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_rsids.to_flip"
    input:
        rules.snpFlipFile.output[0],
        # rules.plinkSplit.output[0],
        # rules.plinkSplit.output[1],
        # rules.plinkSplit.output[2]
        rules.allFix.output[0],
        rules.allFix.output[1],
        rules.allFix.output[2]
    params:
        bfiles_prefix=output_folder+"/00.splitted_input/"+ ref_panel + "/" + cohort_name+"_{chr}_allFix",
        bfiles_flipped_prefix=output_folder+"/02.flipped_input/"+ ref_panel + "/" + cohort_name+"_{chr}_flipped",
        plink=config['tools']['plink']
    shell:
        """
        set +e
        #we need this file and it could be empty, so we will touch it!
        # touch {input[0]}
        {params.plink} --bfile {params.bfiles_prefix} --flip {input[0]} --keep-allele-order --make-bed --out {params.bfiles_flipped_prefix}
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

# convert to vcf file format to use SHAPEIT4
rule plink2vcf:
    output:
        output_folder + "/02.flipped_input/" + ref_panel + "/"+ cohort_name+"_{chr}_flipped.vcf.gz",
        output_folder + "/02.flipped_input/" + ref_panel + "/"+ cohort_name+"_{chr}_flipped.vcf.gz.tbi"
        # strand_rsid=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_rsids.to_flip"
    input:
        rules.snpFlip.output[0],
        rules.snpFlip.output[1],
        rules.snpFlip.output[2]
    params:
        bfiles_allFix_prefix=output_folder + "/02.flipped_input/" + ref_panel + "/"+ cohort_name+"_{chr}_allFix_flipped",
        vcf_flipped_prefix=output_folder + "/02.flipped_input/" + ref_panel + "/"+ cohort_name+"_{chr}_allFix_flipped",
        plink=config['tools']['plink']
    shell:
        """
        set +e
        #we need this file and it could be empty, so we will touch it!
        {params.plink} --bfile {params.bfiles_allFix_prefix} --recode vcf-iid bgz --out {params.vcf_flipped_prefix}
        tabix -p vcf {output[0]}
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
