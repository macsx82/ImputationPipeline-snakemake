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
        o_ped=output_folder+"/00.cleaned_input/"+cohort_name+"_snps_only.ped",
        o_map=output_folder+"/00.cleaned_input/"+cohort_name+"_snps_only.map"
        # o_fam=scatter.split(output_folder+"/00.splitted_input/{scatteritem}_"+cohort_name+".fam")
    input:
        expand(input_prefix+".{ext}", ext=['map','ped'])
    params:
        # scatter_chr= lambda w, output : re.search('(\d+-of-\d+)',output[0]).group(1).split('-of-')[0] ,
        output_prefix=output_folder+"/00.cleaned_input/"+cohort_name+"_snps_only",
        i_prefix=input_prefix,
        plink=config['tools']['plink']
    log:
        stdout=log_folder+"/indelsRemove.o",
        stderr=log_folder+"/indelsRemove.e"
    run:
        cmd="%s --file %s --snps-only 'just-acgt' --recode --out %s" % (params.plink,params.i_prefix,params.output_prefix)
        shell(cmd)

# Update allele positions and chromosomes using 1000G data. We need to apply this one BEFORE splitting the data by chr, since we will have chr and positions moving around
rule mapUpdateExt:
    output:
        o_bim=output_folder+"/00.cleaned_input/"+cohort_name+"_snps_only_mapUpdateExt.bim",
        o_bed=output_folder+"/00.cleaned_input/"+cohort_name+"_snps_only_mapUpdateExt.bed",
        o_fam=output_folder+"/00.cleaned_input/"+cohort_name+"_snps_only_mapUpdateExt.fam"
    input:
        rules.indelsRemove.output[0],
        rules.indelsRemove.output[1]
    params:
        i_prefix=output_folder+"/00.cleaned_input/"+cohort_name+"_snps_only",
        bfiles_out_prefix=output_folder+"/00.cleaned_input/"+cohort_name+"_snps_only_mapUpdateExt",
        plink=config['tools']['plink'],
        update_chr_str=config['paths']['allele_recode_file']+" 1 3 '#'",
        update_map_str=config['paths']['allele_recode_file']+" 2 3 '#'"
    log:
        stdout=log_folder+"/mapUpdateExt.o",
        stderr=log_folder+"/mapUpdateExt.e"
    shell:
        """
        {params.plink} --file {params.i_prefix} --update-chr {params.update_chr_str} --update-map {params.update_map_str} --make-bed --out {params.bfiles_out_prefix}
        """

# align alleles to external reference data and retrieve alleles names removed because monomorphic
rule allFix:
    output:
        output_folder+"/00.cleaned_input/"+ cohort_name+"_snps_only_mapUpdateExt_allFix.bim",
        output_folder+"/00.cleaned_input/"+ cohort_name+"_snps_only_mapUpdateExt_allFix.bed",
        output_folder+"/00.cleaned_input/"+ cohort_name+"_snps_only_mapUpdateExt_allFix.fam",
        output_folder+"/00.cleaned_input/"+ cohort_name+"_snps_only_mapUpdateExt_allFix.log",
        output_folder + "/00.cleaned_input/"+cohort_name+"_snps_only_mapUpdateExt_a1.bim",
        output_folder + "/00.cleaned_input/"+cohort_name+"_snps_only_mapUpdateExt_a1.bed",
        output_folder + "/00.cleaned_input/"+cohort_name+"_snps_only_mapUpdateExt_a1.fam",
        output_folder + "/00.cleaned_input/"+cohort_name+"_snps_only_mapUpdateExt_a1.log"
    input:
        rules.mapUpdateExt.output[0],
        rules.mapUpdateExt.output[1],
        rules.mapUpdateExt.output[2]
    params:
        bfiles_prefix=output_folder+"/00.cleaned_input/"+cohort_name+"_snps_only_mapUpdateExt",
        bfiles_prefix_a1=output_folder + "/00.cleaned_input/"+cohort_name+"_snps_only_mapUpdateExt_a1",
        bfiles_allFix_prefix=output_folder+"/00.cleaned_input/"+cohort_name+"_snps_only_mapUpdateExt_allFix",
        plink=config['tools']['plink'],
        update_a1_str=config['paths']['allele_recode_file']+" 5 3 '#'",
        update_a2_str=config['paths']['allele_recode_file']+" 4 3 '#'"
    log:
        stdout=log_folder+"/allFix.o",
        stderr=log_folder+"/allFix.e"
    shell:
        """
        {params.plink} --bfile {params.bfiles_prefix} --a1-allele {params.update_a1_str} --make-bed --out {params.bfiles_prefix_a1}
        {params.plink} --bfile {params.bfiles_prefix_a1} --keep-allele-order --a2-allele {params.update_a2_str} --make-bed --out {params.bfiles_allFix_prefix}
        """

#generate a file collecting all snps raising an impossible allele assignment warning
rule impossibleAllAssignmentFile:
    output:
        to_flip_rsid=output_folder+"/00.cleaned_input/"+cohort_name+"_impossibleAllAssignmentFile_rsids.to_flip"
    input:
        rules.allFix.output[3],
        rules.allFix.output[7]
    log:
        stdout=log_folder+"/impossibleAllAssignmentFile.o",
        stderr=log_folder+"/impossibleAllAssignmentFile.e"
    run:
        get_not_assigned_snps(output.to_flip_rsid,input[0],input[1])

# flip snps recovered from the previous rule
rule allFixSnpFlip:
    output:
        output_folder+"/00.cleaned_input/"+ cohort_name+"_snps_only_mapUpdateExt_flipped.bim",
        output_folder+"/00.cleaned_input/"+ cohort_name+"_snps_only_mapUpdateExt_flipped.bed",
        output_folder+"/00.cleaned_input/"+ cohort_name+"_snps_only_mapUpdateExt_flipped.fam"
    input:
        rules.impossibleAllAssignmentFile.output[0],
        rules.allFix.output[0],
        rules.allFix.output[1],
        rules.allFix.output[2]
    params:
        bfiles_prefix=output_folder+"/00.cleaned_input/"+cohort_name+"_snps_only_mapUpdateExt",
        bfiles_flipped_prefix=output_folder+"/00.cleaned_input/"+ cohort_name+"_snps_only_mapUpdateExt_flipped",
        plink=config['tools']['plink']
    log:
        stdout=log_folder+"/allFixSnpFlip.o",
        stderr=log_folder+"/allFixSnpFlip.e"
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

# After the flipping, we can continue as before, performing the allele fixing on the splitted dataset
rule plinkSplit:
    output:
        expand(output_folder+"/01.splitted_input/"+cohort_name+"_{chr}.{ext}", ext=['bed','bim','fam'],chr=chrs)
    input:
        rules.allFixSnpFlip.output[0],
        rules.allFixSnpFlip.output[1],
        rules.allFixSnpFlip.output[2]
    params:
        output_prefix=output_folder+"/01.splitted_input/"+cohort_name,
        i_prefix=output_folder+"/00.cleaned_input/"+ cohort_name+"_snps_only_mapUpdateExt_flipped",
        plink=config['tools']['plink']
    log:
        stdout=log_folder+"/plinkSplit.o",
        stderr=log_folder+"/plinkSplit.e"
    run:
        for chr in chrs:
            cmd="%s --bfile %s --chr %s --make-bed --out %s_%s" % (params.plink,params.i_prefix,chr,params.output_prefix,chr)
            shell(cmd)


# Now run the second step of allele alignment to external ref data. Now we should have no warnings at all, or the ones still present are from unrecoverable snps
rule allFixSplitted:
    output:
        output_folder + "/01.splitted_input/" + ref_panel + "/"+ cohort_name+"_{chr}_allFix.bim",
        output_folder + "/01.splitted_input/" + ref_panel + "/"+ cohort_name+"_{chr}_allFix.bed",
        output_folder + "/01.splitted_input/" + ref_panel + "/"+ cohort_name+"_{chr}_allFix.fam",
        output_folder + "/01.splitted_input/"+ ref_panel + "/"+cohort_name+"_{chr}_a1.bim",
        output_folder + "/01.splitted_input/"+ ref_panel + "/"+cohort_name+"_{chr}_a1.bed",
        output_folder + "/01.splitted_input/"+ ref_panel + "/"+cohort_name+"_{chr}_a1.fam"
    input:
        rules.plinkSplit.output[0],
        rules.plinkSplit.output[1],
        rules.plinkSplit.output[2]
    params:
        # bfiles_prefix=output_folder+"/00.splitted_input/"+ ref_panel + "/" + cohort_name+"_{chr}_mapUpdateTGP",
        bfiles_prefix=output_folder + "/01.splitted_input/"+cohort_name+"_{chr}",
        bfiles_prefix_a1=output_folder + "/01.splitted_input/"+ ref_panel + "/"+cohort_name+"_{chr}_a1",
        bfiles_allFix_prefix=output_folder+"/01.splitted_input/"+ ref_panel + "/" + cohort_name+"_{chr}_allFix",
        plink=config['tools']['plink'],
        update_a1_str=config['paths']['allele_recode_file']+" 5 3 '#'",
        update_a2_str=config['paths']['allele_recode_file']+" 4 3 '#'"
    log:
        stdout=log_folder+"/allFixSplitted_{chr}.o",
        stderr=log_folder+"/allFixSplitted_{chr}.e"
    shell:
        """
        {params.plink} --bfile {params.bfiles_prefix} --a1-allele {params.update_a1_str} --make-bed --out {params.bfiles_prefix_a1}
        {params.plink} --bfile {params.bfiles_prefix_a1} --keep-allele-order --a2-allele {params.update_a2_str} --make-bed --out {params.bfiles_allFix_prefix}
        """

#check alele orientation with the genetic map provided
rule snpCheck:
    output:
        # expand("{{output_folder}}/01.refAlign/{{ref_panel}}/{chr}_shapeit_refpanel.alignments.snp.{ext}", ext=['strand','strand.exclude'])
        output_folder+"/02.refAlign/"+ref_panel+"/{chr}_shapeit_"+ref_panel+".alignments.snp.strand",
        output_folder+"/02.refAlign/"+ref_panel+"/{chr}_shapeit_"+ref_panel+".alignments.snp.strand.exclude"
        # config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_shapeit_refpanel.alignments.snp.strand",
        # config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_shapeit_refpanel.alignments.snp.strand.exclude"
    input:
        ug_bed=output_folder + "/01.splitted_input/" + ref_panel + "/"+ cohort_name+"_{chr}_allFix.bed",
        ug_bim=output_folder + "/01.splitted_input/" + ref_panel + "/"+ cohort_name+"_{chr}_allFix.bim",
        ug_fam=output_folder + "/01.splitted_input/" + ref_panel + "/"+ cohort_name+"_{chr}_allFix.fam",
        rp_hap=config["rules"]['snpCheck']["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".hap.gz",
        rp_legend=config["rules"]['snpCheck']["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".legend.gz",
        rp_samples=config["rules"]['snpCheck']["ref_panel_base_folder"]+ "/"+ref_panel+"/{chr}/{chr}."+ ref_panel+".samples",
        g_map=config['rules']['snpCheck']['genetic_map_path']+"/genetic_map_chr{chr}_combined_b37.txt"
    params:
        output_prefix=output_folder+"/02.refAlign/"+ref_panel+"/{chr}_shapeit_"+ref_panel+".alignments",
        shapeit=config['tools']['shapeit']
    log:
        stdout=log_folder+"/snpCheck_{chr}.o",
        stderr=log_folder+"/snpCheck_{chr}.e"
    shell:
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
        strand_rsid=output_folder+"/02.refAlign/"+ref_panel+"/{chr}_shapeit_rsids.to_flip"
    input:
        rules.snpCheck.output[0]
    log:
        stdout=log_folder+"/snpFlipFile_{chr}.o",
        stderr=log_folder+"/snpFlipFile_{chr}.e"
    run:
        get_flippable(input[0],output.strand_rsid)

# flip snps recovered from the previous file
rule snpFlip:
    output:
        output_folder + "/03.flipped_input/" + ref_panel + "/"+ cohort_name+"_{chr}_allFix_flipped.bim",
        output_folder + "/03.flipped_input/" + ref_panel + "/"+ cohort_name+"_{chr}_allFix_flipped.bed",
        output_folder + "/03.flipped_input/" + ref_panel + "/"+ cohort_name+"_{chr}_allFix_flipped.fam"
    input:
        rules.snpFlipFile.output[0],
        rules.allFix.output[0],
        rules.allFix.output[1],
        rules.allFix.output[2]
    params:
        bfiles_prefix=output_folder+"/01.splitted_input/"+ ref_panel + "/" + cohort_name+"_{chr}_allFix",
        bfiles_flipped_prefix=output_folder+"/03.flipped_input/"+ ref_panel + "/" + cohort_name+"_{chr}_allFix_flipped",
        plink=config['tools']['plink']
    log:
        stdout=log_folder+"/snpFlip_{chr}.o",
        stderr=log_folder+"/snpFlip_{chr}.e"
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

# rule to fix the missing alt alleles in plink files before vcf conversion
#here we just need to update the BIM file, no need to generate a new genotype set using plink
#to keep things organised we still will create a copy of the genotypes in the monoRecovered folder
rule recoverMono:
    output:
        output_folder + "/03.flipped_input/" + ref_panel + "/ReMo/"+ cohort_name+"_{chr}_allFix_flipped_ReMo.bim",
        output_folder + "/03.flipped_input/" + ref_panel + "/ReMo/"+ cohort_name+"_{chr}_allFix_flipped_ReMo.bed",
        output_folder + "/03.flipped_input/" + ref_panel + "/ReMo/"+ cohort_name+"_{chr}_allFix_flipped_ReMo.fam"
    input:
        chip_update_allele_file=config['paths']['snp_array_update_allele_file'],
        bim_file=rules.snpFlip.output[0],
        bed_file=rules.snpFlip.output[1],
        fam_file=rules.snpFlip.output[2]
    params:
        bfiles_allFix_prefix=output_folder + "/03.flipped_input/" + ref_panel + "/"+ cohort_name+"_{chr}_allFix_flipped",
        vcf_flipped_prefix=output_folder + "/03.flipped_input/" + ref_panel + "/VCF/"+ cohort_name+"_{chr}_allFix_flipped",
        plink=config['tools']['plink']
    log:
        stdout=log_folder+"/recoverMono_{chr}.o",
        stderr=log_folder+"/recoverMono_{chr}.e"
    run:
        update_mono_snps(input.chip_update_allele_file,input.bim_file,output[0])
        cp_bed_cmd="cp %s %s" %(input.bed_file,output[1])
        cp_fam_cmd="cp %s %s" %(input.fam_file,output[2])
        shell(cp_bed_cmd)
        shell(cp_fam_cmd)

# convert to vcf file format to use SHAPEIT4
rule plink2vcf:
    output:
        output_folder + "/03.flipped_input/" + ref_panel + "/VCF/"+ cohort_name+"_{chr}_allFix_flipped.vcf.gz",
        output_folder + "/03.flipped_input/" + ref_panel + "/VCF/"+ cohort_name+"_{chr}_allFix_flipped.vcf.gz.tbi"
    input:
        rules.recoverMono.output[0],
        rules.recoverMono.output[1],
        rules.recoverMono.output[2]
    params:
        bfiles_allFix_prefix=output_folder + "/03.flipped_input/" + ref_panel + "/ReMo/"+ cohort_name+"_{chr}_allFix_flipped_ReMo",
        vcf_flipped_prefix=output_folder + "/03.flipped_input/" + ref_panel + "/VCF/"+ cohort_name+"_{chr}_allFix_flipped",
        plink=config['tools']['plink']
    log:
        stdout=log_folder+"/plink2vcf_{chr}.o",
        stderr=log_folder+"/plink2vcf_{chr}.e"
    shell:
        """
        set +e
        #we need this file and it could be empty, so we will touch it!
        {params.plink} --bfile {params.bfiles_allFix_prefix} --keep-allele-order --recode vcf-iid bgz --out {params.vcf_flipped_prefix}
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

# fix strand according to external resource. We will be working by rsID, so we hope we can fix most of the issues, if still present.
rule vcfFixRef:
    output:
        output_folder + "/03.flipped_input/" + ref_panel + "/VCF/"+ cohort_name+"_{chr}_fixRef.vcf.gz",
        output_folder + "/03.flipped_input/" + ref_panel + "/VCF/"+ cohort_name+"_{chr}_fixRef_sorted.vcf.gz",
        output_folder + "/03.flipped_input/" + ref_panel + "/VCF/"+ cohort_name+"_{chr}_fixRef_sorted.vcf.gz.tbi"
    input:
        rules.plink2vcf.output[0],
        rules.plink2vcf.output[1]
    params:
        ext_ref_file=config['paths']['ext_ref_annot_file'],
        ref_fasta=config['paths']['ref_fasta'],
        temp=config['rules']['vcfFixRef']['temp'],
        bcftools_bin=config['tools']['bcftools']
    log:
        stdout=log_folder+"/vcfFixRef_{chr}.o",
        stderr=log_folder+"/vcfFixRef_{chr}.e",
    shell:
        """
        # conda activate gatk4170
        {params.bcftools_bin} +fixref {input[0]} -O z -o {output[0]} -- -i {params.ext_ref_file} -f {params.ref_fasta} > {log.stdout} 2> {log.stderr}
        {params.bcftools_bin} sort -T {params.temp} {output[0]} -O z -o {output[1]} >> {log.stdo} 2>> {log.stde}
        tabix -p vcf {output[1]}
        """

# Annotate with rsIds according to external resource
rule vcfAnnotate:
    output:
        output_folder + "/03.flipped_input/" + ref_panel + "/VCF/"+ cohort_name+"_{chr}_fixRef_sorted_pre_rsID.vcf.gz",
        output_folder + "/03.flipped_input/" + ref_panel + "/VCF/"+ cohort_name+"_{chr}_fixRef_sorted_pre_rsID.vcf.gz.tbi",
        output_folder + "/03.flipped_input/" + ref_panel + "/VCF/"+ cohort_name+"_{chr}_fixRef_sorted_rsID.vcf.gz",
        output_folder + "/03.flipped_input/" + ref_panel + "/VCF/"+ cohort_name+"_{chr}_fixRef_sorted_rsID.vcf.gz.tbi"
    input:
        rules.vcfFixRef.output[1],
        rules.vcfFixRef.output[2]
    params:
        ext_ref_file=config['paths']['ext_ref_annot_file'],
        bcftools_bin=config['tools']['bcftools']
    log:
        stdout=log_folder+"/vcfAnnotate_{chr}.o",
        stderr=log_folder+"/vcfAnnotate_{chr}.e"
    shell:
        """
        {params.bcftools_bin} annotate --set-id '%CHROM:%POS\_%REF\_%FIRST_ALT' {input[0]} -O z -o {output[0]} > {log.stdout} 2> {log.stderr}
        tabix -p vcf {output[0]}
        {params.bcftools_bin} annotate -a {params.ext_ref_file} -c ID {output[0]} -O z -o {output[2]} >> {log.stdout} 2>> {log.stderr}
        tabix -p vcf {output[2]}
        """

