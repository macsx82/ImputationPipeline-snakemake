#Snakefile for the Internal Imputation pipeline
#
#Author: Massimiliano Cocca
#
#
# configfile: "config.yaml"

def generate_shapeit_out_files(key):
    # config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/"
    chr_phased= "%s/%s/%s/%s/chr%s.haps.gz" % (config["output_folder"],config["pop"],config["ref_panel"],key,key)
    samples= "%s/%s/%s/%s/chr%s.samples" % (config["output_folder"],config["pop"],config["ref_panel"],key,key)

    return chr_phased,samples

def generate_end_of_pipeline_files(key):
    return "%s/%s/chr%s.pipe.done" % (config["output_folder"],config["pop"],key)

#define parameter useful to cluster job submission
localrules: all

#define rules
rule all:
    input:
        expand(config["output_folder"]+"/"+config["pop"]+"/{chrom}.pipe.done", chrom=config["chr"])
        

#We assume all our data has been already strand oriented with plink
#We will orient the data to match the reference panel orientation, using shapeit
# Input files will be the plink genotypes. We will get them from a config file
rule snp_check:
    input:
        # ug_bed=config["input_folder"] + "/" + config["chr"] + "/" + config["chr"]+ ".bed",
        # ug_bim=config["input_folder"] + "/" + config["chr"] + "/" + config["chr"]+ ".bim",
        # ug_fam=config["input_folder"] + "/" + config["chr"] + "/" + config["chr"]+ ".fam",
        # rp_hap=config["ref_panel_base_folder"]+ "/" + config["ref_panel"]+ "/" + config["chr"]+ "/" + config["chr"] + "." + config["ref_panel"] + ".hap.gz",
        # rp_legend=config["ref_panel_base_folder"]+ "/" + config["ref_panel"]+ "/" + config["chr"]+ "/" + config["chr"] + "." + config["ref_panel"] + ".legend.gz",
        # rp_samples=config["ref_panel_base_folder"]+ "/" + config["ref_panel"]+ "/" + config["chr"]+ "/" + config["chr"] + "." + config["ref_panel"] + ".samples"
        ug_bed=expand(config["input_folder"] + "/{chrom}/{chrom}.bed",chrom=config["chr"]),
        ug_bim=expand(config["input_folder"] + "/{chrom}/{chrom}.bim",chrom=config["chr"]),
        ug_fam=expand(config["input_folder"] + "/{chrom}/{chrom}.fam",chrom=config["chr"]),
        rp_hap=expand(config["ref_panel_base_folder"]+ "/" + config["ref_panel"]+ "/{chrom}/{chrom}." + config["ref_panel"] + ".hap.gz",chrom=config["chr"]),
        rp_legend=expand(config["ref_panel_base_folder"]+ "/" + config["ref_panel"]+ "/{chrom}/{chrom}." + config["ref_panel"] + ".legend.gz",chrom=config["chr"]),
        rp_samples=expand(config["ref_panel_base_folder"]+ "/" + config["ref_panel"]+ "/{chrom}/{chrom}." + config["ref_panel"] + ".samples",chrom=config["chr"])
    params:
        # g_map="/netapp/nfs/resources/1000GP_phase3/impute/genetic_map_chr"+config["chr"]+"_combined_b37.txt",
        # output_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_shapeit_refpanel.alignments"
        g_map=expand("/netapp/nfs/resources/1000GP_phase3/impute/genetic_map_chr{chrom}_combined_b37.txt",chrom=config["chr"]),
        output_prefix=expand(config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/{chrom}/" + config["pop"] + "_shapeit_refpanel.alignments",chrom=config["chr"])
    output:
        # config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_shapeit_refpanel.alignments.snp.strand",
        # config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_shapeit_refpanel.alignments.snp.strand.exclude"
        expand(config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/{chrom}/" + config["pop"] + "_shapeit_refpanel.alignments.snp.strand", chrom=config["chr"]),
        expand(config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/{chrom}/" + config["pop"] + "_shapeit_refpanel.alignments.snp.strand.exclude", chrom=config["chr"]),

    shell:
        # {config[shapeit_path]} -check --input-bed {input.ug_bed} {input.ug_bim} {input.ug_fam} \
        """
        set +e
        {config[shapeit_path]} -check --input-bed {input.ug_bed} {input.ug_bim} {input.ug_fam} \
        -M {params.g_map} \
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
rule snp_flip:
    input:
        rules.snp_check.output[0],
        lambda wildcards: config["chr"][wildcards.chrom],
        ug_bed=expand(config["input_folder"] + "/{chrom}.bed",chrom=config["chr"]),
        ug_bim=expand(config["input_folder"] + "/{chrom}.bim",chrom=config["chr"]),
        ug_fam=expand(config["input_folder"] + "/{chrom}.fam",chrom=config["chr"]),
    output:
        config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/{chrom}/" + config["pop"] + "_flipped.bim",
        config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/{chrom}/" + config["pop"] + "_flipped.bed",
        config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/{chrom}/" + config["pop"] + "_flipped.fam",
        strand_rsid=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/{chrom}/" + config["pop"] + "_rsids.to_flip",
    params:
        bfiles_prefix=config["input_folder"] + "/{chrom}",
        bfiles_flipped_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/{chrom}/" + config["pop"] + "_flipped",
    shell:
        """
        set +e
        fgrep -w "Strand" {input[0]} | cut -f 4 > {output.strand_rsid}
        plink --bfile {params.bfiles_prefix} --flip {output.strand_rsid} --make-bed --out {params.bfiles_flipped_prefix}
        """

rule phase:
    input:
        # expand(config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_flipped" , ext=[".bim",".bed",".fam"]),
        lambda wildcards: config["chr"][wildcards.chrom],
        rules.snp_flip.output[0],
        rules.snp_flip.output[1],
        rules.snp_flip.output[2],
        # config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/{chrom}/" + config["pop"] + "_flipped.bim",
        # config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/{chrom}/" + config["pop"] + "_flipped.bed",
        # config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/{chrom}/" + config["pop"] + "_flipped.fam",
    params:
        # g_map=expand("/netapp/nfs/resources/1000GP_phase3/impute/genetic_map_chr{chrom}_combined_b37.txt",chrom=config["chr"]),
        g_map="/netapp/nfs/resources/1000GP_phase3/impute/genetic_map_chr{chrom}_combined_b37.txt",
        input_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/{chrom}/" + config["pop"] + "_flipped"
    output:
        # generate_shapeit_out_files("{input.chr}")
        touch(config["output_folder"]+"/"+config["pop"]+"/{chrom}.pipe.done"),
        # touch(config["output_folder"]+"/"+config["pop"]+"/"+config["chr"] +".pipe.done"),
        generate_shapeit_out_files("{chrom}")
        # chr_phased=config["output_folder"] +"/" + config["pop"] + "/" + config["ref_panel"] + "/{chrom}/chr{chrom}.haps.gz",
        # samples=config["output_folder"] +"/" + config["pop"] + "/" + config["ref_panel"] + "/{chrom}/chr{chrom}.samples",
    threads: 16
    benchmark:
        config["output_folder"]+"/"+config["pop"]+"/{chrom}.phase_rule.tsv"
    shell:
        # shapeit --input-bed gwas.bed gwas.bim gwas.fam \
        # -M genetic_map.txt \
        # -O gwas.phased
        # "{config[shapeit_path]} -B {params.input_prefix} -M {params.g_map} -O {output.chr_phased} {output.samples} -T {threads}"
        "{config[shapeit_path]} -B {params.input_prefix} -M {params.g_map} -O {output[1]} {output[2]} -T {threads}"


# rule snp_flip:
#     input:
#         rules.snp_check.output[0],
#         ug_bed=config["input_folder"] + "/" + config["chr"]+ ".bed",
#         ug_bim=config["input_folder"] + "/" + config["chr"]+ ".bim",
#         ug_fam=config["input_folder"] + "/" + config["chr"]+ ".fam"
#     output:
#         config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_flipped.bim",
#         config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_flipped.bed",
#         config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_flipped.fam",
#         strand_rsid=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_rsids.to_flip"
#     params:
#         bfiles_prefix=config["input_folder"] + "/" + config["chr"],
#         bfiles_flipped_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_flipped"
#     shell:
#         """
#         set +e
#         fgrep -w "Strand" {input[0]} | cut -f 4 > {output.strand_rsid}
#         plink --bfile {params.bfiles_prefix} --flip {output.strand_rsid} --make-bed --out {params.bfiles_flipped_prefix}
#         """

# rule phase:
#     input:
#         # expand(config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_flipped" , ext=[".bim",".bed",".fam"]),
#         rules.snp_flip.output[0],
#         rules.snp_flip.output[1],
#         rules.snp_flip.output[2]
#     params:
#         g_map="/netapp/nfs/resources/1000GP_phase3/impute/genetic_map_chr"+config["chr"]+"_combined_b37.txt",
#         input_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_flipped"
#     output:
#         # generate_shapeit_out_files("{input.chr}")
#         touch(config["output_folder"]+"/"+config["pop"]+"/"+config["chr"] +".pipe.done"),
#         generate_shapeit_out_files(config["chr"])
#     threads: 16
#     benchmark:
#         config["output_folder"]+"/"+config["pop"]+"/"+config["chr"] +".phase_rule.tsv"
#     shell:
#         # shapeit --input-bed gwas.bed gwas.bim gwas.fam \
#         # -M genetic_map.txt \
#         # -O gwas.phased
#         "{config[shapeit_path]} -B {params.input_prefix} -M {params.g_map} -O {output[1]} {output[2]} -T {threads}"

# rule pipe_finish:
#     input:
#         expand(config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_flipped" , ext=[".bim",".bed",".fam"]),
#         rules.snp_flip.output.strand_rsid
#         # expand(config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_flipped" , ext=[".bim",".bed",".fam"])
#         # expand(config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_shapeit_refpanel.alignments" , ext=[".strand",".strand.exclude"])
#     output:
#         config["output_folder"]+"/"+config["pop"]+"/"+config["chr"] +".pipe.done"
#     shell:
#         "touch {output}"

onsuccess:
    print("The workflow finished without errors!")

onerror:
    print("An error occurred in the current workflow execution!!")