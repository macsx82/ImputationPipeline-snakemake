#Snakefile for the Integrated Prioritization Workflow pipeline
#
#Author: Massimiliano Cocca
#
#
# configfile: "config.yaml"

print (config['chr_to_phase'].keys())


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
        config["output_folder"]+"/"+config["pop"]+"/"+config["chr"] +".pipe.done"

#We assume all our data has been already strand oriented with plink
#We will orient the data to match the reference panel orientation, using shapeit
# Input files will be the plink genotypes. We will get them from a config file
rule snp_check:
    input:
        ug_bed=config["input_folder"] + "/" + config["chr"]+ ".bed",
        ug_bim=config["input_folder"] + "/" + config["chr"]+ ".bim",
        ug_fam=config["input_folder"] + "/" + config["chr"]+ ".fam",
        rp_hap=config["ref_panel_base_folder"]+ "/" + config["ref_panel"]+ "/" + config["chr"]+ "/" + config["chr"] + "." + config["ref_panel"] + ".hap.gz",
        rp_legend=config["ref_panel_base_folder"]+ "/" + config["ref_panel"]+ "/" + config["chr"]+ "/" + config["chr"] + "." + config["ref_panel"] + ".legend.gz",
        rp_samples=config["ref_panel_base_folder"]+ "/" + config["ref_panel"]+ "/" + config["chr"]+ "/" + config["chr"] + "." + config["ref_panel"] + ".samples"
    params:
        g_map="/netapp/nfs/resources/1000GP_phase3/impute/genetic_map_chr"+config["chr"]+"_combined_b37.txt",
        output_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_shapeit_refpanel.alignments"
    output:
        # expand(config["output_folder"] + "/" + config["pop"] + "/" + config["chr"]+ "/chr"+config["chr"]+"_relate_pos_sel{ext}", ext=[".freq",".lin",".sele"])
        # expand(config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_shapeit_refpanel.alignments.snp" , ext=[".strand",".strand.exclude"])
        config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_shapeit_refpanel.alignments.snp.strand",
        config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_shapeit_refpanel.alignments.snp.strand.exclude"
        # /home/cocca/analyses/imputation/03032020/CARL/IGRPv1/18/CARL_shapeit_refpanel.alignments.snp
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
        ug_bed=config["input_folder"] + "/" + config["chr"]+ ".bed",
        ug_bim=config["input_folder"] + "/" + config["chr"]+ ".bim",
        ug_fam=config["input_folder"] + "/" + config["chr"]+ ".fam"
    output:
        config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_flipped.bim",
        config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_flipped.bed",
        config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_flipped.fam",
        strand_rsid=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_rsids.to_flip"
    params:
        bfiles_prefix=config["input_folder"] + "/" + config["chr"],
        bfiles_flipped_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_flipped"
    shell:
        """
        set +e
        fgrep -w "Strand" {input[0]} | cut -f 4 > {output.strand_rsid}
        plink --bfile {params.bfiles_prefix} --flip {output.strand_rsid} --make-bed --out {params.bfiles_flipped_prefix}
        """

rule phase:
    input:
        # expand(config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_flipped" , ext=[".bim",".bed",".fam"]),
        rules.snp_flip.output[0],
        rules.snp_flip.output[1],
        rules.snp_flip.output[2]
    params:
        g_map="/netapp/nfs/resources/1000GP_phase3/impute/genetic_map_chr"+config["chr"]+"_combined_b37.txt",
        input_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_flipped"
    output:
        # generate_shapeit_out_files("{input.chr}")
        touch(config["output_folder"]+"/"+config["pop"]+"/"+config["chr"] +".pipe.done"),
        generate_shapeit_out_files(config["chr"])
    threads: 8
    benchmark:
        config["output_folder"]+"/"+config["pop"]+"/"+config["chr"] +".phase_rule.tsv"
    shell:
        # shapeit --input-bed gwas.bed gwas.bim gwas.fam \
        # -M genetic_map.txt \
        # -O gwas.phased
        "{config[shapeit_path]} -B {params.input_prefix} -M {params.g_map} -O {output[0]} {output[1]} -T {threads}"

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