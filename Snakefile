#Snakefile for the Internal Imputation pipeline
#
#Author: Massimiliano Cocca
#
#
# configfile: "config.yaml"

def get_plink_input_files(key):
    ug_bed="%s/%s/%s.bed" % (config["input_folder"],key,key)
    ug_bim="%s/%s/%s.bim" % (config["input_folder"],key,key)
    ug_fam="%s/%s/%s.fam" % (config["input_folder"],key,key)
    return ug_bed,ug_bim,ug_fam

def get_shapeit_input_files(key):
    rp_hap="%s/%s/%s/%s.%s.hap.gz" % (config["ref_panel_base_folder"],config["ref_panel"],key ,key ,config["ref_panel"]),
    rp_legend="%s/%s/%s/%s.%s.legend.gz" % (config["ref_panel_base_folder"],config["ref_panel"],key ,key ,config["ref_panel"]),
    rp_samples="%s/%s/%s/%s.%s.sample" % (config["ref_panel_base_folder"],config["ref_panel"],key ,key ,config["ref_panel"])
    return rp_hap,rp_legend,rp_samples

def generate_shapeit_out_files(key):
    # config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/"
    chr_phased= "%s/%s/%s/%s/chr%s.haps.gz" % (config["output_folder"],config["pop"],config["ref_panel"],key,key)
    samples= "%s/%s/%s/%s/chr%s.sample" % (config["output_folder"],config["pop"],config["ref_panel"],key,key)
    return chr_phased,samples


def generate_end_of_pipeline_files(key):
    return "%s/%s/chr%s.pipe.done" % (config["output_folder"],config["pop"],key)


def get_flippable(infile,outfile):
    import re 
    # fgrep -w "Strand" {input[0]} | awk 'length($9)==length($10) && $5!="D" && $5!="I"' | awk '{{if($5==$6 && ($5!=$9 && $5!=$10)) print $0;else if($5!=$6){if() print $0} }}' |  | cut -f 4 | sort|uniq -u > {output.strand_rsid}
    # infile="ERBO_shapeit_refpanel.alignments.snp.strand"
    # outfile="ERBO_toflip"
    complement={'A':'T','T':'A','C':'G','G':'C'}
    strand_file=open('%s' %(infile),'r')
    for line in strand_file:
        if re.match("Strand", line.strip().split("\t")[0]):
            c_to_flip=line.strip().split("\t")
            if len(c_to_flip[8]) == len(c_to_flip[9]) and c_to_flip[4] != "D" and c_to_flip[4] != "I" :
                if c_to_flip[4] == c_to_flip[5] and (c_to_flip[4] != c_to_flip[8] and c_to_flip[4] != c_to_flip[9]):
                    print(c_to_flip[3], file=open(outfile,"a"))
                elif (c_to_flip[4] != c_to_flip[5]):
                    #need to check if there are multiallelic sites
                    if (complement.get(c_to_flip[4]) == c_to_flip[8] or complement.get(c_to_flip[4]) == c_to_flip[9]) and (complement.get(c_to_flip[5]) == c_to_flip[8] or complement.get(c_to_flip[5]) == c_to_flip[9]):
                        print(c_to_flip[3], file=open(outfile,"a"))


    

input_f=config["input_folder"]

#define parameter useful to cluster job submission
localrules: all

#define rules
rule all:
    input:
        # config["output_folder"]+"/"+config["pop"]+"/" + config["pop"] + ".pipe.done"
        # lambda wildcards: config["chr"][wildcards.chrom],
        # expand(config["output_folder"]+"/"+config["pop"]+"/{chrom}.pipe.done", chrom=config["chr"])
        config["output_folder"]+"/"+config["pop"]+"/" + config["chr"] + ".pipe.done"
        

#We assume all our data has been already strand oriented with plink
#We will orient the data to match the reference panel orientation, using shapeit
# Input files will be the plink genotypes. We will get them from a config file
###############################################################################################
rule snp_check:
    input:
        ug_bed=config["input_folder"] + "/" + config["chr"] + "/" + config["chr"]+ ".bed",
        ug_bim=config["input_folder"] + "/" + config["chr"] + "/" + config["chr"]+ ".bim",
        ug_fam=config["input_folder"] + "/" + config["chr"] + "/" + config["chr"]+ ".fam",
        rp_hap=config["ref_panel_base_folder"]+ "/" + config["ref_panel"]+ "/" + config["chr"]+ "/" + config["chr"] + "." + config["ref_panel"] + ".hap.gz",
        rp_legend=config["ref_panel_base_folder"]+ "/" + config["ref_panel"]+ "/" + config["chr"]+ "/" + config["chr"] + "." + config["ref_panel"] + ".legend.gz",
        rp_samples=config["ref_panel_base_folder"]+ "/" + config["ref_panel"]+ "/" + config["chr"]+ "/" + config["chr"] + "." + config["ref_panel"] + ".samples"
    params:
        g_map="/netapp/nfs/resources/1000GP_phase3/impute/genetic_map_chr"+config["chr"]+"_combined_b37.txt",
        output_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_shapeit_refpanel.alignments"
    output:
        config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_shapeit_refpanel.alignments.snp.strand",
        config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_shapeit_refpanel.alignments.snp.strand.exclude"
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

rule snp_flip_file:
    input:
        rules.snp_check.output[0],
    output:
        strand_rsid=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_rsids.to_flip"
    params:
        bfiles_prefix=config["input_folder"] + "/" + config["chr"] + "/" + config["chr"],
        bfiles_flipped_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_" + config["chr"] +  "_flipped"
    run:
        get_flippable(input[0],output.strand_rsid)
    
rule snp_flip:
    input:
        rules.snp_flip_file.output[0],
        ug_bed=config["input_folder"] + "/" + config["chr"] + "/" + config["chr"]+ ".bed",
        ug_bim=config["input_folder"] + "/" + config["chr"] + "/" + config["chr"]+ ".bim",
        ug_fam=config["input_folder"] + "/" + config["chr"] + "/" + config["chr"]+ ".fam"
    output:
        config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_" + config["chr"]+ "_flipped.bim",
        config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_" + config["chr"]+ "_flipped.bed",
        config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_" + config["chr"]+ "_flipped.fam",
        # strand_rsid=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_rsids.to_flip"
    params:
        bfiles_prefix=config["input_folder"] + "/" + config["chr"] + "/" + config["chr"],
        bfiles_flipped_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_" + config["chr"] +  "_flipped"
    shell:
        """
        set +e
        
        plink --bfile {params.bfiles_prefix} --flip {input[0]} --make-bed --out {params.bfiles_flipped_prefix}
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

rule phase:
    input:
        # expand(config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_flipped" , ext=[".bim",".bed",".fam"]),
        rules.snp_flip.output[0],
        rules.snp_flip.output[1],
        rules.snp_flip.output[2]
    params:
        g_map="/netapp/nfs/resources/1000GP_phase3/impute/genetic_map_chr"+config["chr"]+"_combined_b37.txt",
        input_prefix=config["output_folder"] + "/" + config["pop"] + "/" + config["ref_panel"] + "/" +config["chr"] + "/" + config["pop"] + "_" + config["chr"] + "_flipped"
    output:
        # generate_shapeit_out_files("{input.chr}")
        touch(config["output_folder"]+"/"+config["pop"]+"/"+config["chr"] +".pipe.done"),
        generate_shapeit_out_files(config["chr"])
    threads: 16
    benchmark:
        config["output_folder"]+"/"+config["pop"]+"/"+config["chr"] +".phase_rule.tsv"
    shell:
        # shapeit --input-bed gwas.bed gwas.bim gwas.fam \
        # -M genetic_map.txt \
        # -O gwas.phased
        "{config[shapeit_path]} -B {params.input_prefix} -M {params.g_map} -O {output[1]} {output[2]} -T {threads}"

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