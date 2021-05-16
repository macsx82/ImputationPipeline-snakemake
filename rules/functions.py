# functions for the imputation pieline

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

# function to split plink input files by chromosome
# def splitPlinkFiles(prefix,cohort_name,chr):
    

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
                else:
                    open(outfile, 'a').close()
