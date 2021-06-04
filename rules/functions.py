# functions for the imputation pieline

def generate_end_of_pipeline_files(key):
    return "%s/%s/chr%s.pipe.done" % (config["output_folder"],config["pop"],key)

# extract snps to flip after strand check vs reference panel
def get_flippable(infile,outfile):
    import re 
    # fgrep -w "Strand" {input[0]} | awk 'length($9)==length($10) && $5!="D" && $5!="I"' | awk '{{if($5==$6 && ($5!=$9 && $5!=$10)) print $0;else if($5!=$6){if() print $0} }}' |  | cut -f 4 | sort|uniq -u > {output.strand_rsid}
    # infile="ERBO_shapeit_refpanel.alignments.snp.strand"
    # outfile="ERBO_toflip"
    complement={'A':'T','T':'A','C':'G','G':'C'}
    strand_file=open('%s' %(infile),'r')

    # get all duplicates ids, since we need to exclude them form the flipping, if we find them only once in the flippable list
    rs_ids=[]
    for line in strand_file:
        if re.match("Strand", line.strip().split("\t")[0]):
            rs_ids.append(line.strip().split("\t")[3])
    unique_rs=list(set(rs_ids))
    duplicates=list(set([rs_id for rs_id in unique_rs if rs_ids.count(rs_id)>1]))
    #now we have the duplciates we can remove from the flippable list 
    strand_file.seek(0)
    flippable=[]
    for line in strand_file:
        if re.match("Strand", line.strip().split("\t")[0]):
            c_to_flip=line.strip().split("\t")
            # check that we are working with snps and not indels
            if len(c_to_flip[8]) == len(c_to_flip[9]) and c_to_flip[4] != "D" and c_to_flip[4] != "I" :
                if c_to_flip[4] == c_to_flip[5] and (c_to_flip[4] != c_to_flip[8] and c_to_flip[4] != c_to_flip[9]):
                    flippable.append(c_to_flip[3])
                    # print(c_to_flip[3], file=open(outfile,"a"))
                elif (c_to_flip[4] != c_to_flip[5]):
                    #need to check if there are multiallelic sites
                    if (complement.get(c_to_flip[4]) == c_to_flip[8] or complement.get(c_to_flip[4]) == c_to_flip[9]) and (complement.get(c_to_flip[5]) == c_to_flip[8] or complement.get(c_to_flip[5]) == c_to_flip[9]):
                        flippable.append(c_to_flip[3])
                # else:
    # remove duplicates from flippable
    unique_flippable=[rs_id for rs_id in flippable if rs_id not in duplicates]
    open(outfile,"w").write("\n".join(unique_flippable))
    # print(unique_flippable, file=open(outfile,"w"))
    open(outfile, 'a').close()

def get_chunk_num(legend,chunk_size):
    import gzip
    import os
    # get start and end of the chromosome
    legend_file=gzip.open(legend) if legend.endswith('.gz') else open(legend)
    all_pos=legend_file.read().splitlines()
    chrom=os.path.basename(legend).split('.')[0]
    start=int(all_pos[1].decode().split(" ")[1])
    end=int(all_pos[-1].decode().split(" ")[1])
    chunk_num=round((end-start+1)/chunk_size)
    return chrom,start,end,chunk_num

# create a dictionary with chroms and the relative chunk number extracted from the selected legend file for the imputation
def get_chunk_by_chr(chr,legend,chunk_size):
    import gzip
    legend_file=gzip.open(legend) if legend.endswith('.gz') else open(legend)
    all_pos=legend_file.read().splitlines()
    start=int(all_pos[1].decode().split(" ")[1])
    end=int(all_pos[-1].decode().split(" ")[1])
    chunk_num=round((end-start+1)/chunk_size)
    return chunk_num

def get_imputation_interval(interval_file):
    import gzip
    interval=gzip.open(interval_file) if interval_file.endswith('.gz') else open(interval_file)
    interval_coord=interval.readlines()[0]
    return interval_coord

def create_chunks(legend,chunk_size,chunk):
    # import gzip
    # get start and end of the chromosome
    chrom,start,end,chunk_num=get_chunk_num(legend,chunk_size)

    # legend_file=gzip.open(legend) if legend.endswith('.gz') else open(legend)
    # all_pos=legend_file.read().splitlines()
    # start=int(all_pos[1].decode().split(" ")[1])
    # end=int(all_pos[-1].decode().split(" ")[1])
    # chunk_num=round((end-start+1)/chunk_size)

    # for chunk in list(range(1,chunk_num+1)):
        # "{:02d}".format(chunk)
    if chunk > 1 :
        chunk_begin=start+(chunk-1)*chunk_size+1
    else:
        chunk_begin=start
    if chunk == chunk_num:
        chunk_end=end
    else:
        chunk_end=start+(chunk*chunk_size)
    interval="%s %s" % (chunk_begin, chunk_end)
    return interval

# function to get rsId of sites for which we couldn't assign alleles based on a refernce file
def get_not_assigned_snps(outfile,*infile):
    import re
    rs_ids=[]
    for inputfile in infile:
        c_input_file=open('%s' %(inputfile),'r')
        for line in c_input_file:
            # if (re.match("Warning: Impossible A1 allele assignment",line.strip())):
            if (re.match("Warning: Impossible A1 allele assignment",line.strip()) or re.match("Warning: Impossible A2 allele assignment",line.strip())):
                print(line.strip())
                rs_ids.append((line.strip().split(" ")[7]).replace('.',''))
    unique_rs=list(set(rs_ids))
    open(outfile,"w").write("\n".join(unique_rs))
    open(outfile, 'a').close()
