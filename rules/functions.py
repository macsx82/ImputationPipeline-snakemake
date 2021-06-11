# functions for the imputation pieline

# define a function to create and open a file in write mode
# got from https://stackoverflow.com/a/30582525
def createAndOpen(filename, mode):
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    return open(filename, mode)

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
                #handle monomorphic sites
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
    # legend='/home/cocca/analyses/test_imputation_20210604/04.phased_data/IGRPv1/Slo_POP_22_phased.vcf.gz'
    # modified 06/06/2021, to read vcf file
    # legend_file=gzip.open(legend) if legend.endswith('.gz') else open(legend)
    # chrom=os.path.basename(legend).split('.')[0]
    with gzip.open(legend) if legend.endswith('.gz') else open(legend) as legend_file:
    # Skip initial comments that starts with #
        while True:
            line = legend_file.readline().decode()
            # break while statement if it is not a comment line
            # i.e. does not startwith #
            if not line.startswith('#'):
                # the first line withouth the comment char is the first line I want
                first = line
                break
        # now we go to the end of the file
        legend_file.seek(os.SEEK_END)
        # go back to read the last line
        last = legend_file.readlines()[-1].decode()

    # all_pos=legend_file.read().splitlines()
    chrom=int(first.split("\t")[0])
    start=int(first.split("\t")[1])
    end=int(last.split("\t")[1])
    chunk_num=round((end-start+1)/chunk_size)
    # all_pos=legend_file.read().splitlines()    
    # start=int(all_pos[1].decode().split(" ")[1])
    # end=int(all_pos[-1].decode().split(" ")[1])
    # chunk_num=round((end-start+1)/chunk_size)
    
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

#function to recover allele names in monomorphic sites from the update allele file
#we need the update alle file for the snp array and the bim file from plink format
def update_mono_snps(allele_update,plink_bim,outfile):
    import re
    # read allele update file
    allele_update_file=open('%s' %(allele_update))
    all_update={}
    for line in allele_update_file:
        variant=line.strip().split('\t')
        rsID=variant[0]
        alleles=variant[2].split(' ')
        all_update[rsID]=alleles

    # we also need to check for flipped alleles
    complement={'A':'T','T':'A','C':'G','G':'C'}
    # read plink bim
    plink_bim_file=open('%s' %(plink_bim),'r')
    # open the stream to the output file
    output_file=open(outfile,'w')

    for bim_line in plink_bim_file:
        # read line
        # now check if we have a monomorphic site for wich we have to update the allele name 
        if bim_line.strip().split('\t')[4]=='0':
            # do stuff to get the correct name
            c_line=bim_line.strip().split('\t')
            c_chr=c_line[0]
            c_rsID=c_line[1]
            c_cm=c_line[2]
            c_pos=c_line[3]
            c_a1=c_line[4]
            c_a2=c_line[5]
            # since we have a cleaned rsID, we need firs to get the right one among update_alleles keys
            rs_key=[x for x in all_update.keys() if re.search(c_rsID+'$',x)]
            # now we have to get the alleles from the update allele dict
            # and check which one we have. Easy case is if there is no flipping
            if c_a2 in all_update[rs_key[0]]:
                new_a1=list(set(all_update[rs_key[0]]) - set(list(c_a2)))[0]
            else:
                # here it could be that we have flipped data, so we need to use the complement for both a2 lookup and a1 retrieve
                new_a1=complement[list(set(all_update[rs_key[0]]) - set(list(complement[c_a2])))[0]]
            _ = output_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' %(c_chr,c_rsID,c_cm,c_pos,new_a1,c_a2))
        else:
            # print line as it is in the output file
            _ = output_file.write(bim_line)
    output_file.close()

# define a function to collect all chunks for a chromosome
# def collect_imputed_chunks(imputed_folder,chrom):
#     from os import listdir
#     import re
#     from os.path import join, isfile
#     # chrom=22
#     # imputed_folder="/home/cocca/analyses/test_imputation_20210604/06.imputed"
#     imputed_pattern=str(chrom)+".\d+.vcf.gz"
#     # read the folder content and return it as a list
#     vcf_files= [join(imputed_folder,str(chrom),vcf_f) for vcf_f in listdir(join(imputed_folder,str(chrom))) if re.search(imputed_pattern,vcf_f)]

#     return vcf_files

def collect_imputed_chunks(wildcards):
    from os import listdir
    import re
    from os.path import join, isfile
    # chrom=22
    # imputed_folder="/home/cocca/analyses/test_imputation_20210604/06.imputed"
    checkpoint_output = checkpoints.chunkIntervalFileGenerator.get(**wildcards).output[0]
    return expand(output_folder+"/06.imputed/{chr}/{chr}.{g_chunk}.{ext}",ext=["vcf.gz","log"],chr=wildcards.chr,g_chunk=glob_wildcards(os.path.join(checkpoint_output, "{chr}.{g_chunk}.vcf.gz")).g_chunk)
