# functions for the imputation pieline

# define a function to create and open a file in write mode
# got from https://stackoverflow.com/a/30582525
def createAndOpen(filename, mode):
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    return open(filename, mode)

def generate_end_of_pipeline_files(key):
    return "%s/%s/chr%s.pipe.done" % (config["output_folder"],config["pop"],key)

#function to retrieve chunk numbers by chromosome after imputation
def getChunkNumByChr(wildcards):
    from os import listdir
    import re
    from os.path import join, isfile
    # chrom=22
    # imputed_folder="/home/cocca/analyses/test_imputation_20210604/06.imputed"
    checkpoint_output = checkpoints.chunkIntervalFileGenerator.get(**wildcards).output[0]
    # return expand(output_folder+"/06.imputed/{chr}/{chr}.{g_chunk}.{ext}",ext=["vcf.gz","log"],chr=wildcards.chr,g_chunk=glob_wildcards(os.path.join(checkpoint_output, "{chr}.{g_chunk}.vcf.gz")).g_chunk)
    chr_chunk_size = len(expand(output_folder+"/06.imputed/{chr}/{chr}.{g_chunk}.vcf.gz",chr=wildcards.chr,g_chunk=glob_wildcards(os.path.join(checkpoint_output, "{chr}.{g_chunk}.int")).g_chunk))
    return chr_chunk_size

#define a conversion table to work on chrX different regions
def getChrForPhasing(wildcards):
    current_chr=wildcards.chr
    if current_chr == "23" :
        converted_chr="X"
    elif current_chr == "24":
        converted_chr="X"
    elif current_chr == "25":
        converted_chr="X"
    else :
        converted_chr=current_chr
    return converted_chr


# extract snps to flip after strand check vs reference panel
def get_flippable(infile,outfile):
    import re 
    # fgrep -w "Strand" {input[0]} | awk 'length($9)==length($10) && $5!="D" && $5!="I"' | awk '{{if($5==$6 && ($5!=$9 && $5!=$10)) print $0;else if($5!=$6){if() print $0} }}' |  | cut -f 4 | sort|uniq -u > {output.strand_rsid}
    # infile="ERBO_shapeit_refpanel.alignments.snp.strand"
    # outfile="ERBO_toflip"
    complement={'A':'T','T':'A','C':'G','G':'C'}
    strand_file=open('%s' %(infile),'r')

    # get all duplicates ids, since we need to exclude them from the flipping, if we find them only once in the flippable list
    rs_ids=[]
    for line in strand_file:
        if re.match("Strand", line.strip().split("\t")[0]):
            rs_ids.append(line.strip().split("\t")[3])
    unique_rs=list(set(rs_ids))
    duplicates=list(set([rs_id for rs_id in unique_rs if rs_ids.count(rs_id)>1]))
    #now we have the duplicates we can remove from the flippable list 
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


def get_imputation_interval(interval_file):
    import gzip
    interval=gzip.open(interval_file) if interval_file.endswith('.gz') else open(interval_file)
    interval_coord=interval.readlines()[0]
    return interval_coord

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
    new_bim=[]
    for bim_line in plink_bim_file:
        # read line
        # now check if we have a monomorphic site for wich we have to update the allele name 
        if bim_line.strip().split('\t')[4]!='0':
            # print line as it is in the output file
            _ = output_file.write(bim_line)
        else:
            # do stuff to get the correct name
            c_line=bim_line.strip().split('\t')
            c_chr=c_line[0]
            c_rsID=c_line[1]
            c_cm=c_line[2]
            c_pos=c_line[3]
            c_a1=c_line[4]
            c_a2=c_line[5]
            # since we have a cleaned rsID, we need first to get the right one among update_alleles keys
            rs_key=[x for x in all_update.keys() if re.search(re.escape(c_rsID)+'$',x)]
            #this object could be empty (it happens for sure for chrX). In this case we need to print the line as we find it
            if rs_key:
                # now we have to get the alleles from the update allele dict
                # and check which one we have. Easy case is if there is no flipping
                if c_a2 in all_update[rs_key[0]]:
                    new_a1=list(set(all_update[rs_key[0]]) - set(list(c_a2)))[0]
                else:
                    # here it could be that we have flipped data, so we need to use the complement for both a2 lookup and a1 retrieve
                    new_a1=complement[list(set(all_update[rs_key[0]]) - set(list(complement[c_a2])))[0]]
                _ = output_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' %(c_chr,c_rsID,c_cm,c_pos,new_a1,c_a2))
            else :
                _ = output_file.write(bim_line)
            # new_bim.append('%s\t%s\t%s\t%s\t%s\t%s\n' %(c_chr,c_rsID,c_cm,c_pos,new_a1,c_a2))
            # new_bim.append(bim_line)
        # if len(new_bim) % 5000 == 0:
        #     print(len(new_bim))
    output_file.close()

# define a function to collect all chunks for a chromosome
def collect_imputed_chunks(wildcards):
    from os import listdir
    import re
    from os.path import join, isfile
    # chrom=22
    # imputed_folder="/home/cocca/analyses/test_imputation_20210604/06.imputed"
    checkpoint_output = checkpoints.chunkIntervalFileGenerator.get(**wildcards).output[0]
    # return expand(output_folder+"/06.imputed/{chr}/{chr}.{g_chunk}.{ext}",ext=["vcf.gz","log"],chr=wildcards.chr,g_chunk=glob_wildcards(os.path.join(checkpoint_output, "{chr}.{g_chunk}.vcf.gz")).g_chunk)
    return expand(output_folder+"/06.imputed/{chr}/{chr}.{g_chunk}.vcf.gz",chr=wildcards.chr,g_chunk=glob_wildcards(os.path.join(checkpoint_output, "{chr}.{g_chunk}.int")).g_chunk)

# define a function to collect all chunks for a chromosome
def collect_splitted_bim(wildcards):
    from os import listdir
    import re
    from os.path import join, isfile
    # chrom=22
    # imputed_folder="/home/cocca/analyses/test_imputation_20210604/06.imputed"
    checkpoint_output = checkpoints.splitBim.get(**wildcards).output[0]
    # return expand(output_folder+"/06.imputed/{chr}/{chr}.{g_chunk}.{ext}",ext=["vcf.gz","log"],chr=wildcards.chr,g_chunk=glob_wildcards(os.path.join(checkpoint_output, "{chr}.{g_chunk}.vcf.gz")).g_chunk)
    return expand(output_folder + "/03.flipped_input/" + ref_panel + "/ReMo/"+ cohort_name+"_{chr}_allFix_flipped_{bim_chunk}_splitBIM_ReMo.bim",chr=wildcards.chr,bim_chunk=glob_wildcards(os.path.join(checkpoint_output, cohort_name+"_{chr}_allFix_flipped_{bim_chunk}_splitBIM.bim")).bim_chunk )
    # return expand(output_folder+"/06.imputed/{chr}/{chr}.{g_chunk}.vcf.gz",chr=wildcards.chr,g_chunk=glob_wildcards(os.path.join(checkpoint_output, "{chr}.{g_chunk}.int")).g_chunk)


# function to define memory requirement based on job rerun attempt, after the fist one
def get_mem_mb(def_mem, attempt):
    import math
    if attempt > 1:
        return floor(int(def_mem) + (int(def_mem) * attempt * 0.5))
    else :
        return int(def_mem)

# function to create a temporary folder
def define_tmp(config_temp):
    import subprocess
    cmd="mktemp -u -d -p %s" % (config_temp)

    exec_cmd=subprocess.Popen(cmd.split(" "),stdout=subprocess.PIPE,stderr=subprocess.PIPE) 
    stdout,stderr=exec_cmd.communicate()
    return stdout.decode().strip()

# function to flatten a list in case it is a list of nested lists
# got it from https://thispointer.com/python-convert-list-of-lists-or-nested-list-to-flat-list/
def flattenNestedList(nestedList):
    ''' Converts a nested list to a flat list '''
    flatList = []
    # Iterate over all the elements in given list
    for elem in nestedList:
        # Check if type of element is list
        if isinstance(elem, list):
            # Extend the flat list by adding contents of this element (list)
            flatList.extend(flattenNestedList(elem))
        else:
            # Append the elemengt to the list
            flatList.append(elem)    
    return flatList

# function to extract a list of duplicate sites to be removed using position
# it is possible that with the update step using the external reference, we will introduce some position duplicates
# we need a rule to take care of them, removing duplicates if they dont have an rsID
# first we want the list of duplicate rsID. We need to look for duplicate positions, then lookup the rsID in the bim file, than check if it has an rsID or a bs name from the array manufcturer.
# we will remove the bs named by using it in our list
def getDupeByPos(infile,outlist):
    import collections
    import re
    # get all duplicates positions, since we need to extract the relevant rsID, to create the exclusion list
    bim_file=open('%s' %(infile),'r')
    rs_ids=collections.defaultdict(list)
    # define a list of ids to remove
    to_rem=[]
    for line in bim_file:
        pos=line.strip().split("\t")[3]
        rsid=line.strip().split("\t")[1]
        rs_ids[pos].append(rsid)

    for k in rs_ids.keys():
        if len(rs_ids[k]) > 1:
            # we want to check if one of the ids starts with rs. If this is the case, we want to remove the other. If we have two rsId, we will remove both of them. If we don't have any rsId, we will remove both of them
            pos_to_rem=[] #we need a variable to temporary store matching ids, so we can check how many of them we have
            for v in rs_ids[k]:
                if not re.match("rs", v) :
                    to_rem.append(v)
                elif re.match("rs", v) :
                    pos_to_rem.append(v)
            if len(pos_to_rem) > 1 : # we could also check if len(pos_to_rem) == len(rs_ids[k]), but if here we have 3 duplicates by position one without rs, and 2 starting with rs, we still have to remove all of them. So we just want one item in this list
                to_rem.append(pos_to_rem)
    # get unique values                
    unique_rs=list(set(flattenNestedList(to_rem)))
    open(outlist,"w").write("\n".join(unique_rs))
    open(outlist, 'a').close()

