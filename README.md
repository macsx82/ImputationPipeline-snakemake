
#03/03/2020

Inizialization of the snakemake imputation pipeline!

* First rule will be the snp check using shapeit to align genotypes to the reference panel
* Second rule will be the update of misaligned data
* third rule will be the Phasing
* Finally we will submit the imputation

test command for rules 1-3:
```bash
snakemake --configfile config_chr_18.yaml --cores 8 
snakemake --configfile config_chr_21_22.yaml --cores 8 
```

We need to submit to the cluster:

```bash
# snakemake --configfile config_chr_18.yaml --cores 8 --cluster-config SGE_cluster.json --cluster "qsub -N ${pop}_${chr}_STEP6 -m ea -M {cluster.user_mail} -pe {cluster.parall_env} {threads} -o {config.output_folder}/\$JOB_ID_{config.pop}_{config.chr}.log -e {config.output_folder}/\$JOB_ID_STEP6_{config.pop}_{config.chr}.e -V -l h_vmem={cluster.mem} -q {cluster.queue}"
# snakemake --jobs 50 --configfile config_chr_18.yaml --cluster-config ~/scripts/pipelines/ImputationPipeline-snakemake/SGE_cluster.json --cluster "qsub -N ${pop}_${chr}_STEP6 -m ea -M {cluster.user_mail} -pe {cluster.parall_env} {threads} -o {config.output_folder}/\$JOB_ID_{config.pop}_{config.chr}.log -e {config.output_folder}/\$JOB_ID_STEP6_{config.pop}_{config.chr}.e -V -l h_vmem={cluster.mem} -q {cluster.queue}"
snakemake --jobs 50 --configfile config_chr_18.yaml --cluster-config ~/scripts/pipelines/ImputationPipeline-snakemake/SGE_cluster.json --cluster "qsub -N {config[pop]}_{config[chr]}_{rule} -m ea -M {cluster.user_mail} -pe {cluster.parall_env} {threads} -o {config[output_folder]}/\$JOB_ID_{config[pop]}_{config[chr]}.log -e {config[output_folder]}/\$JOB_ID_{config[pop]}_{config[chr]}.e -V -l h_vmem={cluster.mem} -q {cluster.queue}"
snakemake --jobs 50 --configfile config_chr_21_22.yaml --cluster-config ~/scripts/pipelines/ImputationPipeline-snakemake/SGE_cluster.json --cluster "qsub -N {config[pop]}_{wildcards.chrom}_{rule} -m ea -M {cluster.user_mail} -pe {cluster.parall_env} {threads} -o {config[output_folder]}/\$JOB_ID_{config[pop]}_{wildcards.chrom}.log -e {config[output_folder]}/\$JOB_ID_{config[pop]}_{wildcards.chrom}.e -V -l h_vmem={cluster.mem} -q {cluster.queue}"
```

qsub -N ${pop}_${chr}_STEP6 -m ea -M massimiliano.cocca@burlo.trieste.it -o ${outfolder}/\$JOB_ID_STEP6_${pop}_${chr}.log -e ${outfolder}/\$JOB_ID_STEP6_${pop}_${chr}.log.e -V -l h_vmem=${m} -q fast

snakemake -j 100 --cluster-config cluster.json --cluster "sbatch -A {cluster.account} --mem={cluster.mem} -t {cluster.time} -c {threads}"


Sample command:

```bash
snakemake -s ~/scripts/pipelines/ImputationPipeline-snakemake/Snakefile -p -r --configfile /home/cocca/analyses/test_imputation_20210604/config_test_2.yaml --cores 10 --keep-going
```

---
#05/03/2020

In order to fix the misaligned data before the imputation and phasing step, we check the genotypes against the reference panel selected.
We will flip all data with the "Strand" warning generated by shapeit, but we will remove data misaligned to INDELs and Multiallelic sites (removing ALL duplicated occurences)

```bash
fgrep -w "Strand" {input[0]} | awk 'length($9)==length($10) && $5!="D" && $5!="I"' | cut -f 4 | sort|uniq -u > {output.strand_rsid}
```

---
#05/05/2020
Create a wrapper to generate the config file

```bash
# snakemake -s ~/scripts/pipelines/ImputationPipeline-snakemake/Snakefile -p -r --jobs 50 --configfile /home/cocca/analyses/test_imputation/config_test_2.yaml --cluster-config ~/scripts/pipelines/ImputationPipeline-snakemake/SGE_cluster.json --cluster "qsub -N {config[cohort_name]}_{wildcards.chr}_{rule} -V -cwd -m ea -M {cluster.user_mail} -pe {cluster.parall_env} {threads} -o {config[paths][log_folder]}/\$JOB_ID_{config[cohort_name]}_{wildcards.chr}.log -e {config[paths][log_folder]}/\$JOB_ID_{config[cohort_name]}_{wildcards.chr}.e -V -l h_vmem={cluster.mem} -q {cluster.queue}"
snakemake -s ~/scripts/pipelines/ImputationPipeline-snakemake/Snakefile -p -r --jobs 50 --configfile /home/cocca/analyses/test_imputation_20210604/config_test_2.yaml --cluster-config ~/scripts/pipelines/ImputationPipeline-snakemake/SGE_cluster.json --keep-going --cluster "qsub -N {config[cohort_name]}_{rule} -V -cwd -m ea -M {cluster.user_mail} -pe {cluster.parall_env} {threads} -o {config[paths][log_folder]}/\$JOB_ID_{config[cohort_name]}_{rule}.log -e {config[paths][log_folder]}/\$JOB_ID_{config[cohort_name]}_{rule}.e -V -l h_vmem={cluster.mem} -q {cluster.queue}"
# echo "snakemake -s ~/scripts/pipelines/ImputationPipeline-snakemake/Snakefile -p -r --configfile /home/cocca/analyses/test_imputation_20210604/config_test_2.yaml --cores 30 --keep-going"
# echo "snakemake -s ~/scripts/pipelines/ImputationPipeline-snakemake/Snakefile -p -r --configfile /home/cocca/analyses/test_imputation/config_test_2.yaml --cores 32 " | qsub -N test_imputation -o ./
```

---
#3/6/2021

We needed to create a resource for allele check and alignment.
We need to use opnl
We used version 154 of dbSNP.

```bash
for chr in {1..22} X Y MT
do
echo "bcftools view -m2 -M2 /netapp/nfs/resources/dbSNP/human_9606_b154_GRCh37p13/VCF/GCF_000001405.25.chr${chr}.vcf.gz |bcftools query -f'%CHROM\t%POS\t%ID\t%REF\t%ALT\n' > /netapp/nfs/resources/dbSNP/human_9606_b154_GRCh37p13/VCF/GCF_000001405.25.chr${chr}.dbSNP154.tab" | qsub -N get_${chr}_data -V -cwd -l h_vmem=15G -q fast
done
```

Merge in a single table

```bash
(for chr in {1..22} X Y MT
do
cat /netapp/nfs/resources/dbSNP/human_9606_b154_GRCh37p13/VCF/GCF_000001405.25.chr${chr}.dbSNP154.tab
done) > /netapp/nfs/resources/dbSNP/human_9606_b154_GRCh37p13/GCF_000001405.25.dbSNP154.tab
```

---
#4/6/2021

The first step of a2 update, gives 47 warning for "Impossible A2 allele assignment".
Working in folder ~/analyses/test_imputation/00.splitted_input.

```bash
fgrep "Impossible A2 allele assignment" *.log |cut -f 8 -d " "|sed 's/\.//g'
```

We want to see what kind of variants are those, in dbSNP data and in our data.

```bash
fgrep -w -f <(fgrep "Impossible A2 allele assignment" *.log |cut -f 8 -d " "|sed 's/\.//g') /netapp/nfs/resources/dbSNP/human_9606_b154_GRCh37p13/GCF_000001405.25.dbSNP154.tab > ../impossible_a2_assign_dbSNP.tab
fgrep -w -f <(fgrep "Impossible A2 allele assignment" *.log |cut -f 8 -d " "|sed 's/\.//g') Slo_POP_snps_only_mapUpdateExt.bim
```

We can see that they are all flipped snps in the snp array dataset. So, actually, it is ok.

We also have 458 warnings for "Impossible A1 allele assignment"

```bash
fgrep "Impossible A1 allele assignment" *.log |cut -f 8 -d " "|sed 's/\.//g'
```

We want to see what kind of variants are those, in dbSNP data and in our data.

```bash
fgrep -w -f <(fgrep "Impossible A1 allele assignment" *.log |cut -f 8 -d " "|sed 's/\.//g') /netapp/nfs/resources/dbSNP/human_9606_b154_GRCh37p13/GCF_000001405.25.dbSNP154.tab > ../impossible_a1_assign_dbSNP.tab
fgrep -w -f <(fgrep "Impossible A1 allele assignment" *.log |cut -f 8 -d " "|sed 's/\.//g') *.bim
```


In the menawhile, we need to rebuild the ref table with dbSNP, removing INDELs

```bash
for chr in {1..22} X Y MT
do
echo "bcftools view -m2 -M2 -v snps /netapp/nfs/resources/dbSNP/human_9606_b154_GRCh37p13/VCF/GCF_000001405.25.chr${chr}.vcf.gz |bcftools query -f'%CHROM\t%POS\t%ID\t%REF\t%ALT\n' > /netapp/nfs/resources/dbSNP/human_9606_b154_GRCh37p13/VCF/GCF_000001405.25.chr${chr}.dbSNP154.tab" | qsub -N get_${chr}_data -V -cwd -l h_vmem=15G -q fast
done
```

and merge all back in a single table

```bash
(for chr in {1..22} X
do
cat /netapp/nfs/resources/dbSNP/human_9606_b154_GRCh37p13/VCF/GCF_000001405.25.chr${chr}.dbSNP154.tab
done) > /netapp/nfs/resources/dbSNP/human_9606_b154_GRCh37p13/GCF_000001405.25.SNPS.dbSNP154.tab
```

One problem with the allele update seems to be that we first updated a2, than a1. We tried also reversing the update, a1 first, than a2. This will still include errors due to the presence of INDELs in the reference file from dbSNP, but we'll fix that later.

We have 40 warnings for "Impossible A1 allele assignment"

```bash
fgrep "Impossible A1 allele assignment" *.log |cut -f 8 -d " "|sed 's/\.//g'
```

We want to see what kind of variants are those, in dbSNP data and in our data.

```bash
fgrep -w -f <(fgrep "Impossible A1 allele assignment" *.log |cut -f 8 -d " "|sed 's/\.//g') /netapp/nfs/resources/dbSNP/human_9606_b154_GRCh37p13/GCF_000001405.25.SNPS.dbSNP154.tab > ../impossible_a1_assign_dbSNP.tab
fgrep -w -f <(fgrep "Impossible A1 allele assignment" *.log |cut -f 8 -d " "|sed 's/\.//g') *.bim
```

We have 409 warnings for "Impossible A2 allele assignment"

```bash
fgrep "Impossible A2 allele assignment" *.log |cut -f 8 -d " "|sed 's/\.//g'
```

We want to see what kind of variants are those, in dbSNP data and in our data.

```bash
fgrep -w -f <(fgrep "Impossible A2 allele assignment" *.log |cut -f 8 -d " "|sed 's/\.//g') /netapp/nfs/resources/dbSNP/human_9606_b154_GRCh37p13/GCF_000001405.25.dbSNP154.tab > ../impossible_a2_assign_dbSNP.tab
fgrep -w -f <(fgrep "Impossible A2 allele assignment" *.log |cut -f 8 -d " "|sed 's/\.//g') Slo_POP_snps_only_mapUpdateExt.bim
```

So, the most likely reason for the impossible allele assignment is a flipping problem. 
We will than use the first rule to get all the unassigned snps, flip them in the original dataset and rerun the rule.
This way, we could still have some unassigned variant, but those should be only due to actual mismatch with the reference data.

rs375065198
rs199994882


#05/06/2020
Sample command to run the pipeline on system queues

```bash
snakemake -s ~/scripts/pipelines/ImputationPipeline-snakemake/Snakefile -p -r --jobs 50 --configfile /home/cocca/analyses/test_imputation_20210604/config_test_2.yaml --cluster-config ~/scripts/pipelines/ImputationPipeline-snakemake/SGE_cluster.json --keep-going --cluster "qsub -N {config[cohort_name]}_{rule} -V -cwd -m ea -M {cluster.user_mail} -pe {cluster.parall_env} {threads} -o {config[paths][log_folder]}/\\$JOB_ID_{config[cohort_name]}_{rule}.log -e {config[paths][log_folder]}/\\$JOB_ID_{config[cohort_name]}_{rule}.e -V -l h_vmem={cluster.mem} -q {cluster.queue}"
```

---
#6/6/2021

Added rules to retrieve monomorphic sites.

We need also to update the imputation rule to use IMPUTE5.

We need to generate a new set of reference panels, using custom IMPUTE5 format.
First, convert to VCF format:

```bash
basepath=/shared/resources/references/IGRPv1
for chr in {1..22}
do

mkdir -p ${basepath}/VCF/${chr}

echo "bcftools convert --haplegendsample2vcf ${basepath}/HAP_LEGEND_SAMPLES/${chr}/${chr}.IGRPv1.hap.gz,${basepath}/HAP_LEGEND_SAMPLES/${chr}/${chr}.IGRPv1.legend.gz,${basepath}/HAP_LEGEND_SAMPLES/${chr}/${chr}.IGRPv1.samples -O z -o ${basepath}/VCF/${chr}/${chr}.IGRPv1.vcf.gz;tabix -p vcf ${basepath}/VCF/${chr}/${chr}.IGRPv1.vcf.gz" |qsub -N convert_${chr}_vcf -V -cwd -l h_vmem=15G -q fast
done
```

We need to clean the VCF files from variants other than SNPs and Indels, since it seems there are issues with characters other than ATCG in ALT field (<CN0>, <INS:...:...>, etc)

```bash
basepath=/shared/resources/references/IGRPv1
for chr in {1..21}
do

mkdir -p ${basepath}/VCF/${chr}

echo "bcftools view -e\"ALT~'<'\" ${basepath}/VCF/${chr}/${chr}.IGRPv1.vcf.gz -O z -o ${basepath}/VCF/${chr}/${chr}.IGRPv1.clean.vcf.gz;tabix -p vcf ${basepath}/VCF/${chr}/${chr}.IGRPv1.clean.vcf.gz" |qsub -N clean_${chr}_vcf -V -cwd -l h_vmem=15G -q fast
done
```


After this first conversion and the cleaning, we can convert all to IMP5 format using the converter by Marchini et al.

```bash
basepath=/shared/resources/references/IGRPv1
# for chr in {7..21}
# for chr in {1..22}
for chr in 22
do
mkdir -p ${basepath}/IMP5/${chr}
# echo "/shared/software/impute5_v1.1.4/imp5Converter_1.1.4_static --h ${basepath}/VCF/${chr}/${chr}.IGRPv1.vcf.gz --r ${chr} --o ${basepath}/IMP5/${chr}/${chr}.IGRPv1.imp5" |qsub -N convert_${chr}_imp5 -V -cwd -l h_vmem=15G -q fast
# echo "/shared/software/impute5_v1.1.5/imp5Converter_1.1.5_static --h ${basepath}/VCF/${chr}/${chr}.IGRPv1.vcf.gz --r ${chr} --o ${basepath}/IMP5/${chr}/${chr}.IGRPv1.imp5" |qsub -N convert_${chr}_imp5 -V -cwd -l h_vmem=15G -q fast
# echo "/shared/software/impute5_v1.1.5/imp5Converter_1.1.5_static --h ${basepath}/VCF/${chr}/${chr}.IGRPv1.clean.vcf.gz --r ${chr} --o ${basepath}/IMP5/${chr}/${chr}.IGRPv1.imp5" |qsub -N convert_${chr}_imp5 -V -cwd -l h_vmem=15G -q fast
echo "/shared/software/impute5_v1.1.5/imp5Chunker_1.1.5_static --h ${basepath}/VCF/${chr}/${chr}.IGRPv1.clean.vcf.gz --r ${chr} --g /home/cocca/analyses/test_imputation_20210604/04.phased_data/IGRPv1/Slo_POP_22_phased.vcf.gz --o /home/cocca/analyses/test_imputation_20210604/04.phased_data/IGRPv1/Slo_POP_22_phased.coordinates.txt" |qsub -N chunk_${chr}_imp5 -V -cwd -l h_vmem=15G -q fast
done
```


Chunk ID / chromosome ID / Buffered region / Imputation region /Length / Number of target markers / Number of reference markers
0       22      22:16050036-26192765    22:16050036-25941503    9883102 2032    318194
1       22      22:25688138-35225787    22:25941504-34849841    8905956 2031    288410
2       22      22:34550905-44869030    22:34849842-44596965    9744107 2031    333267
3       22      22:44340904-51244237    22:44596966-51244237    6616406 2030    273868
