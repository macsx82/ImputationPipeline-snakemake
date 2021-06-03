
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
snakemake -s ~/scripts/pipelines/ImputationPipeline-snakemake/Snakefile -p -r --jobs 50 --configfile /home/cocca/analyses/test_imputation/config_test_2.yaml --cluster-config ~/scripts/pipelines/ImputationPipeline-snakemake/SGE_cluster.json --cluster "qsub -N {config[cohort_name]}_{wildcards.chr}_{rule} -V -cwd -m ea -M {cluster.user_mail} -pe {cluster.parall_env} {threads} -o {config[paths][log_folder]}/\$JOB_ID_{config[cohort_name]}_{wildcards.chr}.log -e {config[paths][log_folder]}/\$JOB_ID_{config[cohort_name]}_{wildcards.chr}.e -V -l h_vmem={cluster.mem} -q {cluster.queue}"
# echo "snakemake -s ~/scripts/pipelines/ImputationPipeline-snakemake/Snakefile -p -r --configfile /home/cocca/analyses/test_imputation/config_test_2.yaml --cores 32 " | qsub -N test_imputation -o ./
```

