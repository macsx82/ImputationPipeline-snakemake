# Imputation Pipeline

Implementation of an imputation workflow using a Snakemake pipeline.

## Setting things up

In order to run the pipeline, there are some requirements to fullfill and some set up needs to be perfomed.
In this current version, the pipeline is tested and configured to be run on the APOLLO cluster
It can be deployed also on the [ORFEO cluster](https://orfeo-documentation.readthedocs.io/en/latest/), but it will require to manually specify the location of all software binaries in the provided config file.

### Required Software

The following software has to be installed system-wide, in a user-defined Conda environment or using the modules architecture (ORFEO cluster).
Some of the software is not available as conda package, but compiled binaries or source code can be downloaded via the provided links.

+ awk
+ sed
+ python3
+ R
+ git
+ snakemake
+ [bcftools](http://www.htslib.org/doc/)
+ [shapeit v2](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#download)
+ [shapeit v4](https://odelaneau.github.io/shapeit4/#documentation)
+ [impute v5](https://jmarchini.org/software/#impute-5)
+ [plink v1.9](https://www.cog-genomics.org/plink/1.9/)

**Before switching to a new version of each software/module, a test run should be performed to check that the expected output files are generated and that they are consistent with the previous production version.**

### Required python packages

In order to run the pipeline, the following python packages have to be installed in your conda environment:

+ errno
+ gzip 
+ import
+ io
+ itertools
+ matplotlib
+ multiprocessing
+ numpy
+ os
+ pandas
+ pathlib
+ psutil
+ re
+ scipy
+ snakemake
+ sys

### ORFEO/general set up
1. Install Snakemake via conda ([link](https://snakemake.readthedocs.io/en/stable/getting\_started/installation.html));
    ```bash
    conda create -c conda-forge -c bioconda -n snakemake snakemake pandas
    ```
2. Activate the environment you created
    ```bash
    conda activate snakemake
    ```

### Apollo set up

1. Add the global snakemake environment to your environment list:
    ```bash
    conda config --append envs_dirs /shared/software/conda/envs
    conda config --prepend envs_dirs ~/.conda/envs
    ```

2. Check that the environment is available to you (you should see an entry "snakemake_g" in the list)
    ```bash
    conda env list
    ```
3. Load the environment
    ```bash
    conda activate snakemake_g
    ```
---


## Resources setup

In order to run the pipeline on a new system, there are some preparatory steps to perform, in order to retrieve or generate the resources needed.

At the moment, all resources needed are already available on the APOLLO cluster, at the following locations:

```
  genetic_map_path: "/netapp/nfs/resources/gen_map/shapeit4" #path containing genetic maps for IMPUTATION and PHASING (SHAPEIT4 and IMPUTE5 version)
  ref_panel_base_folder: "/shared/resources/references/VCF" #folder path containing reference panel files
  allele_recode_file: "/netapp/nfs/resources/dbSNP/human_9606_b154_GRCh37p13/GCF_000001405.25.SNPS.dbSNP154.tab" #folder path containing an allele recode file
  ext_ref_annot_file: "/netapp/nfs/resources/dbSNP/human_9606_b154_GRCh37p13/GCF_000001405.25.vcf.gz" #external reference file to perform annotation on VCF files
  ref_fasta: "/shared/resources/hgRef/hg19/hg19_nochr.fasta" #reference fasta file


  genetic_map_path: "/netapp/nfs/resources/1000GP_phase3/impute" #shapeit2 genetic map format path, to be used in this rule ONLY
  genetic_map: "/shared/software/eagle/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz"
```

**All the resources and the reference panels for this pipeline, are aligned to GRCh37.**

---

## Pipeline configuration

Below the details on the content of the config file the user should personalize in order to run the pipeline.

### GENERAL INFO section

In this section the user should provide information about the cohort and the imputation panel used:

```
cohort_name: "Sample_POP" #study population name
pop_group: "Pop_group" #population ancestry [can be one of EUR, AFR, ASN, EAS, SAS]
ref_panel: "ref_pane_name" #name of the reference panel. Can be one of [IGRPv1, TGP3]

```

### CHROMOSOMES section

In this section the user define the chromosomes to impute. The pipeline will work splitting each chromosome in chunks, using the parameters defined in the rules section.

### PATHS section

In this section, the user has to define input files and output folders. 
Since the main input files are required to be in **PLINK MAP and PED format**, in the first part, the user has to specify:

+ the *input folder path*
+ the *input file prefix* with the **complete path** (i.e. /home/cocca/input_files/GENETIC_DATA_TEST, assuming that in the folder /home/cocca/input_files are present the files GENETIC_DATA_TEST.map, GENETIC_DATA_TEST.ped)
+ the *output folder*, in which the pipeline will create all files needed to generate the final imputation
+ the *release folder* path, in which the pipeline will copy the final imputation results
+ a *log folder*, usually with a path relative to the pipeline execution folder

```
  input_folder: "input_folder" #base folder containing input data
  input_file_prefix: "input file prefix" #this is the prefix with the complete path of the PLINK file to be used.
  output_folder: "output_folder" #base folder for output
  release_folder: "release_folder" #base folder for data release. This can be the same as the output folder
  log_folder: "log_folder" #base folder for output
```

In the second section, the user has to specify the path for some resources needed for the pipeline to run:

```
  genetic_map_path: "/netapp/nfs/resources/gen_map/shapeit4" #path containing genetic maps for IMPUTATION and PHASING (SHAPEIT4 and IMPUTE5 version)
  ref_panel_base_folder: "/shared/resources/references/VCF" #folder path containing reference panel files
  allele_recode_file: "/netapp/nfs/resources/dbSNP/human_9606_b154_GRCh37p13/GCF_000001405.25.SNPS.dbSNP154.tab" #folder path containing an allele recode file
  ext_ref_annot_file: "/netapp/nfs/resources/dbSNP/human_9606_b154_GRCh37p13/GCF_000001405.25.vcf.gz" #external reference file to perform annotation on VCF files
  ref_fasta: "/shared/resources/hgRef/hg19/hg19_nochr.fasta" #reference fasta file
  snp_array_update_allele_file: "snp_array_update_allele_file" #absolute path containing the update allele file specific for the snp array used
  scripts: "scripts_folder_path" #folder containing scripts needed in the workflow
```

It is important to specify the correct SNP array update allele file, with the parameter *snp_array_update_allele_file*, using the absolute path.
Resources for most of the SNP arrays used are located, on the APOLLO cluster, in :

```
/shared/resources/genotyping/
```

The **scripts** parameter has to be filled by the user, and it should be the absolute path of the "scripts" folder available in the pipeline main folder.
This parameter is needed to run some external scripts from the pipeline.


### RULES section

In this section are defined rules specific parameters. Most of the parameters are set to default values according to each software documentation.

The user should define absolute paths for the **temp** folders for rules **concatImputed** and **vcfFixRef**, while all other resources are already set to default values on the Apollo cluster.


### TOOLS section

Here the user can define the path of each binary used, if it is not present in the $PATH.

At the moment, the template config file contains path specific for the usage on Apollo cluster.

```
tools:
  shapeit: "/share/apps/bio/bin/shapeit"
  phasing_tool: "/shared/software/eagle/Eagle_v2.4.1/eagle"
  impute: "/shared/software/impute5_v1.1.5/impute5_1.1.5_static"
  chunker_tool: "/shared/software/impute5_v1.1.5/imp5Chunker_1.1.5_static"
  plink: "plink"
  bcftools: "/share/apps/bio/bin/bcftools"
```

---
## Running the pipeline

There are different ways to run the pipeline: **Local mode**, **Cluster mode** or **Single node mode**

**ALL THE EXAMPLES BELOW ARE TAILORED TO THE APOLLO CLUSTER SUBMISSION SYSTEM, FOR WHICH YOU HAVE TO SPECIFY THE NAME OF THE PARALLEL ENVIRONMENT AVAILABLE FOR A SELECTED QUEUE**

On the **Apollo cluster**, due to latency problems witht the NetApp storages, it is advisable to add the option **--latency-wait 100** to the snakemake command, to avoid incurring in random execution errors.

### Local mode


In Local mode, the pipeline is executed in an interactive shell session (locally or on a cluster) and all the rules are treated as processes that can be run sequentially or in parallel, depending on the resources provided. One example of a Local execution is:

```bash
conda activate snakemake

base_cwd=/<USER_DEFINED_PATH>/Imputation_run
log_folder=${base_cwd}/Log
mkdir -p ${log_folder}
cd ${base_cwd}

snakefile=/<USER_DEFINED_PATH>/ImputationPipeline-snakemake/Snakefile
configfile=/<USER_DEFINED_PATH>/Imputation_run/ImputationPipeline_pipeline.yaml
cores=24

snakemake -p -r -s ${snakefile} --configfile ${configfile} --keep-going --cores ${cores} --latency-wait 100

```

In this example we assumed we had 24 CPU available for our calculation



### Cluster mode

In cluster mode, the pipeline runs on a interactive shell (**screen or tmux**) and each rule is submitted as a job on the cluster.
One example of a Cluster execution, on the **APOLLO cluster**, is:

```bash
conda source snakemake

base_cwd=/<USER_DEFINED_PATH>/Imputation_run
log_folder=${base_cwd}/Log
mkdir -p ${log_folder}
cd ${base_cwd}

snakefile=/<USER_DEFINED_PATH>/ImputationPipeline-snakemake/Snakefile
configfile=/<USER_DEFINED_PATH>/Imputation_run/ImputationPipeline_pipeline.yaml

log_name=Imputation_pipeline_r01.log
stderr_name=Imputation_pipeline_r01.err
cores=24
queue=fast
parall_env=smp

snakemake -p -r -s ${snakefile} --configfile ${configfile} --keep-going --cluster "qsub -N {rule}_{config[cohort_name]} -V -cwd -pe ${parall_env} {threads} -o {log.stdout} -e {log.stderr} -l h_vmem={resources.mem_mb} -q ${queue}" --latency-wait 100 -j ${cores} 1> ${log_name} 2> ${stderr_name}

```

In this example we defined also the name for two additional log files, which will help to keep track of the pipeline execution. In this case, the **-j** option will define how many concurrent jobs are submitted on the cluster.

Since the cluster used for this example is the Apollo cluster (SGE queue manager), the user has also to define the name of the parallel environment available. In this specific case, for the queue **fast** the parallel environment name to be specified with the **-pe** option, is **"smp"**

### Single node mode

In Single node mode, the pipeline runs as a job on the cluster and all rules are treated as processes that can be run sequentially or in parallel, depending on the resources provided. Similar to the Local execution mode.
One example of a single node mode execution, on the **APOLLO cluster**, is:

```bash
conda source snakemake

base_cwd=/<USER_DEFINED_PATH>/Imputation_run
log_folder=${base_cwd}/Log
mkdir -p ${log_folder}
cd ${base_cwd}

snakefile=/<USER_DEFINED_PATH>/ImputationPipeline-snakemake/Snakefile
configfile=/<USER_DEFINED_PATH>/Imputation_run/ImputationPipeline_pipeline.yaml

log_name=Imputation_pipeline_r01.log
stderr_name=Imputation_pipeline_r01.err

cores=32
queue=fast
mem=250G
parall_env=smp

echo "cd ${base_cwd};conda source snakemake; snakemake -p -r -s ${snakefile} --configfile ${configfile} --cores ${cores} --keep-going --latency-wait 100" | qsub -N {rule}_{config[cohort_name]} -V -cwd -pe ${parall_env} ${cores} -o {log_name} -e {stderr_name} -l h_vmem=${mem} -q ${queue}
```

In this example we selected an entire cluster node on the "fast" queue of the Apollo cluster, defining the number of CPU (32), the parallel environment (smp) and the total amount of RAM required to run the pipeline (250G). We defined also the name for the two additional log files, to keep track of the pipeline execution.

---

## Known issues

### Chromosome X imputation

At the moment, this pipeline covers only imputation of autosomal chromosomes. In order to perform imputation on chrX it is advisable to use one of the imputation server available from [Michighan](https://imputationserver.sph.umich.edu/index.html#!) or [Sanger](https://imputation.sanger.ac.uk/). 

### PDF report creation

At the moment, the creation of the PDF reports is triggered at the beginning of the pipeline execution, resulting in rules exiting with error status. This issue will happen as long as the pipeline generates the files needed for the PDF report creation (all imputation results for all chromosomes) .

Until this issue is fixed, it is advisable to perform a final re-run of the pipeline after removing the folder **07.stats** .
