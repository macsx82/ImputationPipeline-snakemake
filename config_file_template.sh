#!/usr/bin/env bash
#
# Wrapper to create the config fle for snakemake imputataion pipeline
# This wersion handle a single chr at a the time
# next version will generate a template with all chr to process

#Define the template file creation function
function build_template(){
pop=$1
pop_group=$2
chr=$3
input_folder=$4
output_folder=$5
ref_panel=$6 #"IGRPv1"
ref_panel_base_dir=$7 #"/netapp/dati/WGS_REF_PANEL/REFERENCES"

cat << EOF
#config file in yaml format for the IPF pipeline
#
#
pop:
  "${pop}"
pop_group:
  "${pop_group}"
#we need to start with the vcf files input for the phasing step
chr:
  "${chr}"

input_folder:
  # "/netapp/dati/INGI_WGS/18112015/CARL/12112015_FILTERED_REL/30092016_CONV_ID/07022018_VEP90_CSQ"
  #"/home/cocca/analyses/imputation/03032020/ITT/18"
  "${input_folder}"
genetic_map_chr:
  "/netapp/nfs/resources/1000GP_phase3/impute/genetic_map_chr{chr_to_phase}_combined_b37.txt"
output_folder:
  "${output_folder}"
ref_panel:
  "${ref_panel}"
ref_panel_base_folder:
  "${ref_panel_base_dir}"

shapeit_path:
  "/share/apps/bio/bin/shapeit"

EOF

}

#generate the template with command line options

script_dir=`dirname $0`

if [ $# -lt 1 ]
then
    echo "#########################"
    echo "WRONG argument number!"
    echo "Usage:"
    echo "config_file_template.sh -i <input_file_folder> -t <template_folder> -o <output_folder> -c <chr> -p <pop_name> -g <pop_group> -r <ref_panel_name> -f <ref_panel_base_folder>"
    echo "Execution options:-i <input_file_folder>"
    echo "					-t <template_folder>"
    echo " 					-o <output_folder>"
    echo "  				-c <chr>"
    echo "  				-p <pop_name>"
    echo "  				-g <pop_group>"
    echo "  				-r <ref_panel_name>"
    echo "  				-f <ref_panel_base_folder>"
    echo "                  -h: this help message "
    echo "#########################"
    exit 1
fi

#ARGS
# template_dir=$1
# out_dir=$2

suffix=`date +"%d%m%Y%H%M%S"`
runner_mode=()

echo "${@}"
while getopts ":i:t:o:c:p:g:r:f:h" opt ${@}; do
  case $opt in
    i)
      echo ${OPTARG}
      input_file_folder=${OPTARG}
    ;;
    t)
      echo ${OPTARG}
      template_dir=${OPTARG}
      ;;
    o)
      echo ${OPTARG}
      out_dir=${OPTARG}
      ;;
    c)
      echo ${OPTARG}
      chr=${OPTARG}
      ;;
    p)
      echo ${OPTARG}
      pop=${OPTARG}
      ;;  
    g)
      echo ${OPTARG}
      pop_group=${OPTARG}
    ;;
    r)
      echo ${OPTARG}
      ref_panel=${OPTARG}
    ;;
    f)
      echo ${OPTARG}
      ref_panel_base_folder=${OPTARG}
    ;;
    h)
    echo "#########################"
    echo "WRONG argument number!"
    echo "Usage:"
    echo "config_file_template.sh -i <input_file_folder> -t <template_folder> -o <output_folder> -c <chr> -p <pop_name> -g <pop_group> -r <ref_panel_name> -f <ref_panel_base_folder>"
    echo "Execution options:-i <input_file_folder>"
    echo "					-t <template_folder>"
    echo " 					-o <output_folder>"
    echo "  				-c <chr>"
    echo "  				-p <pop_name>"
    echo "  				-g <pop_group>"
    echo "  				-r <ref_panel_name>"
    echo "  				-f <ref_panel_base_folder>"
    echo "                  -h: this help message "
    echo "#########################"
    exit 1
	;;
    *)
      echo $opt
    ;;
  esac

done
#Create folders
mkdir -p ${template_dir} ${out_dir}

#generate the template
build_template ${pop} ${pop_group} ${chr} ${input_file_folder} ${out_dir} ${ref_panel} ${ref_panel_base_folder} > ${template_dir}/config_file_${pop}_${pop_group}_${chr}.yml