#!/usr/bin/env bash
#service script to merge vcfs post imputation and generate relevant stats, on Apollo cluster
# B1_basef=$1
# B2_basef=$2
# B1_vcf=${B1_basef}/${chr}/${chr}.vcf.gz
# B2_vcf=${B2_basef}/${chr}/${chr}.vcf.gz
B1_vcf=$1
B2_vcf=$2
outfolder=$3
chr=$4
prefix=$5

mkdir -p ${outfolder}/CALLRATE95
mkdir -p ${outfolder}/MERGED
mkdir -p ${outfolder}/03.IMPUTED/VCF
mkdir -p ${outfolder}/03.IMPUTED/BIMBAM
mkdir -p ${outfolder}/04.STATS

#we need a conda env containing bcftools v1.14 and to correctly setup the plugin folder!
source activate bioinfo

outmerged=${outfolder}/MERGED/${prefix}_${chr}.vcf.gz


echo "Merging batches..."
bcftools merge -m none -i INFO:join ${B1_vcf} ${B2_vcf} -O z -o ${outmerged}
bcftools index -t ${outmerged}
echo "Merging done"

echo "Filtering call rate 95%"
bcftools +fill-tags ${outmerged} -- -t all,F_MISSING | bcftools view -e "F_MISSING>=0.05"| bcftools view -e "HWE <= 0.000001" -O z -o ${outfolder}/CALLRATE95/${prefix}_${chr}.vcf.gz
bcftools index -t ${outfolder}/CALLRATE95/${prefix}_${chr}.vcf.gz
echo "Filtering done"

#Generate info stats with qctools: calculate the impute info scores
echo "Calculate stats and index"
#we need to reformat a little bit the resulting table
(echo -e "#CHROM\tPOS\tID\tREF\tALT\tR2";fgrep -v "#" ${outfolder}/CALLRATE95/${prefix}_${chr}.info_stats| fgrep -v "alternate_ids" | cut -f 2-6,18 | awk 'BEGIN{FS=" ";OFS="\t"}{if($6=="NA") print $2,$3,$1,$4,$5,1;else print $2,$3,$1,$4,$5,sprintf("%.3f",$6)}') | bgzip -c > ${outfolder}/CALLRATE95/${prefix}_${chr}.info.tab.gz
tabix -f -s 1 -b 2 -e 2 ${outfolder}/CALLRATE95/${prefix}_${chr}.info.tab.gz
echo "stats done"

#add the info stats as R2 field
echo "Annotate and index"
#create file for header annotation
echo "##INFO=<ID=R2,Number=.,Type=Float,Description=\"Imputation Info score calculated using qctool software after merging different batches\">" > ${outfolder}/CALLRATE95/${prefix}_${chr}.header
bcftools annotate -c CHROM,POS,-,REF,ALT,R2 -h ${outfolder}/CALLRATE95/${prefix}_${chr}.header -a ${outfolder}/CALLRATE95/${prefix}_${chr}.info.tab.gz ${outfolder}/CALLRATE95/${prefix}_${chr}.vcf.gz -O z -o ${outfolder}/03.IMPUTED/VCF/${chr}.vcf.gz
bcftools index -t ${outfolder}/03.IMPUTED/VCF/${chr}.vcf.gz
rm ${outfolder}/CALLRATE95/${prefix}_${chr}.header
echo "Done"

echo "Generate bcftools stats for comparison with the splitted imputation results"
bcftools stats ${outfolder}/03.IMPUTED/VCF/${chr}.vcf.gz > ${outfolder}/03.IMPUTED/VCF/${chr}.merge_stats
echo "Done"

echo "Generate table style stats after merging"
bcftools query -H -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%INFO/R2\t%INFO/IMP\n" ${outfolder}/03.IMPUTED/VCF/${chr}.vcf.gz| awk 'BEGIN{FS="\t";OFS="\t"}{if($8==".") print $1,$2,$3,$4,$5,$6,$7,0;else print $0}' | gzip -c > ${outfolder}/04.STATS/{chr}.info_stats.gz
echo "Done"

echo "Create the bimbam version"
#We just need to format the vcf file and extract the DS field already present in the VCF, for the main bimbam
bcftools query -f'%ID,%REF,%ALT[,%DS]\\n' ${outfolder}/03.IMPUTED/VCF/${chr}.vcf.gz | gzip --best -c > ${outfolder}/03.IMPUTED/BIMBAM/${chr}.bimbam.gz
#We just need to format the vcf file and extract the position
bcftools query -f'%ID,%POS,%CHROM\\n' ${outfolder}/03.IMPUTED/VCF/${chr}.vcf.gz -o ${outfolder}/03.IMPUTED/BIMBAM/${chr}.pos
echo "Done"