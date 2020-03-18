#!/bin/bash
VCF_FILES=data/vcf_files/*.vcf
VCFGZ_FILES=data/vcf_files/*.vcf.gz
for f in $VCF_FILES
do
	if [ "$f" != "$VCF_FILES" ]; then
		bcftools view -q 0.05:minor "$f" | bgzip -c > ./data/temp_vcfs/"$(basename $f)"
		bcftools index -t ./data/temp_vcfs/"$(basename $f)"
	fi
done
for f in $VCFGZ_FILES
do
	if [ "$f" != "$VCFGZ_FILES" ]; then
		bcftools view -q 0.05:minor "$f" | bgzip -c > ./data/temp_vcfs/"$(basename $f)"
		bcftools index -t ./data/temp_vcfs/"$(basename $f)"
	fi
done
bcftools concat ./data/temp_vcfs/*.vcf.gz -o ./data/merged_vcf.vcf.gz -O z
bcftools index -f -t ./data/merged_vcf.vcf.gz
gatk SelectVariants -V ./data/merged_vcf.vcf.gz -fraction 0.01 -O ./data/final.vcf.gz
bcftools index -f -t ./data/final.vcf.gz
