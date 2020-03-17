#!/bin/bash
VCF_FILES=data/vcf_files/*.vcf
VCFGZ_FILES=data/vcf_files/*.vcf.gz
for f in $VCF_FILES
do
	if [ "$f" != "$VCF_FILES" ]; then
		bcftools view -q 0.05:minor "$f" | bgzip -c > ./data/final_vcfs/"$(basename $f)"
		bcftools index -t ./data/final_vcfs/"$(basename $f)"
	fi
done
for f in $VCFGZ_FILES
do
	if [ "$f" != "$VCFGZ_FILES" ]; then
		bcftools view -q 0.05:minor "$f" | bgzip -c > ./data/final_vcfs/"$(basename $f)"
		bcftools index -t ./data/final_vcfs/"$(basename $f)"
	fi
done
bcftools concat ./data/final_vcfs/*.vcf.gz -o ./data/merged_vcf.vcf -O v
bcftools index -t ./data/merged_vcf.vcf
