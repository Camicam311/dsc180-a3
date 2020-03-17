#!/bin/bash
VCF_FILES=data/vcf_files/*.vcf
VCFGZ_FILES=data/vcf_files/*.vcf.gz
for f in $VCF_FILES
do
	bcftools view -q 0.05:minor "$f" | bgzip -c > ./data/final_vcfs/"$(basename $f)"
	bcftools index -t ./data/final_vcfs/"$(basename $f)"
done
for f in $VCFGZ_FILES
do
	bcftools view -q 0.05:minor "$f" | bgzip -c > ./data/final_vcfs/"$(basename $f)"
	bcftools index -t ./data/final_vcfs/"$(basename $f)"
done
bcftools merge ./data/final_vcfs/*.vcf.gz | bgzip -c > ./data/merged_vcf.vcf
