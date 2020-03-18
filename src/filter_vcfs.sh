#!/bin/bash
bcftools concat ./data/vcf_files/*.vcf.gz -o ./data/merged_vcf.vcf.gz -O z
bcftools index -f -t ./data/merged_vcf.vcf.gz
gatk SelectVariants -V ./data/merged_vcf.vcf.gz -fraction 0.01 -O ./data/final.vcf.gz
bcftools index -f -t ./data/final.vcf.gzx
