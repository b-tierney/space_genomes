#!/bin/bash

### analyze snp data

mkdir snp_freqs

# get allele frequences
for file in */snps.vcf;

do 

filename=$(cat $file | echo 's/\//-/g')
vcftools --vcf snps.vcf --freq --out snp_freqs/"$filename"

# convert to bed file


# merge 