#!/bin/bash

# merge snp data

ls *earth/snps.vcf | xargs vcfcombine > ref-earth_snps.vcf
ls *typestrain/snps.vcf | xargs vcfcombine > ref-typestrain_snps.vcf
ls *mt1/snps.vcf | grep -v GCA_013450165 | xargs vcfcombine > ref-mt1_snps.vcf
ls *otheriss/snps.vcf | xargs vcfcombine > ref-otheriss_snps.vcf

for file in ref-*; do

while read p; 
do 
		refsource=$(echo $file | cut -f1 -d_ | cut -f2 -d-) 
		val1=$(echo $p | cut -f1 | cut -f1 -d' ') 
		val2=$(echo $p | cut -f2 | cut -f2 -d' ') 
		sed -i "s/"${val2}"/"${val1}"/g" $file 
		sed -i "s/_"${refsource}"//g" $file 
done<mapping.tsv

done