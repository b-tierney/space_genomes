#!/bin/bash

### analyze snp data

mkdir snp_freqs

# get allele frequences
for file in */snps.vcf;

do 

filename=$(echo $file | sed 's/\//-/g')
#echo $filename
vcftools --vcf $file --freq --out snp_freqs/"$filename"
awk -F'\t' '{print $1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' snp_freqs/"$filename" > snp_freqs/"$filename".bed

done 

# convert to bed file


# merge 





# or convert merged vcf to bcf
bgzip -c ref-mt1_snps.vcf > ref-mt1_snps.vcf.gz
tabix -p vcf ref-mt1_snps.vcf.gz
bcftools view -O b -o ref-mt1_snps.bcf ref-mt1_snps.vcf.gz
bcftools query -H -f '%CHROM %POS[\t%DP\t%RO\t%AO\t]\n' ref-mt1_snps.bcf > ref-mt1_snps_for_heatmap.tsv
