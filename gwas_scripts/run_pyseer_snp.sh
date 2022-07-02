#!/bin/bash

# run pyseer and parse for phandango

bcftools norm -m - out.filtered.vcf > out.split.filtered.vcf

pyseer --phenotypes pyseer_snpmetadata.tsv --vcf out.split.filtered.vcf --lmm --similarity sim  --cpu 4 > mgwas_assoc

cut -f1 -d_ mgwas_assoc | sed '1d' > contigs
cut -f3- -d_ mgwas_assoc | cut -f1 sed '1d'  > snp
cut -f2 -d_ mgwas_assoc |sed '1d' > bp 
cut -f4 mgwas_assoc | sed '1d' | awk '{print -log($1)/log(10),-log($1)/log(10),"0"}' | sed -e 's/ /\t/g'> log10pval


echo "#CHR SNP BP minLOG10(P) log10(p) r^2" > names

paste contigs bp log10pval > foo
cat names foo > mgwas_assoc.plot
rm foo names contigs bp
