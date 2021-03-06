#!/bin/bash

pyseer --lmm --phenotypes pyseer_metadata.tsv --phenotype-column EARTH_SPACE_BINARY --kmers unitig_output/unitigs.txt.gz --similarity distances_with_bins --output-patterns kmer_patterns_apit_isolation_source.txt --cpu 10 > apit_earth_v_space_kmers.txt
pyseer --lmm --phenotypes ../pyseer_metadata.tsv --phenotype-column FLIGHT_BINARY --vcf ref-mt1_snps.vcf --similarity ../distances_names --cpu 10 > apit_mt1-mt2_snps.txt &>mt1_mt2_earth_space_log
pyseer --lmm --phenotypes ../pyseer_metadata.tsv --phenotype-column EARTH_SPACE_BINARY --vcf ref-earth_snps.vcf --similarity ../distances_names  --cpu 10 > apit_earth_v_space_snps.txt &> snp_earth_space_log
pyseer --lmm --phenotypes ../pyseer_metadata.tsv --phenotype-column EARTH_SPACE_BINARY --vcf ref-typestrain_snps.vcf --similarity ../distances_names --cpu 10 > apit_typestrain_v_space_snps.txt &> typestrain_earth_space_log


