#!/bin/bash

pyseer --lmm --phenotypes metadata_filtered.tsv --phenotype-column Isolation_Source --kmers output2.txt.gz --similarity distances_apit --output-patterns kmer_patterns_apit_isolation_source.txt --cpu 10 > apit_isolation_source_kmers.txt

pyseer --lmm --phenotypes metadata_filtered.tsv --phenotype-column Project --kmers output2.txt.gz --similarity distances_apit --output-patterns kmer_patterns_apit_isolation_source.txt --cpu 10 > apit_project_kmers.txt

