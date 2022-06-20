#!/bin/bash

### run freebayes

bamfile=$1
reference=$2
referencebasename=$3

samtools faidx $reference

freebayes-parallel <(fasta_generate_regions.py "${reference}".fai 100000) 5 -f ${reference} --ploidy 1 --report-monomorphic --bam-list $bamfile --vcf "$referencebasename".vcf
