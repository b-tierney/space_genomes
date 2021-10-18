#!/bin/bash

# generate contig bins

contigpath=$1
rawreads1=$2
rawreads2=$3

# run bowtie2

bowtie2-build "$contigpath" "$contigpath"

bowtie2 -p 5 -x "$contigpath" --very-sensitive-local -1 "$rawreads1" -2 "$rawreads2" | samtools view -bS | samtools sort > "$contigpath"_bowtieout.bam
samtools index "$contigpath"_bowtieout.bam

# run metabat2
metabat2 -i "$contigpath"  -a "$contigpath"_bowtieout.bam -o "$contigpath"_binned --threads 5 --minContig 1500