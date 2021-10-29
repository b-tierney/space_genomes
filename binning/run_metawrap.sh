#!/bin/bash

basename=$1
r1=$2
r2=$3

cp "$r1" "$basename"_1.fastq.gz 
gunzip -f "$basename"_1.fastq.gz 

cp "$r2" "$basename"_2.fastq.gz 
gunzip -f "$basename"_2.fastq.gz 

metawrap binning --run-checkm -a $1 --metabat2 --maxbin2 --concoct -o "$basename"_binning "$basename"_1.fastq "$basename"_2.fastq

rm "$basename"_1.fastq "$basename"_2.fastq