#!/bin/bash

# need to index ref genome first

basename=$1
read1=$2
read2=$3
reference=$4
referencebasename=$5

bwa mem -M -t 10 $reference $read1 $read2 | samtools view -buS - > "${basename}"_"${referencebasename}".bam

bam="${basename}"_"${referencebasename}".bam
name="${basename}"_"${referencebasename}"

java -jar ~/picard/picard.jar SortSam \
      I=${bam} \
      O=${bam%.*}-sorted.bam \
      SORT_ORDER=coordinate

java -jar ~/picard/picard.jar MarkDuplicates \
      I=${bam%.*}-sorted.bam \
      O=${bam%.*}-sorted-md.bam \
      M=${bam%.*}-md-metrics.txt

java -jar ~/picard/picard.jar AddOrReplaceReadGroups \
      I=${bam%.*}-sorted-md.bam \
      O=${bam%.*}-sorted-md-rg.bam \
      RGID=${bam%.*} \
      RGLB=${name} \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=${name}

samtools index ${bam%.*}-sorted-md-rg.bam

