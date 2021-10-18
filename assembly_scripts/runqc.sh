#!/bin/bash

# qc paired end, non interleaved metagenomes with bbtools

# you'll arguments are the base file name (eg sample_1) followed by the two paired end files.

# note absolute paths in bowtie2 and bbduk commands.

# adapters.fa file comes from the "atlas metagenome" workflow, available on github

basename=$1
read1=$2
read2=$3

echo "$basename"_log
echo $1

touch "$basename"_log

echo $1 > "$basename"_log

echo 'STARTING QC' >> "$basename"_log

clumpify.sh in=$read1 \
            in2=$read2 \
            -Xmx60G \
            out="${basename}"_R1_clumped.fastq.gz \
            out2="${basename}"_R2_clumped.fastq.gz \
            dedupe=t \
            dupesubs=2 \
            qin=33 \
            threads=10 \
            overwrite=true \
            optical=f 2>> "$basename"_log

echo 'CLUMPIFY DONE' >> "$basename"_log

bbduk.sh in="${basename}"_R1_clumped.fastq.gz \
           -Xmx60G \
           in2="${basename}"_R2_clumped.fastq.gz  \
           out="${basename}"_R1_clumped_bbduk.fastq.gz \
           out2="${basename}"_R2_clumped_bbduk.fastq.gz \
           interleaved=f \
           stats="${basename}"_stats_bbduk \
           overwrite=true \
           qout=33 \
           trd=t \
           hdist=1 \
           k=27 \
           ktrim="r" \
           mink=8 \
           overwrite=true \
           trimq=10 \
           qtrim='rl' \
           threads=10 \
           minlength=51 \
           maxns=-1 \
           minbasefrequency=0.05 \
           ecco=f \
           prealloc=t \
           ref=/athena/masonlab/scratch/users/btt4001/atlas_databases/adapters.fa 2>> "$basename"_log


rm -rf "${basename}"_R1_clumped.fastq.gz "${basename}"_R2_clumped.fastq.gz

echo 'BBDUK DONE' >> "$basename"_log

echo 'RUNNING BOWTIE' >> "$basename"_log
bowtie2 -p 10 \
        -x /athena/masonlab/scratch/users/btt4001/hg38/GRCh38_noalt_as/GRCh38_noalt_as \
        -1 "${basename}"_R1_clumped_bbduk.fastq.gz \
        -2 "${basename}"_R2_clumped_bbduk.fastq.gz \
        --al-conc-gz "${basename}"_human_reads \
        --very-sensitive-local \
        --un-conc-gz "${basename}"_SAMPLE_hg38removed > "${basename}"_SAMPLE_mapped_and_unmapped.sam 2>> "$basename"_log

mv "${basename}"_SAMPLE_hg38removed.1 "${basename}"_R1_clumped_bbduk_hg38removed.fastq.gz
mv "${basename}"_SAMPLE_hg38removed.2 "${basename}"_R2_clumped_bbduk_hg38removed.fastq.gz

echo 'BOWTIE DONE' >> "$basename"_log

rm -rf "${basename}"_R1_clumped_bbduk.fastq.gz "${basename}"_R2_clumped_bbduk.fastq.gz

repair.sh in="${basename}"_R1_clumped_bbduk_hg38removed.fastq.gz \
          in2="${basename}"_R2_clumped_bbduk_hg38removed.fastq.gz \
          out1="${basename}"_R1_clumped_bbduk_hg38removed_repaired.fastq.gz \
          out2="${basename}"_R2_clumped_bbduk_hg38removed_repaired.fastq.gz 2>> "$basename"_log

rm -rf "${basename}"_R1_clumped_bbduk_hg38removed.fastq.gz "${basename}"_R2_clumped_bbduk_hg38removed.fastq.gz

echo 'REPAIRING DONE' >> "$basename"_log

tadpole.sh prealloc=1 \
           -Xmx60G \
           in="${basename}"_R1_clumped_bbduk_hg38removed_repaired.fastq.gz \
           in2="${basename}"_R2_clumped_bbduk_hg38removed_repaired.fastq.gz \
           out="${basename}"_R1_clumped_bbduk_hg38removed_repaired_tadpole.fastq.gz \
           out2="${basename}"_R2_clumped_bbduk_hg38removed_repaired_tadpole.fastq.gz \
           mode=correct \
           threads=10 \
           ecc=t \
           ecco=t 2>> "$basename"_log

rm -rf "${basename}"_R1_clumped_bbduk_hg38removed_repaired.fastq.gz "${basename}"_R2_clumped_bbduk_hg38removed_repaired.fastq.gz

echo 'TADPOLE DONE' >> "$basename"_log

echo 'QC COMPLETE' >> "$basename"_log
