#!/bin/bash

# download all sra files for a given set of assembly ids

while read p <&3;

do 

sra=$(esearch -db assembly -query $p | elink -target biosample | elink -target sra | efetch -format docsum | xtract -pattern DocumentSummary -ACC @acc -block DocumentSummary -element "&ACC"  | cut -f5)

echo $sra
prefetch $sra
fastq-dump --split-files $sra

done 3<$1