#!/bin/bash

# build prokka id mapping file

find . -type d | cut -f2 -d/ | grep prokka  > prokkadirs

rm prokka_id_mapping
while read p;

do 

name=$(echo $p | sed 's/.fa_prokka//g');
prokkaid=$(head -2 "${p}"/*tsv  | tail -n 1 | cut -f1 -d_);

echo -e ""$name"\t"$prokkaid"" >> prokka_id_mapping

done<prokkadirs