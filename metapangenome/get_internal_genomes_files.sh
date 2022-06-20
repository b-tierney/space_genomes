#!/bin/bash

prefix=$1
collection="Genomes"

echo -e "name\tbin_id\tcollection_id\tprofile_db_path\tcontigs_db_path" > internal_genomes_"${prefix}".txt

for bin in $(cat genome_ids); do

echo -e "$bin\t$bin\t$collection\tap-MERGED-${prefix}/PROFILE.db\tap_isolate_genomes-"${prefix}".db" >> internal_genomes_"${prefix}".txt

done