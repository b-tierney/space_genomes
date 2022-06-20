#!/bin/bash

prefix="ap"
readsStart=11  # which column in anvio layers table is the start of the metagenomes

sites=(MT1ENV MT2ENV MT2CREW TWINS)  # important that td is first if that will be the one onto which others are added

#mkdir temp
#cp *MT1ENV* temp

#anvi-delete-misc-data -p temp/$prefix-PAN-$site/$prefix-PAN-$site-PAN.db  -t layers --keys-to-remove $(head -1 ap-MT1ENV-layers.txt | cut -f$readsStart- - | tr '\t' ',')

for site in ${sites[*]}; do

anvi-export-misc-data -p $prefix-PAN-$site/$prefix-PAN-$site-PAN.db -t items -o $prefix-$site-items.txt
anvi-export-misc-data -p $prefix-PAN-$site/$prefix-PAN-$site-PAN.db -t layers -o $prefix-$site-layers.txt

# delete the unprocessed columns from TD so don't end up with two sets of TD metagenomes (formatted + unformatted)

# select metapan rings from items data and add oral habitat as prefix
awk -F"\t" -v site=$site 'BEGIN{site=toupper(site)}; NR==1{print $1 FS site"-"$9 FS site"-"$10 FS site"-"$11 FS site"-"$12} NR>1 {print $1 FS $9 FS $10 FS $11 FS $12}' $prefix-$site-items.txt > $prefix-$site-items.tmp
# prepend capitalized site ID before each metagenome's coverage of the genomes
cut -d' ' -f1,$readsStart- $prefix-$site-layers.txt | sed -e "s/^/\U$site-/g" > $prefix-$site-layers.tmp

# import this data back into the mt1 data
anvi-import-misc-data -p temp/ap-PAN-MT1ENV/ap-PAN-MT1ENV-PAN.db -t items $prefix-$site-items.tmp
anvi-import-misc-data -p temp/ap-PAN-MT1ENV/ap-PAN-MT1ENV-PAN.db  -t layers $prefix-$site-layers.tmp

done