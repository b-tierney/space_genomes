# anvio pipeline

for file in bacterial_genomes/J*; do

anvi-script-reformat-fasta $file -o "$file"_reformated --simplify-names
anvi-gen-contigs-database -f "$file"_reformated -o "$file".db -n space_genomes -T 3
anvi-run-hmms -c "$file".db

###add line to get functional annotations cogs
done

anvi-gen-genomes-storage -e external_genomes_subset -o ap-subset-GENOMES.db

anvi-pan-genome -g ap-subset-GENOMES.db -n space_genomes_subset_anvi

### scratch work


#!/bin/bash
#anvi-script-reformat-fasta $1 -o "$1"_reformated --simplify-names
#anvi-gen-contigs-database -f "$1"_reformated -o "$1".db -n space_genomes -T 1

