# notes on anvio metapangenome test run

var=( {a..z}{a..z}{a..z} )
i=0

for file in *fna;

do

i=$((i+1))

alphavar=${var[i]}

anvi-script-reformat-fasta --seq-type NT --simplify-names --prefix $alphavar -o "$file"_reformatted $file  

echo $file $alphavar >> metapangenome_contig_prefix_mapping

done

cat *_reformatted  > all_contigs.fna

anvi-gen-contigs-database -T 8 -f all_contigs.fna -o contigs.db -n 'PDS08 isolate genomes'

# run hmms
anvi-run-hmms -c contigs.db --num-threads 4

# run cogs
anvi-run-ncbi-cogs -c contigs.db -T 4

# generate bowtie index
bowtie2-build all_contigs.fna contigs_bowtie

# run bowtie_samples.sh 
while read p; do
    sample=$(echo $p | cut -f1 -d' ')
    r1=$(echo $p | cut -f2 -d' ')
    r2=$(echo $p | cut -f3 -d' ')    
    if [ "$sample" == "sample" ]; then continue; fi
    echo $sample
    bowtie2 --threads 4  -x contigs_bowtie -1 "$r1" -2 "$r2" --no-unal -S "$sample".sam
    samtools view -F 4 -bS "$sample".sam > "$sample"-RAW.bam
    samtools sort "$sample"-RAW.bam -o "$sample".bam
    samtools index "$sample".bam
    rm $sample.sam "$sample"-RAW.bam
done<samples.txt

# anvio_profile.sh
for sample in `awk '{print $1}' samples.txt`
do
    if [ "$sample" == "sample" ]; then continue; fi

    anvi-profile -c contigs.db \
                 -i $sample.bam \
                 -M 100 \
                 --skip-SNV-profiling \
                 --num-threads 10 \
                 -o $sample
done

### get split genome mapping
for split_name in `sqlite3 contigs.db 'select split from splits_basic_info;'`
do
    GENOME=`echo $split_name | cut -f1 -d_`
    echo -e "$split_name\t$GENOME"
done > contigs.txt

cut -f2 contigs.txt | sort | uniq > genome_ids


### create 4 merged groups, one for each meatgenome class
anvi-merge */PROFILE.db \
           -o pds08-MERGED\
           -c contigs.db


# IMPORT -- SUMMARIZE -- PANGENOME -- METAPANGENOME

anvi-import-collection contigs.txt \
                       -c contigs.db \
                       -p pds08-MERGED/PROFILE.db \
                       -C Genomes

anvi-summarize -c contigs.db \
               -p pds08-MERGED/PROFILE.db \
               -C Genomes \
               --init-gene-coverages \
               -o pds08-SUMMARY

./get_internal_genome_files.sh

anvi-gen-genomes-storage -i internal_genomes.txt \
                         -o pds08-PAN-GENOMES.db

anvi-pan-genome -g pds08-PAN-GENOMES.db \
                --use-ncbi-blast \
                --minbit 0.5 \
                --mcl-inflation 10 \
                --project-name pds08-PAN \
                --num-threads 20

anvi-meta-pan-genome -i internal_genomestxt \
                                           -p pds08-PAN/pds08-PAN-PAN.db \
                                           -g pds08-PAN-GENOMES.db \
                                           --fraction-of-median-coverage 0.25   


### lookin at genes across database

while read p; do 

anvi-script-gen-distribution-of-genes-in-a-bin -p pds08-MERGED/PROFILE.db \
                                               -c contigs.db \
                                               -C Genomes \
                                               --fraction-of-median-coverage 0.25 \
                                               -b "$p" 

done<genome_ids

# generate gene-level distribution data for MT1 genome:
anvi-script-gen-distribution-of-genes-in-a-bin -p pds08-MERGED/PROFILE.db \
                                               -c contigs.db \
                                               -C Genomes \
                                               --fraction-of-median-coverage 0.25 \
                                               -b aag 

# do the same for MT2 genome:
anvi-script-gen-distribution-of-genes-in-a-bin -p pds08-MERGED/PROFILE.db \
                                               -c contigs.db \
                                               -C Genomes \
                                               --fraction-of-median-coverage 0.25 \
                                               -b aas

# and type strain

# do the same for MT2 genome GCA_017166275:
anvi-script-gen-distribution-of-genes-in-a-bin -p pds08-MERGED/PROFILE.db \
                                               -c contigs.db \
                                               -C Genomes \
                                               --fraction-of-median-coverage 0.25 \
                                               -b aab      


anvi-interactive -d aab-GENE-COVs.txt \
                 -A aab-ENV-DETECTION.txt \
                 --manual \
                 -p GENES-PROFILE.db \
                 --title "PDS08-1 genes MT1/MT2 metagenomes"


anvi-interactive -d aag-GENE-COVs.txt \
                 -A aag-ENV-DETECTION.txt \
                 --manual \
                 -p GENES-PROFILE.db \
                 --title "AP MT1 genes MT1/MT2 metagenomes"

anvi-interactive -d aas-GENE-COVs.txt \
                 -A aas-ENV-DETECTION.txt \
                 --manual \
                 -p GENES-PROFILE.db \
                 --title "AP MT2 genes MT1/MT2 metagenomes"














