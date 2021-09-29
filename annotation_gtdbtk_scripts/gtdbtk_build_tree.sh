#!/bin/bash

# build inferred tree

# may need to change gtdbtk data path, requires pyseer in working directory

GTDBTK_DATA_PATH=/athena/masonlab/scratch/users/btt4001/gtdbtk_refs/release202

genomes_to_grab=$1
bacterial_genome_locs=$2

grep -Ff $genomes_to_grab $bacterial_genome_locs > genome_config_file_"${genomes_to_grab}"
paste genome_config_file_"${genomes_to_grab}" $genomes_to_grab > temp
mv temp genome_config_file_"${genomes_to_grab}"

gtdbtk identify --cpus 10 --batchfile genome_config_file_"${genomes_to_grab}"  --out_dir gtdbk_identify_"${genomes_to_grab}"
gtdbtk align --cpus 10 --identify_dir gtdbk_identify_"${genomes_to_grab}" --out_dir gtdbk_align_"${genomes_to_grab}"
gtdbtk infer --cpus 10 --msa_file gtdbk_align_"${genomes_to_grab}"/gtdbtk.bac120.user_msa.fasta --out_dir gtdbk_infer_"${genomes_to_grab}"

python pyseer/scripts/phylogeny_distance.py gtdbk_infer_"${genomes_to_grab}"/gtdbtk.unrooted.tree > distances_"${genomes_to_grab}"
