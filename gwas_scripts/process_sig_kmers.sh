#!/bin/bash

# get prokka output
while read p; do
	prokka --cpus 5 --outdir "$1"_prokka $1
done<$1

find . -name '*_prokka/*faa' | xargs cat > all_genes

# build reference dataset of genes from prokka output
mmseqs/bin/mmseqs easy-cluster all_genes clusterRes tmp --min-seq-id 1.0 --threads 10 -c 0.9 --cov-mode 1

# build diamond database
diamond makedb --db apit_seqs --in clusterRes_rep_seq.fasta

# align kmers to clustered refseqs
diamond blastx -q apit_earth_v_space_kmers_sig_seqs_pos.fa --threads 5 --db apit_seqs --out aligned_kmers_pos.txt --outfmt 6 
diamond blastx -q apit_earth_v_space_kmers_sig_seqs_neg.fa --threads 5 --db apit_seqs --out aligned_kmers_neg.txt --outfmt 6 

cut -f2 aligned_kmers_pos.txt > pos_genes
cut -f2 aligned_kmers_neg.txt > neg_genes

grep -Ff pos_genes prokka_annotation_data > apit_earth_v_space_kmers_sig_seqs_pos_annotations.tsv
grep -Ff neg_genes prokka_annotation_data > apit_earth_v_space_kmers_sig_seqs_neg_annotations.tsv

