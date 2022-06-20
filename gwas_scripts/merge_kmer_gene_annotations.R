# merge kmer - gene - annotation output

library(tidyverse)
library(Biostrings)


parse_kmer_findings <- function(association_output,alignment_output,annotation_data,sequences){
	data_pos = read.table(association_output,header=T,sep='\t')
	data_pos = data_pos %>% select(variant,af,beta,by)
	alignment_pos = read.table(alignment_output,header=F,sep='\t')
	alignment_pos = alignment_pos %>% select(V1,V2)
	colnames(alignment_pos) = c('kmer','gene_id')
	annotation_data = read.csv(annotation_data,sep='\t',header=T)
	fastaFile <- readDNAStringSet(sequences)
	seq_name = names(fastaFile)
	sequence = paste(fastaFile)
	seq_mapping <- data.frame(seq_name, sequence) %>% mutate(length = nchar(sequence)) %>% filter(nchar >= 21)

	data_pos = left_join(data_pos,seq_mapping,by=c('variant'='sequence'))
	alignment_pos = left_join(alignment_pos,data_pos,by=c('kmer'='seq_name'))
	merged_data = left_join(alignment_pos,annotation_data,by=c('gene_id'='locus_tag'))

	merged_data_means_gene = merged_data %>% group_by(gene_id,product,COG,EC_number,gene) %>% summarise(beta = mean(beta),af = mean(af), by = mean(by))  %>% arrange(by,desc(beta),desc(af))
	merged_data_COG = merged_data %>% group_by(COG)%>% summarise(beta = mean(beta),af = mean(af), by = mean(by)) %>% arrange(by,desc(beta),desc(af))
	merged_data_EC_number = merged_data %>% group_by(EC_number)%>% summarise(beta = mean(beta),af = mean(af), by = mean(by))  %>% arrange(by,desc(beta),desc(af))
	merged_data_product = merged_data %>% group_by(product)%>% summarise(beta = mean(beta),af = mean(af), by = mean(by))  %>% arrange(by,desc(beta),desc(af))
	merged_data_gene = merged_data %>% group_by(gene)%>% summarise(beta = mean(beta),af = mean(af), by = mean(by))  %>% arrange(by,desc(beta),desc(af))
	return(list(merged_data_means_gene,merged_data_COG,merged_data_EC_number,merged_data_product,merged_data_gene))
}


positive_associations = parse_kmer_findings('apit_earth_v_space_kmers_sig_pos.tsv','aligned_kmers_pos.txt','prokka_annotation_data', 'apit_earth_v_space_kmers_sig_seqs_pos.fa')
negative_associations = parse_kmer_findings('apit_earth_v_space_kmers_sig_neg.tsv','aligned_kmers_neg.txt','prokka_annotation_data', 'apit_earth_v_space_kmers_sig_seqs_neg.fa')



