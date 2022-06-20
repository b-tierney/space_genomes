# parse kmer association output

library(tidyverse)

data = read.table('apit_earth_v_space_kmers.txt',sep='\t',header=T)

seqs = data %>% select(variant) %>% unlist %>% unname

output = list()
count = 0
for(seq in seqs_pos){
	count = count + 1
	output[[count]] = paste('>kmer_',count,sep='')
	count = count + 1
	output[[count]] = seq
}

output = output %>% as.data.frame %>% t
write.table(output,'all_kmers.fa',row.names=F,col.names=F,quote=F)

data_sub = data %>% mutate(by = p.adjust(lrt.pvalue,method='BY')) %>% filter(by <0.05)
write.table(data_sub %>% filter(beta>=0),'apit_earth_v_space_kmers_sig_pos.tsv',row.names=F,sep='\t')
write.table(data_sub %>% filter(beta<0),'apit_earth_v_space_kmers_sig_neg.tsv',row.names=F,sep='\t')

seqs_pos = data_sub %>% filter(beta>=0) %>% select(variant) %>% unlist %>% unname
seqs_neg = data_sub %>% filter(beta<0)%>% select(variant) %>% unlist %>% unname

output = list()
count = 0
for(seq in seqs_pos){
	count = count + 1
	output[[count]] = paste('>kmer_',count,sep='')
	count = count + 1
	output[[count]] = seq
}

output = output %>% as.data.frame %>% t
write.table(output,'apit_earth_v_space_kmers_sig_seqs_pos.fa',row.names=F,col.names=F,quote=F)

output = list()
count = 0
for(seq in seqs_neg){
	count = count + 1
	output[[count]] = paste('>kmer_',count,sep='')
	count = count + 1
	output[[count]] = seq
}

output = output %>% as.data.frame %>% t
write.table(output,'apit_earth_v_space_kmers_sig_seqs_neg.fa',row.names=F,col.names=F,quote=F)