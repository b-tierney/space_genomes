# parse kmer association output

library(tidyverse)

data = read.table('apit_earth_v_space_kmers.txt',sep='\t',header=T)
data_sub = data %>% mutate(by = p.adjust(lrt.pvalue,method='BY')) %>% filter(by <0.05)
saveRDS(data_sub,'apit_earth_v_space_kmers_sig.rds')

seqs_pos = data_sub %>% filter(beta>=0) %>% select(variant) %>% unlist %>% unname
seqs_neg = data_sub %>% filter(beta<=0)%>% select(variant) %>% unlist %>% unname

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