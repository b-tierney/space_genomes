# parse kmer association output

library(tidyverse)

data = read.table('apit_earth_v_space_kmers.txt',sep='\t',header=T)
data = data %>% mutate(by = p.adjust(lrt-pvalue,method='BY')) %>% filter(by <0.05)

data %>% select(XXX) %>% unlist %>% unname

output = list()
count = 0
for(seq in data){
	count = count + 1
	output[[count]] = paste('>kmer_',count,sep='')
	count = count + 1
	output[[count]] = seq
}

output = output %>% as.data.frame
write.table(output,row.names=F,header=F)