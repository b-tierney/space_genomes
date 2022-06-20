
library(tidyverse)

files = list.files()
files = files[grep('tsv',files)]

d= read.table(files[[1]],header=F)%>% mutate(V3 = V3/V2) %>% select(V1,V3)
d_nn =read.table(files[[1]],header=F) %>% select(V1,V3)
colnames(d)[2]=gsub('representative_viral_sequences_','',gsub('.tsv','',files[[1]]))
colnames(d_nn)[2]=gsub('representative_viral_sequences_','',gsub('.tsv','',files[[1]]))

for(f in files[2:length(files)]){
	d2=read.table(f,header=F) %>% mutate(V3 = V3/V2) %>% select(V1,V3)
#	d2$V3 = d2$V3/sum(d2$V3)
	d2_nn=read.table(f,header=F) %>% select(V1,V3)
	colnames(d2)[2]=gsub('representative_viral_sequences_','',gsub('.tsv','',f))
	colnames(d2_nn)[2]=gsub('representative_viral_sequences_','',gsub('.tsv','',f))
	d=left_join(d,d2)
	d_nn=left_join(d_nn,d2_nn)
}

d = d %>% column_to_rownames("V1")
viralcounts =  colSums(d) %>% data.frame %>% rownames_to_column('V1')
colnames(viralcounts)[2] = 'counts'
d=t(d) %>% as.data.frame 

d_nn = d_nn %>% column_to_rownames("V1")
d_nn=t(d_nn) %>% as.data.frame 
write.table(d_nn,'assembled_viral_abundances_viral_counts.csv',sep=',')

d = d %>% rownames_to_column('V1')

counts=read.table('sample_read_counts',header=F)
colnames(counts) = c('total_reads','V1')
counts$total_reads = as.numeric(counts$total_reads)
counts$V1 = gsub('GLDS_252_metagenomics_','GLDS-252_metagenomics_',counts$V1)

out = left_join(d,counts)
out = out %>% column_to_rownames("V1")
out = out/out$total_reads
out = out %>% select(-total_reads)
write.table(out,'assembled_viral_abundances_total_reads.csv',sep=',')

out = left_join(d,viralcounts)
out = out %>% column_to_rownames("V1")
out = out/out$counts
out = out %>% select(-counts)
write.table(out,'assembled_viral_abundances_viral_reads.csv',sep=',')