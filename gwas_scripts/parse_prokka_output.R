# parse prokka output for associations

library(tidyverse)
library(reshape2)

prokka_mapping = read.csv('prokka_id_mapping',sep='\t',header=F)
colnames(prokka_mapping) = c('Sample_ID','prokkaid')
data = read.csv('apit_functional_data',sep='\t')
data = data %>% filter(locus_tag!='locus_tag',grepl('_',locus_tag))
data = data %>% mutate(prokkaid = strsplit(locus_tag,'_') %>% map_chr(1))

data = left_join(data,prokka_mapping,by='prokkaid')
data = data %>% relocate(Sample_ID) %>% select(-locus_tag,-ftype,-length_bp) 

data_cogs = data %>% filter(COG!='') %>% dcast(Sample_ID ~ COG) %>% mutate_at(vars(-Sample_ID),function(x) if_else(x>=1, 1, 0))%>% column_to_rownames('Sample_ID') %>% t %>% as.data.frame %>% rownames_to_column('COG')
write.table(data_cogs,'apit_cog_data.tsv',row.names=F,quote=F,sep='\t')
data_prods = data %>% filter(product!='') %>% dcast(Sample_ID ~ product) %>% mutate_at(vars(-Sample_ID),function(x) if_else(x>=1, 1, 0)) %>% column_to_rownames('Sample_ID') %>% t %>% as.data.frame %>% rownames_to_column('PRODUCT')
write.table(data_prods,'apit_prod_data.tsv',row.names=F,quote=F,sep='\t')
data_ec = data %>% filter(EC_number!='') %>% dcast(Sample_ID ~ EC_number) %>% mutate_at(vars(-Sample_ID),function(x) if_else(x>=1, 1, 0)) %>% column_to_rownames('Sample_ID') %>% t %>% as.data.frame %>% rownames_to_column('EC')
write.table(data_ec,'apit_ec_data.tsv',row.names=F,quote=F,sep='\t')
data_gene = data %>% filter(gene!='') %>% dcast(Sample_ID ~ gene) %>% mutate_at(vars(-Sample_ID),function(x) if_else(x>=1, 1, 0)) %>% column_to_rownames('Sample_ID') %>% t %>% as.data.frame %>% rownames_to_column('GENE')
write.table(data_gene,'apit_gene_data.tsv',row.names=F,quote=F,sep='\t')


