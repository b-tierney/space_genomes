# parse and translate pyseer output

library(tidyverse)

# load mapping dictionaries 
cogmap = read.csv('cog_fx_map')
ecmap = read.csv('ec_mapping',header=F,sep='\t') %>% rename(variant=V1,name=V2)

# space vs earth
cog_output = read.csv('pyseer_EARTH_SPACE_BINARY_apit_cog_data.tsv',sep='\t') %>% mutate(VARTYPE = 'COG') %>% mutate(by = p.adjust(lrt.pvalue,method='BY'))
ec_output = read.csv('pyseer_EARTH_SPACE_BINARY_apit_ec_data.tsv',sep='\t') %>% mutate(VARTYPE = 'EC')  %>% mutate(by = p.adjust(lrt.pvalue,method='BY'))
gene_output = read.csv('pyseer_EARTH_SPACE_BINARY_apit_gene_data.tsv',sep='\t') %>% mutate(VARTYPE = 'GENE')  %>% mutate(by = p.adjust(lrt.pvalue,method='BY'))
prod_output = read.csv('pyseer_EARTH_SPACE_BINARY_apit_prod_data.tsv',sep='\t') %>% mutate(VARTYPE = 'PRODUCT') %>% mutate(by = p.adjust(lrt.pvalue,method='BY'))

cog_output_merged = left_join(cog_output,cogmap,by=c('variant'='COG'))
ec_output_merged = left_join(ec_output,ecmap,by=c('variant'))
prod_output = prod_output %>% mutate(name = variant)
gene_output = gene_output %>% mutate(name = variant) 

bound_data = bind_rows(cog_output_merged,ec_output_merged,prod_output,gene_output)
bound_data = bound_data %>% relocate(name) %>% select(-variant)

write.csv(bound_data,'earth_space_cleaned_fxnl_output.csv')

# flight data bound
cog_output = read.csv('pyseer_FLIGHT_BINARY_apit_cog_data.tsv',sep='\t') %>% mutate(VARTYPE = 'COG') %>% mutate(by = p.adjust(lrt.pvalue,method='BY'))
ec_output = read.csv('pyseer_FLIGHT_BINARY_apit_ec_data.tsv',sep='\t') %>% mutate(VARTYPE = 'EC')  %>% mutate(by = p.adjust(lrt.pvalue,method='BY'))
gene_output = read.csv('pyseer_FLIGHT_BINARY_apit_gene_data.tsv',sep='\t') %>% mutate(VARTYPE = 'GENE')  %>% mutate(by = p.adjust(lrt.pvalue,method='BY'))
prod_output = read.csv('pyseer_FLIGHT_BINARY_apit_prod_data.tsv',sep='\t') %>% mutate(VARTYPE = 'PRODUCT') %>% mutate(by = p.adjust(lrt.pvalue,method='BY'))

cog_output_merged = left_join(cog_output,cogmap,by=c('variant'='COG'))
ec_output_merged = left_join(ec_output,ecmap,by=c('variant'))
prod_output = prod_output %>% mutate(name = variant)
gene_output = gene_output %>% mutate(name = variant) 

bound_data = bind_rows(cog_output_merged,ec_output_merged,prod_output,gene_output)
bound_data = bound_data %>% relocate(name) %>% select(-variant)

write.csv(bound_data,'flight_binary_cleaned_fxnl_output.csv')