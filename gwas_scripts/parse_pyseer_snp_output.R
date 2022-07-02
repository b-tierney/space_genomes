# parse pyseer output for SNPs

library(tidyverse)

setwd('~/Dropbox (Mason Lab)/space_genomes/ap_subproject/revisions/snp_calling_spandx/')

regression_output = read.csv('mgwas_assoc',sep='\t') %>% dplyr::rename(term=variant) %>% mutate(ADJ = -log10(p.adjust(lrt.pvalue,method='bonferroni')))

regression_output = left_join(regression_output,merged %>% rownames_to_column('term'))

regression_output = left_join(regression_output,prokkamap_prod)
regression_output = left_join(regression_output,prokkamap_cog)