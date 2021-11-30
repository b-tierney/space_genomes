#!/usr/bin/env Rscript

# convert mmseqs output to binary sparse matrix

args = commandArgs(trailingOnly=TRUE)

library(tidyverse)

clusterfile = args[[1]]
prokka_id_mapping = args[[2]]

clusterfile = 'clusters_90perc_cluster.tsv'
prokka_id_mapping = 'prokka_id_mapping'

data = read.csv(clusterfile,sep='\t',header=F)
colnames(data) = c('congene','rawgene')
prokkamap = read.csv(prokka_id_mapping,sep = '\t',header=F)
colnames(prokkamap) = c('sample_id','prokka_id')

data = data %>% mutate(prokka_id = strsplit(rawgene,'_') %>% map_chr(1))
data = left_join(data,prokkamap,by='prokka_id')
data = data %>% select(-prokka_id,-rawgene) %>% relocate(sample_id) %>% mutate(val=1)

data_wide = dcast(data, sample_id ~ congene, value.var="val")
genecols = colnames(data_wide %>% select(-sample_id))

# convert to binary if...that's what you want...? idk man
data_wide[data_wide>0] = 1

# merge in metadata if it's there
mdat = read.csv('pyseer_metadata.tsv',sep='\t')

data_wide = left_join(data_wide,mdat,by = c('sample_id'='ASSEMBLY'))

data %>% select([whateverthemdatcolumnis],all_of(genecols))