# revisions analysis apitt

### NEW ANALYSIS TO INCLUDE

# ROARY
  # CHI SQUARED ANALYSIS + REPORT SIG HITS SOMEHOW, NEW UNIREF ALIGNMENTS, PROBABLY JUST PLOT NON-ANNOTATED
# ANI FOR SPACE CLADE
# SNPS FOR SPADE CLADE
# MANY NEW TREE TYPES
# CIRCULAR ALIGNMENT SHOWING SNPS, GWAS OUTPUT

library(ComplexHeatmap)
library(tidyverse)
squish_trans <- function(from, to, factor) {
  
  trans <- function(x) {
    
    if (any(is.na(x))) return(x)
    
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    
    return(x)
  }
  
  inv <- function(x) {
    
    if (any(is.na(x))) return(x)
    
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    
    return(x)
  }
  
  # return the transformation
  return(trans_new("squished", trans, inv))
}

setwd('~/Dropbox (Mason Lab)/space_genomes/ap_subproject/revisions/roary_out_1655231395/')

data = read.table('gene_presence_absence.Rtab',header=T)
data = data %>% column_to_rownames('Gene')

tokeep = rowSums(data) %>% as.data.frame %>% filter(. != 313) %>% filter(.>17) %>% rownames

data_sub = data[tokeep,]
data_sub_t = t(data_sub)

metadata = read.csv('../../metadata_20211015.csv') %>% select(ASSEMBLY,HUMAN_FLIGHT_OTHER) %>% column_to_rownames('ASSEMBLY')

# save for pyseer
pyseerdata = data_sub %>% as.data.frame %>% rownames_to_column("Gene")
pyseermetadata = metadata %>% mutate(earth_space = if_else(grepl('MT',HUMAN_FLIGHT_OTHER),1,0)) %>% select(-HUMAN_FLIGHT_OTHER) %>% rownames_to_column('sample')

write.table(pyseerdata,'pyseer_genedata.tsv',quote=F,row.names = F,sep='\t')
write.table(pyseermetadata,'pyseer_metadata.tsv',quote=F,row.names = F,sep='\t')

metadata_sorted = metadata[order(rownames(metadata)), , drop=F]
data_sub_t_sorted = data_sub_t[order(rownames(data_sub_t)),]

rowanno = rowAnnotation(Isolation_Source = metadata_sorted$HUMAN_FLIGHT_OTHER)

# compute fraction of genomes of a given type a gene was found in
space_ids = metadata %>% filter(grepl('MT',HUMAN_FLIGHT_OTHER)) %>% rownames
earth_ids = metadata %>% filter(!grepl('MT',HUMAN_FLIGHT_OTHER)) %>% rownames

data_sub_t_space = data_sub_t[space_ids,]
data_sub_t_space = colSums(data_sub_t_space) %>% data.frame 
data_sub_t_space = data_sub_t_space/22
data_sub_t_space = data_sub_t_space %>% rownames_to_column('gene')
colnames(data_sub_t_space) = c('gene','space_fraction')

data_sub_t_earth = data_sub_t[earth_ids,]
data_sub_t_earth = colSums(data_sub_t_earth) %>% data.frame 
data_sub_t_earth = data_sub_t_earth/291
data_sub_t_earth = data_sub_t_earth %>% rownames_to_column('gene')
colnames(data_sub_t_earth) = c('gene','earth_fraction')

merged = inner_join(data_sub_t_space,data_sub_t_earth,by='gene') %>% column_to_rownames('gene')
merged = merged[colnames(data_sub_t_sorted), , drop=F]

colanno = columnAnnotation(space_fraction = merged$space_fraction, earth_fraction = merged$earth_fraction)

pdf('roary_presence_absence_heatmap.pdf',width=18,height=9)
Heatmap(data_sub_t_sorted,show_column_names = F,show_row_names = F,right_annotation = rowanno,bottom_annotation = colanno)
dev.off()

# do a logistic regression and report the results

data_sub_t_annos = data_sub_t %>% as.data.frame %>% rownames_to_column('sample')
metadata_anno = metadata %>% as.data.frame %>% rownames_to_column('sample')

data_sub_t_annos = inner_join(metadata_anno,data_sub_t_annos,by='sample')
 
pcdata = read.table('a_pit_distances.tsv',sep='\t',header=T,row.names=1) %>% select(-c('JC.1293_assembled','JC.1294_assembled','JC.1295_assembled')) 
pcs = prcomp(pcdata,scale=TRUE)$x
pcs = pcs %>% as.data.frame %>% rownames_to_column('sample')

data_sub_t_annos = inner_join(data_sub_t_annos,pcs,by='sample')
data_sub_t_annos = data_sub_t_annos %>% mutate(earth_space = if_else(grepl('MT',HUMAN_FLIGHT_OTHER),1,0)) %>% select(-sample,-HUMAN_FLIGHT_OTHER)

regression_output = map(colnames(data_sub_t) ,function(x) glm(data = data_sub_t_annos, family=binomial(),earth_space ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10) %>% tidy %>% filter(term != '(Intercept)',!grepl('PC',term)) %>% mutate(term = x)) %>% bind_rows()
regression_output = regression_output %>% mutate(BY = p.adjust(p.value,method='BY'))

### volcano plot

### CHECK NO PCS vs 5 PCS vs 10 PCS TO DETERMINE IF POP STRUCTURE HAS IMPACT

# upset plot showing overlap in significant findings

# upset plot showing overlap in significant findings with high prevalence 

# load in prokka -- gene name mapping for more interpretable annotations
prokkamap_prod = read.csv('~/Dropbox (Mason Lab)/space_genomes/ap_subproject/gene_catalogs/gene_product_mappping.tsv',sep='\t',header=F)
colnames(prokkamap_prod) = c('term','product')

prokkamap_cog = read.csv('~/Dropbox (Mason Lab)/space_genomes/ap_subproject/gene_catalogs/gene_cog_mapping.tsv',sep='\t',header=F)
colnames(prokkamap_cog) = c('term','COG')

regression_output = read.csv('pyseer_output_5pcs.tsv',sep='\t') %>% rename(term=variant) %>% mutate(ADJ = p.adjust(lrt.pvalue,method='bonferroni'))
regression_output = regression_output %>% mutate(odds_ratio = exp(beta))
regression_output = regression_output %>% mutate(log_odds_ratio = log10(odds_ratio))

regression_output = left_join(regression_output,merged %>% rownames_to_column('term'))

regression_output = left_join(regression_output,prokkamap_prod)
regression_output = left_join(regression_output,prokkamap_cog)
regression_output = regression_output %>% mutate(term = if_else(is.na(product),term,product))

ggplot(data = regression_output,aes(x = beta,color=space_fraction,y = -log10(ADJ))) + geom_point(size=2) + theme_bw() +geom_hline(yintercept = -log(0.05,10)) + scale_color_viridis_c(option = "viridis") + geom_label_repel(data = regression_output %>% filter(!grepl('group',term)) %>% filter(ADJ<0.05) %>% arrange(ADJ) %>% head(50),aes(label = term),color='black',box.padding   = 0.1, point.padding = 0.1,alpha=.8,max.overlaps=25,size=4,segment.color = 'grey50') + scale_x_continuous(trans = squish_trans(50,200,50),limits = c(-50,230),labels =c(-50,0,25,50,100,150,200)) 
ggsave('volcano_5pcs.pdf',width=16,height=16)






