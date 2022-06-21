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
library(broom)
library(ggtree)
library(scales)
library(tidyverse)
library(ggrepel)
library(ggnewscale)
library(reshape2)
library(RColorBrewer)
library(pheatmap)
library(cowplot)

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

#regression_output = map(colnames(data_sub_t) ,function(x) glm(data = data_sub_t_annos, family=binomial(),earth_space ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10) %>% tidy %>% filter(term != '(Intercept)',!grepl('PC',term)) %>% mutate(term = x)) %>% bind_rows()
#regression_output = regression_output %>% mutate(BY = p.adjust(p.value,method='BY'))

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

### snp analysis

setwd('~/Dropbox (Mason Lab)/space_genomes/ap_subproject/revisions/snp_calling/')

ts_dist = read.table('typestrain_dists.tab',sep='\t',header=T) %>% filter(Reference=='Reference') %>% mutate(Reference = 'Type Strain')
colnames(ts_dist) = c('r1','r2','dist')
ts_dist$r1 = gsub('_type-strain_2','',ts_dist$r1)
ts_dist$r2 = gsub('_type-strain_2','',ts_dist$r2)

bonomo1_dist = read.table('bonomo1_dists.tab',sep='\t',header=T) %>% filter(Reference=='Reference') %>% mutate(Reference = 'earth_strain_1')
colnames(bonomo1_dist) = c('r1','r2','dist')
bonomo1_dist$r1 = gsub('_bonomo1_2','',bonomo1_dist$r1)
bonomo1_dist$r2 = gsub('_bonomo1_2','',bonomo1_dist$r2)

bonomo2_dist = read.table('bonomo2_dists.tab',sep='\t',header=T)%>% filter(Reference=='Reference') %>% mutate(Reference = 'earth_strain_2')
colnames(bonomo2_dist)  = c('r1','r2','dist')
bonomo2_dist$r1 = gsub('_bonomo2_2','',bonomo2_dist$r1)
bonomo2_dist$r2 = gsub('_bonomo2_2','',bonomo2_dist$r2)

# load and parse metadata

mapping = read.table('assembly_biosample_mapping')
colnames(mapping) = c('ASSEMBLY','sample')

ts_dist = left_join(ts_dist,mapping,by=c('r2'='sample')) 
bonomo1_dist = left_join(bonomo1_dist,mapping,by=c('r2'='sample')) 
bonomo2_dist = left_join(bonomo2_dist,mapping,by=c('r2'='sample')) 

ts_dist = left_join(ts_dist,metadata_anno,by=c('ASSEMBLY'='sample')) %>% filter(grepl('MT',HUMAN_FLIGHT_OTHER)) %>% select(HUMAN_FLIGHT_OTHER,r1,dist)
bonomo1_dist = left_join(bonomo1_dist,metadata_anno,by=c('ASSEMBLY'='sample')) %>% filter(grepl('MT',HUMAN_FLIGHT_OTHER))%>% select(HUMAN_FLIGHT_OTHER,r1,dist)
bonomo2_dist = left_join(bonomo2_dist,metadata_anno,by=c('ASSEMBLY'='sample')) %>% filter(grepl('MT',HUMAN_FLIGHT_OTHER)) %>% select(HUMAN_FLIGHT_OTHER,r1,dist)

snp_counts = bind_rows(ts_dist,bonomo1_dist,bonomo2_dist)

ggplot(data = snp_counts,aes(x = HUMAN_FLIGHT_OTHER, y = dist)) + theme_bw() + geom_boxplot() + facet_wrap(facets=vars(r1),scales = 'free') + xlab('') +ylab('SNP COUNT')
ggsave('snp_distances.pdf',width=6,height=3)

# roary tree

setwd('~/Dropbox (Mason Lab)/space_genomes/ap_subproject/revisions/tree_full_dataset/')

tree <- read.tree("raxml_run1/RAxML_bestTree.roary")
metadata = read.table('../../metadata_20211015.csv',header=T,sep=',') %>% mutate(Assembly = strsplit(ASSEMBLY,'\\.') %>% map_chr(1))
rownames(metadata)=metadata$Assembly
metadata = metadata %>% select(-Assembly,ISOLATION_SOURCE_CLEANED,HUMAN_FLIGHT_OTER)

metadata$ISOLATION_SOURCE_CLEANED = as.factor(metadata$ISOLATION_SOURCE_CLEANED)
metadata$HUMAN_FLIGHT_OTHER = as.factor(metadata$HUMAN_FLIGHT_OTHER)

new_tip_labels = list()
for(i in seq(tree$tip.label)){
  val = tree$tip.label[[i]]
  new_tip_labels[[i]] = as.character(metadata[val,'ISOLATION_SOURCE_CLEANED'])
}

new_tip_labels=unlist(unname(new_tip_labels))

labelmap = data.frame(label = tree$tip.label, newlabel = new_tip_labels)

p = ggtree(tree)+ theme_tree2() + scale_x_ggtree()

p2=gheatmap(p + geom_tiplab(align=T,label='none'), metadata %>% select(HUMAN_FLIGHT_OTHER), offset=0, width=.1, colnames=FALSE, legend_title="Isolation Source")+ scale_fill_viridis_d(option="D", name="General isolation Source")  

ggsave('apit_tree_313_raxml.pdf',height=15,width=15)

p3 = p2 %<+% labelmap + geom_tiplab(offset = 0,aes(label=newlabel),size = 2, parse = F)
tips = labelmap %>% filter(grepl('MT',newlabel)) %>% select(label) %>% unlist %>% unname

viewClade(p2, MRCA(p2, .node1=c(tips,'GCA_001415895')))%<+% labelmap+ geom_tiplab(offset = -.03,aes(label=newlabel),size = 3, parse = F,align = T)
ggsave('zoomed_clade.pdf')

p = ggtree(tree,branch.length = 'none',layout='circular') 

p2=gheatmap(p,metadata %>% select(HUMAN_FLIGHT_OTHER), offset=0, width=.1, colnames=FALSE, legend_title="Isolation Source")+ scale_fill_viridis_d(option="D", name="General isolation Source")  

p3 = p2 %<+% labelmap + geom_tiplab(offset = 6,aes(label=newlabel),size = 2, parse = F)

ggsave('apit_tree_313_raxml_nobranchlength.pdf',height=10,width=10)

# SNP tree

setwd('~/Dropbox (Mason Lab)/space_genomes/ap_subproject/revisions/snp_tree//')

tree = read.tree('clean.core.tree')
metadata = read.table('../../metadata_20211015.csv',header=T,sep=',') %>% mutate(Assembly = strsplit(ASSEMBLY,'\\.') %>% map_chr(1))
rownames(metadata)=metadata$Assembly
metadata = metadata %>% select(-Assembly,ISOLATION_SOURCE_CLEANED,HUMAN_FLIGHT_OTHER)

metadata$ISOLATION_SOURCE_CLEANED = as.factor(metadata$ISOLATION_SOURCE_CLEANED)
metadata$HUMAN_FLIGHT_OTHER = as.factor(metadata$HUMAN_FLIGHT_OTHER)

mapping = read.table('../snp_calling/assembly_biosample_mapping')
colnames(mapping) = c('ASSEMBLY','sample')

new_tip_labels = list()
for(i in seq(tree$tip.label)){
  i = tree$tip.label[[i]]
  if(i == 'Reference'){
    new_tip_labels[[i]] = 'REFERENCE'
  }
  if(i != 'Reference'){
  i2 = gsub('_type-strain_2','',i)
  mapped = mapping$ASSEMBLY[mapping$sample==i2]
  new_tip_labels[[i]] = as.character(metadata[mapped,'ISOLATION_SOURCE_CLEANED'])
  }
}

new_tip_labels=unlist(unname(new_tip_labels))

labelmap = data.frame(label = tree$tip.label, newlabel = new_tip_labels)

p = ggtree(tree) + theme_tree2() + scale_x_ggtree() 

p2 = p %<+% labelmap + geom_tiplab(offset = 0,aes(label=newlabel),size = 5, parse = F)

ggsave('apit_gubbins_snptree.pdf',height=6,width=6)
