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
library(ggbeeswarm)
library(tidyverse)
library(ggrepel)
library(ggnewscale)
library(reshape2)
library(RColorBrewer)
library(pheatmap)
library(cowplot)
library(tidyverse)
#library(rtracklayer)
library(ComplexHeatmap)
library(Biostrings)

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

# make spandx variable

pyseersnipmetadata = read.table('~/Dropbox (Mason Lab)/space_genomes/ap_subproject/revisions/snp_calling_spandx/spandx_samples_of_interest')
colnames(pyseersnipmetadata) = c('sample','sample2')
pyseersnipmetadata = left_join(pyseermetadata,pyseersnipmetadata) %>% select(-sample) %>% dplyr::rename(sample=sample2) %>% filter(!is.na(sample)) %>% select(sample,earth_space)

write.table(pyseersnipmetadata,'~/Dropbox (Mason Lab)/space_genomes/ap_subproject/revisions/snp_calling_spandx/pyseer_snpmetadata.tsv',quote=F,row.names = F,sep='\t')
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

regression_output = read.csv('pyseer_output_5pcs.tsv',sep='\t') %>% dplyr::rename(term=variant) %>% mutate(ADJ = p.adjust(lrt.pvalue,method='BY'))
regression_output = regression_output %>% mutate(odds_ratio = exp(beta))
regression_output = regression_output %>% mutate(log_odds_ratio = log10(odds_ratio))

regression_output = left_join(regression_output,merged %>% rownames_to_column('term'))

regression_output = left_join(regression_output,prokkamap_prod)
#regression_output = left_join(regression_output,prokkamap_cog)
#regression_output = regression_output %>% mutate(term = if_else(is.na(product),term,product))

p = ggplot(data = regression_output,aes(x = beta,color=space_fraction,y = -log10(ADJ))) + theme_bw() + geom_point(size=4) +geom_hline(yintercept = -log(0.05,10)) + scale_color_viridis_c(option = "viridis") + geom_label_repel(data=regression_output %>% filter(!grepl('group',term)) %>% filter(ADJ<0.05) %>% arrange(ADJ) %>% filter(!duplicated(product)) %>% head(50),aes(label=product),size=8,alpha=.7,min.segment.length = .1,direction='both',color='black') + scale_x_continuous(trans = squish_trans(50,200,50),limits = c(-50,230),labels =c(-50,0,25,50,100,150,200))  + theme(axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.x = element_text(size = 30),axis.title.y = element_text(size = 30))
ggsave(plot = p, 'volcano_5pcs.pdf',width=30,height=30)

# build out gene significant location map
regression_output = read.csv('pyseer_output_5pcs.tsv',sep='\t') %>% dplyr::rename(term=variant) %>% mutate(ADJ = p.adjust(lrt.pvalue,method='BY'))
regression_output = left_join(regression_output,merged %>% rownames_to_column('term'))

mappinginfo = read.csv('gene_presence_absence_space.csv',row.names=1)
gffstoload = mappinginfo %>% select(-Gene) %>% colnames %>% unlist %>% unname
for(g in gffstoload){
  gff = readGFF(paste('space_gffs/',g,'.gff',sep=''), version=0,columns=NULL, tags=NULL, filter=NULL, nrows=-1,raw_data=FALSE) %>% as.data.frame %>% dplyr::rename(CHR=seqid) %>% select(CHR,locus_tag,start,end,gene,product)
  mappinginfo_sub = mappinginfo %>% select(Gene,all_of(g)) %>% dplyr::rename(locus_tag=g,term=Gene)
  temp = left_join(gff,mappinginfo_sub)
  temp = left_join(temp,regression_output %>% select(term,ADJ,beta,space_fraction,earth_fraction)) %>% mutate(significant = if_else(ADJ<0.05,"Significant","0"))
  temp$significant[temp$significant=='0']='Not significant'
  temp$significant[temp$product == 'hypothetical protein']='Significant hypothetical protein'
  forplotsub = temp %>% arrange(CHR,start) %>% mutate(label =if_else(significant == 'Significant',product,''),axis = if_else(label!= '',as.character(start),'')) %>% select(CHR,beta,start,end,axis,significant,label,ADJ,gene,product,space_fraction,earth_fraction)
  contigs = forplotsub %>% select(CHR,end) %>% group_by(CHR) %>% slice_max(end,n=1) %>% dplyr::rename(length=end)
  contigs = contigs%>% ungroup %>% mutate(start_position = cumsum(length)) %>% select(-length)
  contigs = contigs %>% mutate(start_position =start_position - min(start_position))
  forplotsub = left_join(forplotsub,contigs) %>% filter(!is.na(significant))
  forplotsub$`Genomic Coordinate` =   forplotsub$start +  forplotsub$start_position 
  ggplot(forplotsub,aes(x=`Genomic Coordinate`,color=beta,y=-log10(ADJ)))+ geom_point() +geom_label_repel(data = forplotsub %>% filter(!grepl('hypothetical',product)) %>% filter(ADJ<0.05,space_fraction == 1,earth_fraction<.25) %>% arrange(ADJ) ,aes(label = gene),color='black',box.padding   = 0.1, point.padding = .1,alpha=.8,size=5,segment.color = 'grey50',label.padding = .1,min.segment.length = .1) + theme(axis.text.x = element_blank()) +theme_bw()+ geom_hline(yintercept = -log10(0.05)) + scale_color_viridis_c(option = 'plasma')
  ggsave(paste('space_gffs/',g,'.pdf',sep=''),width=20,height=5)
}

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
metadata = metadata %>% select(-Assembly,ISOLATION_SOURCE_CLEANED,HUMAN_FLIGHT_OTHER)

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

# get nodes for snp tree

snptreenodes= offspring(bar,MRCA(p2, .node1=c(tips,'GCA_001415895'))) %>% filter(!is.na(label)) %>% select(label) %>% filter(!(label %in% space_ids))
write.table(as.data.frame(snptreenodes),'snp_tree_nodes',quote=F,sep='\t',row.names=F)

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

####mGWAS OUTPUT
setwd('~/Dropbox (Mason Lab)/space_genomes/ap_subproject/revisions/snp_calling_spandx/')

# get reference contigs and lengths

contig_lengths = read.table('contig_lengths',sep='\t',header=F)
colnames(contig_lengths) = c('CHR','length')
contig_lengths = contig_lengths %>% dplyr::mutate(contig_end = cumsum(length),contig_start = contig_end -length +1 ) 
contig_lengths = contig_lengths %>%  group_by(CHR) %>% mutate(contig_start = contig_end-length)

# load and parse gff
gff = readGFF('bonomo1.gff', version=0,columns=NULL, tags=NULL, filter=NULL, nrows=-1,raw_data=FALSE) %>% as.data.frame %>% dplyr::rename(CHR=seqid) %>% select(CHR,locus_tag,start,end,gene,product)

gff = left_join(gff,contig_lengths)
gff$gene_start =   gff$start +  gff$contig_start 
gff = gff %>% mutate(gene_length = end - start, gene_end = gene_start + gene_length) %>% dplyr::rename(internal_gene_start = start,internal_gene_end = end)

# load and parse reg output
regression_output = read.csv('mgwas_assoc',sep='\t') %>% dplyr::rename(term=variant) %>% mutate(ADJ = (p.adjust(lrt.pvalue,method='BY')),ADJ_l10 = -log10(p.adjust(lrt.pvalue,method='BY')))

regression_output = regression_output %>% mutate(CHR = strsplit(term,'_') %>% map_chr(1),contig_loc = strsplit(term,'_') %>% map_chr(2),reference =  strsplit(term,'_') %>% map_chr(3), variant = strsplit(term,'_') %>% map_chr(4), type =if_else(str_length(reference)>1 | str_length(variant) > 1,'indel','snp'))

regression_output = left_join(regression_output ,contig_lengths %>% select(CHR,contig_start,contig_end,length)) %>% mutate(POS = as.numeric(contig_loc) + contig_start)

merged = dplyr::left_join(gff %>% select(-contig_start,-contig_end,-length), regression_output, by = c("CHR" = "CHR")) %>% filter(POS < gene_end & POS > gene_start)  %>% dplyr::rename(contig_length = length)
# find integenic snps
regression_output_int = regression_output %>% filter(!(term %in% merged$term)) 

regression_output_with_coords = bind_rows(merged,regression_output_int) %>% mutate(Intergenic = if_else(is.na(locus_tag),'Intergenic','In CDS'))

regression_output_with_coords = regression_output_with_coords %>% mutate(type = if_else(type == 'indel' & str_length(variant)<str_length(reference),'Deletion',if_else(type == 'indel' & str_length(reference)<str_length(variant),'Insertion',type))) %>% mutate(type = if_else(str_length(reference)==str_length(variant),'SNP',type))

write.csv(regression_output_with_coords,'mgwas_output.csv',quote=F)

# manhattan plot
ggplot(data = regression_output_with_coords,aes(x=POS,y=-log10(ADJ),color=type)) + facet_grid(cols = vars(CHR))+ geom_quasirandom(alpha=.6,size=4)  + theme(axis.text.x = element_blank()) +theme_bw()+ geom_hline(yintercept = -log10(0.05)) + xlab('Position') + theme(axis.text.x= element_blank(),legend.position = 'bottom',legend.title = element_blank(),strip.text = element_blank(),strip.background = element_blank()) + scale_color_manual(values=brewer.pal(n=2,'Set1'))
ggsave('manhattan_snp.pdf',width=18,height=5)

# number significant indels vs snps 

si = regression_output_with_coords%>% mutate(Significant = if_else(ADJ<0.05,'Yes','No')) %>% select(type,Significant) %>% table %>% data.frame
colnames(si) = c('type','Significant','freq')
si$type = toupper(si$type)
ggplot(si,aes(x=type,y=freq,fill=Significant,group=Significant)) + geom_bar(stat='identity',position='dodge') + xlab('') + ylab('Number of features') + theme_bw() + theme(axis.text.x = element_text(hjust =1,angle = 60)) + geom_text(aes(label = freq,group = Significant), vjust = 0,position = position_dodge(.9)) + scale_fill_manual(values=c('grey','black')) + theme(legend.position = 'bottom')
ggsave('indels.pdf',width=2,height=4)

# number intergenic vs genic

ig = regression_output_with_coords %>% mutate(Significant = if_else(ADJ<0.05,'Yes','No'))%>% select(Intergenic,Significant) %>% table %>% data.frame
colnames(ig) = c('type','Significant','freq')
ig$type = toupper(ig$type)
ggplot(ig,aes(x=type,y=freq,fill=Significant,group=Significant)) + geom_bar(stat='identity',position='dodge') + xlab('') + ylab('Number of features') + theme_bw() + theme(axis.text.x = element_text(hjust =1,angle = 60)) + geom_text(aes(label = freq,group = Significant), vjust = 0,position = position_dodge(.9)) + scale_fill_manual(values=c('grey','black')) + theme(legend.position = 'bottom')
ggsave('intergenic.pdf',width=3,height=4)

# merged

both = regression_output_with_coords %>% mutate(Significant = if_else(ADJ<0.05,'Yes','No'))%>% select(type,Intergenic,Significant) %>% table %>% data.frame
colnames(both) = c('Type','Location','Significant','freq')
#both$type = toupper(both$type)
ggplot(both,aes(x=Location,y=freq,fill=Significant,group=Significant),color='black') + geom_bar(stat='identity',position='dodge') + facet_wrap(facets = vars(Type)) + xlab('') + ylab('Number of features') + theme_bw() + theme(axis.text.x = element_text(hjust =1,angle = 45)) + geom_text(aes(label = freq,group = Significant), vjust = 0,position = position_dodge(.9)) + scale_fill_manual(values=c('grey','black')) + theme(legend.position = 'bottom') 
ggsave('merged_snpdata.pdf',width=9,height=6)

# find genes with the most significant mutations
freqplot = regression_output_with_coords %>% filter(Intergenic == 'In CDS')%>% filter(ADJ<0.05,beta>0) %>% select(CHR,internal_gene_start,internal_gene_end,locus_tag,product,gene_start,gene_end) %>% mutate(gene_pos = paste(comma(internal_gene_start),comma(internal_gene_end),sep='-'))  %>% group_by(CHR,locus_tag,product,gene_pos,gene_start) %>% count %>% ungroup %>% arrange(desc(n)) %>% mutate(label = paste(product,' (',gene_pos,')',sep='')) %>% arrange(gene_start) %>% filter(n>1)
freqplot$label = fct_reorder(freqplot$label,freqplot$gene_start)
freqplot$CHR = as.factor(freqplot$CHR)

ggplot(freqplot,aes(x = label, y= n)) + theme_bw() +facet_grid(cols=vars(CHR),scales='free')+ geom_bar(stat='identity',position='dodge') + ylab('Number of significant mutations') + xlab('')  +  theme(axis.text.x = element_text(hjust =1,angle = 45))
ggsave('mutation_freqplot.pdf',width = 9,height=6)

# hypothetical protein = MNFGEKDL_02574, let's blast it and check it out

# get the information for the 4 gene sequences of interest specifically and plot

fa = readDNAStringSet("../snp_calling_snippy/bonomo1.ffn")
names(fa) = strsplit(names(fa),' ') %>% map_chr(1)

lociofint = unique(freqplot$locus_tag)

outputseqs=list()
for(l in lociofint){
  seq = fa[[l]]
  seq = strsplit(as.character(fa[[l]]),'')[[1]]
  sub = regression_output_with_coords %>% filter(locus_tag == l) %>% select(variant,type,gene_length,product,gene_start,gene_end,POS)
  prod = unique(sub$product)
  start = unique(sub$gene_start)
  end = unique(sub$gene_end)
  length = sub$gene_length %>% unique
  sub = sub %>% select(-gene_length)
  temp = data.frame(seq(start,end),seq) 
  colnames(temp) = c('POS','REFSEQ')
  sub_exp = left_join(temp,sub,by='POS') %>% mutate(product = prod) %>% mutate(variant = if_else(is.na(variant),REFSEQ,variant), type = if_else(is.na(type),'None',type))
  sub_exp$type = factor(sub_exp$type,levels= c('None','SNP','Insertion','Deletion'))
  sub_exp = sub_exp  %>% mutate(variation = if_else(type =='None','Reference','Variant'))
  sub_exp$variation = factor(sub_exp$variation,levels= c('Reference','Variant'))
  ggplot(sub_exp,aes(x=POS,y=variation,color=type)) + geom_quasirandom(groupOnX = F,width = .07) + theme_bw() + scale_color_manual(values = c('gray','blue','red','darkgreen')) + xlab('Open Reading Frame position') + ylab('') + ggtitle(prod) + theme(axis.text.x = element_blank())
  ggsave(paste(l,'variant_map.pdf',sep=''),width=8,height=2)
  variantseq=gsub(', ','',toString(sub_exp$variant))
  refseq=gsub(', ','',toString(sub_exp$REFSEQ))
  outputseqs[[paste(l,'ref')]] = paste(sep='','>REF_',l,' ',prod)
  outputseqs[[paste(l,'refseq')]] = refseq
  outputseqs[[paste(l,'var')]] = paste(sep='','>VAR_',l,'',prod)
  outputseqs[[paste(l,'varseq')]] = variantseq
}

writeLines(paste(outputseqs,sep='\n'),'seqs_snps_references_bonomo1.fa')




