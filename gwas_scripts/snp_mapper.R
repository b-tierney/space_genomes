# snapper

library(tidyverse)
library(rtracklayer)
library(CompleHeatmap)
# creating positional heatmaps from snippy output

setwd('~/Dropbox (Mason Lab)/space_genomes/ap_subproject/revisions/snp_calling/')

# load in snippy core output tab file

snp_id = read.table('typestrain.tab',header=T,sep='\t')

# load in file containing gene/fnxl mapping information

annotation_data = readGFF('typestrain.gff', version=0,columns=NULL, tags=NULL, filter=NULL, nrows=-1,raw_data=FALSE) %>% as.data.frame %>% dplyr::rename(CHR=seqid)

# load metadata 

metadata = read.csv('../../metadata_20211015.csv') %>% select(ASSEMBLY,HUMAN_FLIGHT_OTHER) %>% column_to_rownames('ASSEMBLY')
metadata_anno = metadata %>% as.data.frame %>% rownames_to_column('ASSEMBLY')
mapping = read.table('assembly_biosample_mapping')
colnames(mapping) = c('ASSEMBLY','sample')
metadata_anno = inner_join(metadata_anno,mapping) %>% filter(grepl('MT',HUMAN_FLIGHT_OTHER))

# join and keep only snps in genes
merged = dplyr::left_join(snp_id, annotation_data, by = c("CHR" = "CHR")) %>% filter(POS < end & POS > start)  %>% select(all_of(colnames(snp_id)),start,end,gene,product,locus_tag)
freq = merged$locus_tag %>% table %>% data.frame
colnames(freq) = c('locus_tag','number_of_loci')
merged = left_join(merged,freq,by='locus_tag') %>% mutate(length = end-start)
merged$length_norm = merged$length/median(merged$length)

# plot gene product by number of mutations
forfreqplot = merged %>% select(product,number_of_loci) %>% distinct %>% filter(number_of_loci > 100)  %>% arrange(desc(number_of_loci))
forfreqplot$product = as.factor(forfreqplot$product)
forfreqplot$product = fct_reorder(forfreqplot$product,forfreqplot$number_of_loci)
forfreqplot$temp = seq(1,nrow(forfreqplot))
ggplot(forfreqplot,aes(x = temp, y= number_of_loci)) + theme_bw() + geom_bar(stat='identity')  + scale_x_continuous(breaks = forfreqplot$temp,labels=forfreqplot$product)+ theme(axis.text.x = element_text(angle = 60,hjust=1)) + ylab('Number of SNPs') + xlab('')
ggsave('typestrain_freq.pdf',width=12,height=6)

# generate sample by sample plots for genes with logs of snps
loci = merged %>% filter(number_of_loci > 100,product != 'hypothetical protein') %>% select(locus_tag,number_of_loci) %>% distinct %>% arrange(desc(number_of_loci)) %>% head(10)%>% select(locus_tag)  %>% unique %>% unlist %>% unname 
ht_list = NULL

for(l in loci){
  sub = merged %>% filter(locus_tag == l) %>% select(all_of(colnames(snp_id)),length,product,start,end) %>% select(-CHR)
  prod = unique(sub$product)
  start = unique(sub$start)
  end = unique(sub$end)
  length = sub$length_norm %>% unique
  sizelist[[l]] = length
  sub = sub %>% select(-length)
  temp = data.frame(seq(start,end)) 
  colnames(temp) = c('POS')
  sub_exp = left_join(temp,sub) %>% mutate(product = prod)
  sub_exp$REF[is.na(sub_exp$REF)] = 'NO SNP FOUND'
  sub_exp[is.na(sub_exp)] = 'NO SNP FOUND'
  sub_exp = sub_exp %>% column_to_rownames('POS')
  sub_anno_col = sub_exp %>% select(product,REF)
  sub_exp_hm = sub_exp %>% select(-REF,-product) %>% t %>% as.data.frame%>% rownames_to_column('sample')
  sub_exp_hm$sample = gsub('_type.strain_2','',sub_exp_hm$sample)
  sub_exp_hm = inner_join(sub_exp_hm,metadata_anno) %>% select(-sample) %>% column_to_rownames('ASSEMBLY') %>% arrange(desc(HUMAN_FLIGHT_OTHER))
  sub_anno_row = sub_exp_hm %>% select(HUMAN_FLIGHT_OTHER)
  sub_exp_hm = sub_exp_hm %>% select(-HUMAN_FLIGHT_OTHER)
  colors =  structure(1:5, names = c('NO SNP FOUND','A','T','G','C'))
  colors2 =  structure(6:7, names = c('MT1','MT2'))
  row_anno = rowAnnotation(FLIGHT = sub_anno_row$HUMAN_FLIGHT_OTHER,col=list(FLIGHT = colors2))
  sub_anno_col = sub_anno_col  %>% mutate(label = (rownames(sub_anno_col)))
  sub_anno_col = sub_anno_col %>% mutate(label = if_else(as.numeric(label)%%500 == 0,label,''))
  anno_col = columnAnnotation(BASE = sub_anno_col$REF,col=list(BASE = colors))
  anno_col = columnAnnotation(BASE = sub_anno_col$REF,col=list(BASE = colors))
  h=Heatmap(sub_exp_hm,cluster_columns=F,cluster_rows=F,width = length,show_heatmap_legend = F,show_row_names = F,column_title_rot = 270,show_column_names = T,col=colors,border = T,left_annotation = row_anno,column_labels = sub_anno_col$label,column_title = prod,bottom_annotation = anno_col)
  pdf(paste(l,'typestrain.pdf',sep=''),width=8,height=3)
  print(Heatmap(sub_exp_hm,cluster_columns=F,cluster_rows=F,width = length,show_heatmap_legend = F,show_row_names = F,column_title_rot = 0,height=.1,show_column_names = T,col=colors,border = T,left_annotation = row_anno,column_labels = sub_anno_col$label,column_title = prod,bottom_annotation = anno_col))
  dev.off()
  ht_list = ht_list + h
}
