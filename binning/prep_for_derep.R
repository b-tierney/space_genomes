# for each sample

library(tidyverse)

#binning output made by ls */*assembled/scaffolds.fasta_binning/ | grep : | cut -f1 -d: > binning_output
binning_output_locs = read.csv('binning_output',header=F) %>% unlist %>% unname

for(b in binning_output_locs){
	statsfiles = list.files(b)[list.files(b) %>% grep('stats',.)]
	if(length(statsfiles)!=0){
		binning_stats = list()
		for(f in statsfiles){
			bintype = strsplit(f,'\\.') %>% map_chr(1)
			binning_stats[[f]] = read.csv(paste(b,f,sep=''),sep='\t') %>% mutate(path=paste(b,bintype,'/',bin,'.fa',sep=''),basename = strsplit(b,'/') %>% map_chr(2) %>% gsub('_assembled','',.))
		}
		binning_stats = bind_rows(binning_stats)
		binfold = paste(b,'bins_non_derep',sep='')
		binning_stats_filt = binning_stats %>% filter(completeness > 90,contamination<5)
		if(nrow(binning_stats_filt)!=0){
			system(paste('mkdir',binfold))
			for(p in unlist(unname(binning_stats_filt %>% select(path)))){
				binname = gsub('/','-',p)
				system(paste("cp ",p," ",binfold,"/",binname,sep=''))
			}
			derepoutdir = paste(b,"DEREPPED",sep='')
			system(paste("dRep dereplicate ",derepoutdir," -g ",binfold,"/*fa",sep=''))
		}
	}
}

