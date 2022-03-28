

setwd('C:/Users/mariav/Asaf_analysis/analysis_10f/scripts')

f_counts_ECsPerTaxa='../user_input/Route_4_input_genus_BjSa_NTC_new.txt'
f_counts_ECsPerTaxa_design='../user_input/Mazzola_metadata.txt'
ECsPerTaxa_EC="EC"
ECsPerTaxa_taxa="taxa"
ECsPerTaxa_design_sample="sample"
ECsPerTaxa_design_group="group"
out1='../output/route4'

dominance_cutoff=0.4
select_n_taxa=200

##############################

library(data.table)
library(vegan)
library(ggplot2)

# load data:

t_counts_ECsPerTaxa_1 = read.csv(f_counts_ECsPerTaxa,sep="\t",stringsAsFactors=F)
t_counts_ECsPerTaxa_design_1 = read.csv(f_counts_ECsPerTaxa_design,sep="\t",stringsAsFactors=F)

# check columns definitions:
procede=T
if(!(ECsPerTaxa_design_sample %in% names(t_counts_ECsPerTaxa_design_1))){cat("undefined sample column in design\n");procede=F}
if(!(ECsPerTaxa_design_group %in% names(t_counts_ECsPerTaxa_design_1))){cat("undefined group column in design\n");procede=F}

samples=t_counts_ECsPerTaxa_design_1[,ECsPerTaxa_design_sample]
group0=t_counts_ECsPerTaxa_design_1[,ECsPerTaxa_design_group]
groups=unique(group0)

if(sum(samples %in% names(t_counts_ECsPerTaxa_1))!=length(samples)){cat("undefined EC column\n");procede=F}
if(!(ECsPerTaxa_EC %in% names(t_counts_ECsPerTaxa_1))){cat("undefined EC column\n");procede=F}
if(!(ECsPerTaxa_taxa %in% names(t_counts_ECsPerTaxa_1))){cat("undefined taxa column\n");procede=F}

if(!procede){stop("stopped: wrong columns definitions\n")}

# data statisitcs:

t_counts_ECsPerTaxa_2=t_counts_ECsPerTaxa_1[,c(ECsPerTaxa_EC,ECsPerTaxa_taxa,samples)]
names(t_counts_ECsPerTaxa_2)[1]='EC'
names(t_counts_ECsPerTaxa_2)[2]='taxa'
stat1=list()
for(s in 1:length(samples)){
  cat(samples[s],"\n")
  abundance1 = data.table(t_counts_ECsPerTaxa_2[,c('EC','taxa',samples[s])])
  names(abundance1)=c('EC','taxa','counts')
  stat1[[s]] = abundance1[,.(sum=sum(counts),
                             shannon=diversity(counts,index='shannon'),
                             simpson=diversity(counts,index='simpson'),
                             dominance=1-diversity(counts,index='simpson'),
                             top_taxa=taxa[max(counts)==counts][1],
                             top_taxa_count=counts[max(counts)==counts][1]),
                          by = .(EC)]
  gg1=ggplot(data.frame(stat1[[s]]),aes(x=dominance,y=shannon))+geom_point()+xlab('dominance (1-simpson)')+ylab('diversity (shannon)')
  pdf(paste0(out1,'.',samples[s],'.dominance_vs_diversity.pdf'))
  print(gg1)
  dev.off()
  temp=data.frame(stat1[[s]])
  names(temp)[2:ncol(temp)]=sapply(names(temp)[2:ncol(temp)],function(x)paste0(samples[s],'.',x))
  if(s==1){ stat2=temp
  }else if(sum(stat2$EC==temp$EC)==nrow(stat2)){
    stat2=cbind(stat2,temp[,2:ncol(temp)])
  }else{cat("error\n")}
}

# ECs selection:

dominance1=list()
topTaxa1=list()
ecs1=stat1[[1]]$EC
#for(g in 1:length(groups)){#4 groups
g=2
  cat(groups[g],"\n")
  group_name=gsub('\\s+|\\.','_',groups[g],perl=T)
  # write all ECs:
  sink(file=paste0(out1,'.',group_name,'.ECs_with_dominant_taxon.txt'),append=F)
  cat(paste0("all ",paste(ecs1,collapse=" "),"\n"))
  sink()
  sink(file=paste0(out1,'.',group_name,'.ECs_with_dominant_taxon_knockout.txt'),append=F)
  cat(paste0("all ",paste(ecs1,collapse=" "),"\n"))
  sink()
  # select experiment-specific ECs:
  samples1=samples[group0 %in% groups[g]]
  idxs1=(1:length(samples))[group0 %in% groups[g]]
  dominance1[[g]]=data.frame(EC=ecs1)
  topTaxa1[[g]]=data.frame(EC=ecs1)
  for(s in idxs1){
    if(sum(stat1[[s]]$EC == ecs1) == length(ecs1)){
      dominance1[[g]][,samples[s]] = stat1[[s]][,dominance] >= dominance_cutoff
      topTaxa1[[g]][,samples[s]] = stat1[[s]][,top_taxa]
    }else{cat("error\n")}
  }
  count_dom=apply(dominance1[[g]][,2:ncol(dominance1[[g]])],1,sum)
  taxa_dom = topTaxa1[[g]][count_dom > 0,]
  taxa_dom$most_freq= apply(taxa_dom[,2:ncol(taxa_dom)],1,function(x) ifelse(length(x)>1,names(sort(table(x),decreasing=T))[1],x[1]))
  taxa_dom_freq = data.frame(sort(table(taxa_dom$most_freq),decreasing=T))
  taxa_dom_freq$taxa=rownames(taxa_dom_freq)
  #taxa_dom_freq=taxa_dom_freq[,c(2,1)]
  names(taxa_dom_freq)=c('taxa','cases_dominant')
  ECs_with_dominant_taxa=list()
  # write experiment-specific ECs:
  if(nrow(taxa_dom_freq) > 0 ){
    n=1:min(select_n_taxa,nrow(taxa_dom_freq))
    for(j in n){
      taxa1=taxa_dom_freq$taxa[j]
      taxa1_name=gsub('\\s+|\\.','_',taxa1,perl=T)
      ECs_with_dominant_taxa[[j]]=taxa_dom[taxa_dom$most_freq == taxa1,'EC']
      ecs_knockout = ecs1[!(ecs1 %in% ECs_with_dominant_taxa[[j]])]
      sink(file=paste0(out1,'.',group_name,'.ECs_with_dominant_taxon.txt'),append=T)
      cat(paste0(group_name,"_",taxa1_name," ",paste(ECs_with_dominant_taxa[[j]],collapse=" "),ifelse((g==length(groups))&(j==length(n)),"","\n")))
      sink()
      sink(file=paste0(out1,'.',group_name,'.ECs_with_dominant_taxon_knockout.txt'),append=T)
      cat(paste0(group_name,"_",taxa1_name," ",paste(ecs_knockout,collapse=" "),ifelse((g==length(groups))&(j==length(n)),"","\n")))
      sink()
    }
  }else{cat("no taxa was found to be dominant\n")}

