#load packages
library("tidyverse")
library("edgeR")


args = commandArgs(trailingOnly=TRUE)

#Input count table
input = args[1]
if (is.na(input)) {print("Input file required as first positional argument") ; quit()} else {print(paste("Input is ", input, "", sep = ' '))}

#Metadata file
metadata = args[2]
if (is.na(metadata)) {print("Metedata file required as second positional argument") ; quit()} else {print(paste("Metadata is ", metadata, "", sep = ' '))}

#Output name
output = args[3]
if (is.na(output)) {print("Output name required as third positional argument") ; quit()} else {print(paste("Output is ", output, "", sep = ' '))}

#Name of column in the metadata file that divide the samples in the groups with biological replicates
num_groupcolumn = args[4]
if (is.na(num_groupcolumn)) {print("Group column required as fourth positional argument") ; quit()} else {print(paste("Group columns is ", num_groupcolumn, "", sep = ' '))}

#List of contrasts names
contrasts_file = args[5]
if (is.na(contrasts_file)) {print("Contrasts file required as fifth positional argument") ; quit()} else {print(paste("Contrasts file is ", contrasts_file, "", sep = ' '))}

#Threshold for significant FDR
threshold = args[6]
if (is.na(threshold)) {print("Significativity threshold required as sixth positional argument") ; quit()} else {print(paste("Significativity threshold is ", threshold, "", sep = ' '))}

print("starting to work")
#Prepare workig dataframe
data = read.csv(input, sep = "\t")
targets = read.csv(metadata, sep = "\t")
Group <- factor(targets[, num_groupcolumn])
cbind(targets,Group=Group)
#Now I estimate the dispersion
design <- model.matrix(~0+Group)
colnames(design) <- levels(Group)
contrasts_list <- scan(contrasts_file, what="", sep="\n")
datedge <- DGEList(counts=data[,-1], group=Group, genes = data[1])
datedge <- calcNormFactors(datedge)
datedge <- estimateDisp(datedge, design)
fit <- glmQLFit(datedge, design)

#Now we repeat something for each contrast

for(i in contrasts_list) {
  name = strsplit(i, " = ")[[1]][1]
  contr = strsplit(i, " = ")[[1]][2]
  contrast <- makeContrasts(contr, levels = design)
  
  #Now we look for the DEGs!
  
  qlf <- glmQLFTest(fit, contrast=contrast)
  
  #export plot
  png(paste(output, name, ".png", sep = ''))
  plotMD(qlf)
  abline(h=c(-1, 1), col="blue")
  dev.off()
  
  ##Add association
  tab <- topTags(qlf, n=Inf)$table
  tab$association <- 'Not_associated'
  contrast_elements <- strsplit(name, "VS")[[1]]
  tab <- tab %>% mutate(association = ifelse(round(FDR, digits = 3) < threshold & logFC > 0, contrast_elements[1], tab$association))
  tab <- tab %>% mutate(association = ifelse(round(FDR, digits = 3) < threshold & logFC < 0, contrast_elements[2], tab$association))
  
  ##export DEG table
  write.table(tab, file=paste(output, name, ".txt", sep = ''), row.names=FALSE, sep = "\t", quote = FALSE)
  
}
