library(data.table)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(dplyr)
library(gplots)# contains the heatmap.2 package
library(pheatmap)
library(RColorBrewer)


#load the data
#df <- fread("/Users/polina/Documents/EPFL/TA_GB_2020/projects/4_cancer/main/cell2015/brca_tcga_pub2015/data_RNA_Seq_v2_expression_median.txt")
#df <- fread("/Users/polina/Documents/EPFL/TA_GB_2020/projects/4_cancer/main/cell2015/brca_tcga_pub2015/data_RNA_Seq_v2_mRNA_median_Zscores.txt")
df = fread("/Users/tariqhicham/Documents/EPFL/Master/genomics and bioinformatics/Project/brca_tcga_pub2015 (1)/data_RNA_Seq_v2_expression_median.txt")

# How many genes? How many samples?
nrow(df)
# 20440 genes studied (all the genes in the human genome, as it is obtained by whole exome sequencing), but some of them are not expressed in any sample, such as SSX9
# --> Need to filter those genes, in order to enhance the time processing

ncol(df[, grep("TCGA", names(df)), with = FALSE] )
# there are 817 samples, the two first columns display the names of the genes and their ids 



# Preprocess data
#remove genes that are not expressed in any of the samples
df <- df[rowSums(df[,3:819])>0,]
#295 genes were filtered out

#Get gene lengths
#install.packages("rentrez")
library(rentrez)

# make "for" loop for all genes to get gene length

genes <- entrez_summary(db="gene", id=391343 ) # id is from df$Entrez_Gene_Id        100134869
genes_length <- abs(genes$genomicinfo$chrstart - genes$genomicinfo$chrstop)  # this is rough estimate of gene length

# id 8225 is located in the chromosome Y and X -->it is counted twice 
#id 391343 is weird --> no length

#Get all the ids of the genes expressed by the samples
genes_ids <-df$Entrez_Gene_Id


  length <- function(id){
  genes <- entrez_summary(db="gene",id)
  genes_lengths <-abs(genes$genomicinfo$chrstart - genes$genomicinfo$chrstop)
  return (genes_lengths)
}



genes_lengths = sapply(genes_ids,length)

#it took 3h to get the length so do not remove it 
#Some genes are expressed in allosomes and have 2 lengths ==> we need to remove one of them


L <- genes_lengths
for(i in which(lengths(L)==2)) {
  L[[i]] <- L[[i]][1]
}

#Other have no size because there is no genomic info in the esummary. This is 
for(k in which(lengths(L)==0) ){
  genes <- entrez_summary(db ="gene",id=genes_ids[k])
  L[[k]] <- abs(strtoi(genes$locationhist$chrstart[1]) - strtoi(genes$locationhist$chrstop[1]) )  #for gene 9500, chrstop was a str so I had to change it to convert it to int
}

#very very long


# 3 lengths are still empty : PRINS (id:100169750, 13945 in df), GAGE4 (id:2576, 6762 in df) and GAGE8 (id: 100101629, 6763 in df).
#PRINS : long non coding RNA. Based on genbank, 2199bp
#GAGE4 : is apparently an alias of GAGE 1 (id: 2543) --≥ 528bp
#GAGE8 : GenBank: AF055473.1 -->528 bp

L[c(6762,6763,13945)]=c(528,528,2199) 

Q<-data.frame(Hugo_symbol = df$Hugo_Symbol)

Q$lengths <- unlist(L)

write.table(Q,'gene_length.txt',row.names = F)

Q<- read.table('gene_length.txt',header = TRUE)



