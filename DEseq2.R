library(limma)
library(edgeR)
library(dplyr)
library(RColorBrewer)
library(gplots)
library(Rtsne)

df = fread("/Users/tariqhicham/Documents/EPFL/Master/genomics and bioinformatics/Project/brca_tcga_pub2015 (1)/data_RNA_Seq_v2_expression_median.txt")
df <- df[rowSums(df[,3:819])>0,] #remove genes that are not expressed in all the samples
names(df)






DGE <- DGEList(df[,3:819])


#add cancer type from metadata

metadf <- fread("/Users/tariqhicham/Documents/EPFL/Master/genomics and bioinformatics/Project/brca_tcga_pub2015 (1)/data_clinical_sample.txt",skip =4)
metadf <- metadf[order(metadf$SAMPLE_ID),] 


#there is an issu in the metadf : TCGA-BH-A1ES patient has 2 sample_id and only one correspond to one gene_count (TCGA-BH-A1ES-06)
#I am going to remove the element  TCGA-BH-A1ES-01 that is not in df
metadf= subset(metadf, metadf$SAMPLE_ID != 'TCGA-BH-A1ES-01')




type = as.factor(metadf$ONCOTREE_CODE)

DGE$samples$type <- type

#ILC_samples = metadf[metadf$ONCOTREE_CODE == "ILC"]$SAMPLE_ID
#ILC_sample2 <- metadf[col]$ONCOTREE_CODE

ILC_samples = metadf[type == "ILC"]$SAMPLE_ID


#Add the information about the genes in the DGE object in a section genes : ENTREZ ID,GID, gene length

#extract entrez Id from df
ENTREZ_ID = df$Entrez_Gene_Id


GID <- df$Hugo_Symbol


#Extract the genes lengths from 'gene_length.txt'
gene_length <- fread('gene_length.txt')


Genes = data.frame(ENTREZ_ID,GID,gene_length$lengths)
DGE$genes <- Genes


#extract the ILC patients 
ILC = DGE[,DGE$samples$type=='ILC']#,keep.lib.size = FALSE]

#Attribute to each ILC patient its subtype 
subtype_assignment <- fread('subtype_assignment.txt')
#Sort the subtype asignment 
subtype_assignment<-subtype_assignment[order(subtype_assignment$`Sample ID`),]
subtype = as.factor(subtype_assignment$`Class Assignment`)
ILC$samples$subtype <- subtype



#Data pre-processing 
#here I used RPKM instead of TPM, as there is a function in edgeR
#But I could directly used my TPM.txt instead of df
#I will also take the log2 of the RPKM


#filtering
#Genes that are not expressed at all are already removed

#keep.exprs <- filterByExpr(ILC)#,group = ILC$samples$type) #not sure if this was the right group + 14298 genes
#ILC <- ILC[keep.exprs,, keep.lib.sizes=FALSE]



#test

L <- mean(ILC$samples$lib.size) * 1e-6
M <- median(ILC$samples$lib.size) * 1e-6
lcpm.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(ILC)
col <- brewer.pal(nsamples, "Paired")
lcpm <- cpm(ILC,log = TRUE)

par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-CPM")
abline(v=lcpm.cutoff, lty=3)

for (i in 2:10){#only check the density of the tenth first patient
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
#gene filtering : remove genes that are lowly expressed
keep.exprs <- filterByExpr(ILC,group = ILC$samples$subtype)#,group = ILC$samples$type) #not sure if this was the right group + 14298 genes
ILC <- ILC[keep.exprs,, keep.lib.sizes=FALSE]
lcpm =cpm(ILC,log =TRUE)


plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:10){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}






#Normalizing the data to prevent any batch effect
#Normalization method is trimmed mean of M-values (TMM).




par(mfrow =c(1,2))
#boxplot non-normalized
boxplot(lcpm[,1:5],names =NULL)


#Normalization
ILC<- calcNormFactors(ILC, method = "TMM")
ILC$samples$norm.factors

#boxplot normalized
lcpm =cpm(ILC,log=TRUE)
boxplot(lcpm[,1:5],names=NULL) #normalization does not really change the boxplot because the batch effect 
#was not very strong (factor1)


col.group <- ILC$samples$subtype
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")

par(mfrow=(c(1,1)))
plotMDS(lcpm, col=col.group,pch=1)
title(main="Dimensionality reduction of ILC types")
legend("bottomright",levels(ILC$samples$subtype),col = levels(col.group),pch=1)





ILC$samples$type <- factor(ILC$samples$type)
design <- model.matrix(~0+ILC$samples$subtype)
colnames(design) <- c('ImmuneRelated',"Proliferative","ReactiveLike")

design

contr.matrix <- makeContrasts(
  ImmuneRelated -Proliferative ,
  ImmuneRelated - ReactiveLike, 
  Proliferative - ReactiveLike,
  levels = colnames(design))
contr.matrix


par(mfrow=c(1,2))
v <- voom(ILC, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit)

summary(decideTests(efit))
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)


par(mfrow=c(2,2))
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], xlim=c(-8,13))
plotMD(tfit, column=2, status=dt[,2], main=colnames(tfit)[2], xlim=c(-8,13))
plotMD(tfit, column=3, status=dt[,3], main=colnames(tfit)[3], xlim=c(-8,13))


Immune_vs_proliferative = topTreat(tfit,coef = 1,n=Inf,p.value = 0.05)
top100.I_P = Immune_vs_proliferative[1:100,]
Immune_vs_Reactive = topTreat(tfit,coef = 2,n=Inf,p.value = 0.05)
top100.I_R = Immune_vs_Reactive[1:100,]
Proliferative_vs_Reactive = topTreat(tfit,coef = 3,n=Inf,p.value = 0.05)
top100.I_R = Immune_vs_Reactive[1:100,]


Immune_vs_proliferative_topgenes <- Immune_vs_proliferative$ENTREZ_ID[1:20]
Immune_vs_Reactive_topgenes <- Immune_vs_Reactive$ENTREZ_ID[1:20]
Proliferative_vs_Reactive_topgenes_ID <- Proliferative_vs_Reactive$ENTREZ_ID[1:20]


topgenes <- c(Immune_vs_proliferative_topgenes,Immune_vs_Reactive_topgenes,Proliferative_vs_Reactive_topgenes_ID)
i <- which((v$genes$ENTREZ_ID %in% Immune_vs_proliferative_topgenes) | (v$genes$ENTREZ_ID%in% Immune_vs_Reactive_topgenes) | (v$genes$ENTREZ_ID %in% Proliferative_vs_Reactive_topgenes_ID))




#with 60 genes classifier
ilc_subtype_genes <- fread("ILC 60 Gene Classifier.txt")

mycol <- brewer.pal(3, "Dark2")
f <- factor(ILC$samples$subtype)

j= which(v$genes$ENTREZ_ID %in%ilc_subtype_genes$UID)

lcpm <- cpm(ILC,log =TRUE)
heatmap.2(lcpm[j,], scale="row",
          labRow=v$genes$GID[j],labCol =FALSE ,
          ColSideColors = mycol[f],
          symbreaks=FALSE, key=TRUE, symkey=FALSE,
          col=rev(redgreen(250)), margin=c(5,18), trace="none", density.info="none",cexRow = 0.6,
          lhei=c(1,5), dendrogram="row",hclustfun = function(x) hclust(x, method='ward.D'))
legend('right',legend = c("Reactive-Like","Immune-Related","Proliferative"), fill =c('#7570B3',"#1B9E77","#D95F02"),
       cex=.7,title ="ILC subtypes")




#with DE genes found




lcpm <- cpm(ILC,log =TRUE)
heatmap.2(lcpm[i,], scale="row",
          labRow=v$genes$GID[i],labCol =FALSE ,
          ColSideColors = mycol[f],
          symbreaks=FALSE, key=TRUE, symkey=FALSE,
          col=rev(redgreen(250)), margin=c(5,18), trace="none", density.info="none",cexRow = 0.6,
          lhei=c(1,5), dendrogram="row",hclustfun = function(x) hclust(x, method='ward.D2'))
legend('right',legend = c("Reactive-Like","Immune-Related","Proliferative"), fill =c('#7570B3',"#1B9E77","#D95F02"),
       cex=.7,title ="ILC subtypes")

###########################
# T-SNE for ILC subtypes  #
##########################

colors = rainbow(length(unique(ILC$samples$subtype)))
names(colors) = unique(ILC$samples$subtype)

tSNE_subtype=Rtsne(t(lcpm[i,]), dims = 2,perplexity = 10, verbose=TRUE)
plot(tSNE_subtype$Y,col=colors[ILC$samples$subtype],main ="t-SNE projection for ILC subtypes", xlab ="t-SNE dimension 1",ylab = "t-SNE dimension 2")
legend("topright",c('Reactive-like',"Immune-related","Proliferative"),col = colors,pch=1,cex=1)











