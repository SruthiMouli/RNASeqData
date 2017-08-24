source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2")
source('analyzeDE.R')
install.packages('gplots')
library("RColorBrewer")
##Read csv
data_xls= read.csv("Lund_Neb2test_annotated_counts.csv",sep='\t')
## Read files from excel


######################################################

############################################################
data_mat = data_xls[,c(1,2:4,6:8)]
write.csv(data_mat,sep='\t',file = "Data_for_AnalyzeDE.csv")
#Format in DESeq readable matrix
##Select required data and read to matrix
DESeq_res = analyzeDE(fileName = 'Data_for_AnalyzeDE.csv',c("D25minus","D25minus","D25minus","D25plus","D25plus","D25plus"),"D25minus","D25plus",0.05,1)
##Get the last 6 columns
res = DESeq_res[,(ncol(DESeq_res)-6+1):ncol(DESeq_res)]
hmcol = colorRampPalette(brewer.pal(11, "RdYlGn"))(100)
res=as.matrix(res)
## Heat map for all six values
heatmap.2(res,Colv=F,col=hmcol,scale='row',dendrogram='row',key=T,density.info='none',trace='none')
## Take row means of two conditions
Exp_Mean= cbind(rowMeans(x[,c(1:3)]),rowMeans(x[,c(4:6)]))
heatmap.2(Exp_Mean,Colv=F,col=hmcol,scale='row',dendrogram='row',key=T,density.info='none',trace='none')