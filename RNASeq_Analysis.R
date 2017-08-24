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
DESeq_res=analyzeDE(fileName = 'Data_for_AnalyzeDE.csv',c("D25minus","D25minus","D25minus","D25plus","D25plus","D25plus"),"D25minus","D25plus",0.05,1)

## Replace row.names from ensembl Ids to external names
result_Data = merge(DESeq_res,data_xls,by = "names.dds.")
result_Data = result_Data[order(result_Data$log2FoldChange,decreasing = TRUE),]
vsd_exprs= cbind(result_Data$vsd_ctrl,result_Data$vsd_case)
row.names(vsd_exprs) = result_Data$external_gene_name
colnames(vsd_exprs) = c("vsd_ctrl","vsd_case")
heatmap.2(vsd_exprs, trace="none", margin=c(10, 6),scale='row')
