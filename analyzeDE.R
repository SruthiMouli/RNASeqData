#' Analyze Sequencing counts for differential expression
#'
#' This function allows you to analyze Sequencing counts using DESeq package.
#' @param fileName Seq counts csv file
#' @param condVec Vector of conditions in the same order as in the seq counts file
#' @param controlName Name of the control sample
#' @param caseName Name of the case sample
#' @param adjPVal Adjusted p-value to be used for cut-off
#' @param logFC Log fold-change values to be used for cut-off
#' @keywords DESeq Sequencing
#' @return This function will analyze Seq counts data for differential expression.
#' @export
#' @import BiocGenerics Biobase DESeq locfit
#' @importFrom utils read.csv write.csv
#' @examples
#' analyzeDE('Example_SeqFile.txt',c('Case1','Case1','Case1','Control','Control','Control'),'Control','Case1',0.05,1)
analyzeDE <- function(fileName,condVec,controlName,caseName,adjPVal,logFC) {
    a <- utils::read.csv(fileName,header = T)
    b <- a[,2:ncol(a)]
    if(is.vector(condVec) == FALSE){print('Conditions should be in a vector format. Pls use c()')}
    if(length(condVec) != ncol(b)){print('Conditions vector should have the same length as the conditions in the Seq data')}
    row.names(b) <- a[,1]
    countTable <- as.matrix(b[rowSums(b) != 0, ])
    storage.mode(countTable) = 'integer'
    cds <- DESeq::newCountDataSet(countTable, factor(condVec))
    cds <- DESeq::estimateSizeFactors(cds)
    cds <- DESeq::estimateDispersions(cds, method = 'pooled', sharingMode = 'maximum', fitType = 'local')
    vsd <- DESeq::getVarianceStabilizedData(cds)
    mod_lfc <- (rowMeans(vsd[,conditions(cds) == caseName, drop = FALSE]) - rowMeans(vsd[,conditions(cds) == controlName, drop = FALSE]))
    res <- DESeq::nbinomTest(cds, controlName, caseName)
    res <- cbind(res, mod_lfc)
    res2 <- res[which(res$padj < adjPVal),]
    res3 <- res2[which(abs(res2$mod_lfc) > logFC),]
    vsd1 <- vsd[res3$id,conditions(cds) %in% c(controlName,caseName)]
    res4 <- cbind(res3,vsd1)
    write.csv(res4,paste0('DE_',caseName,'-',controlName,'.csv'),row.names = FALSE)
    return(res4)}