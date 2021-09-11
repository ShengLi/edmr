#' Empirical differentially methylated regions
#' 
#' Comprehensive DMR analysis based on bimodal normal distribution model and 
#' weighted cost function for regional methylation analysis optimization. 
#' It captures the regional methylation modification by taking the spatial 
#' distribution of CpGs into account for the enrichment DNA methylation 
#' sequencing data so as to optimize the definition of the empirical regions. 
#' Combined with the dependent adjustment for regional p-value combination.
#' 
#' @param myDiff a \code{data.frame} object created by calculateDiffMeth from methylKit package and converted into data.frame. Required.
#' @param step a numeric variable for calculating auto-correlation, default: 100.
#' @param dist distance cutoff to call a gap for DMR, default: "none", which will be automatically determined by the bimodal normal distribution, default: 100.
#' @param DMC.qvalue qvalue cutoff for DMC definition, default: 0.01
#' @param DMC.methdiff methylation difference cutoff for DMC definition, default: 25.
#' @param num.DMCs cutoff of the number DMCs in each region to call DMR, default: 1.
#' @param num.CpGs cutoff of the number of CpGs, default: 3.
#' @param DMR.methdiff cutoff of the DMR mean methylation difference, default=20.
#' @param plot plot the bimodal normal distribution fitting or not, default=FAlSE.
#' @param main the title of the plot, if plot=TRUE. Default=FALSE.
#' @param mode the mode of call DMRs. 1: using all CpGs together. 2: use unidirectional CpGs to call DMRs. default: 1.
#' @param ACF p-value combination test with (TRUE, default) or without (FALSE) dependency adjustment. 
#' @param fuzzypval p-value cutoff for raw regions definition.
#' @importFrom data.table data.table
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges findOverlaps
#' @importFrom IRanges Rle
#' @importFrom IRanges IRanges
#' @return \code{GRanges}
#' @export
#' @examples
#' library(GenomicRanges)
#' library(IRanges)
#' library(mixtools)
#' library(data.table)
#' data(edmr)
#' mydmr=edmr(myDiff[1:5000,], mode=1, ACF=FALSE)
#' mysigdmr=filter.dmr(mydmr)
edmr=function(myDiff, step=100, dist="none", DMC.qvalue=0.01, DMC.methdiff=25, num.DMCs=1, num.CpGs=3, DMR.methdiff=20, plot=FALSE, main="", mode=1, ACF=TRUE, fuzzypval=1){
  myDiff=as.data.frame(myDiff)
  if(mode==1){
    myDMR=eDMR.sub(myDiff, step=step, dist=dist, DMC.qvalue=DMC.qvalue, DMC.methdiff=DMC.methdiff, num.DMCs=num.DMCs, num.CpGs=num.CpGs, 
                   DMR.methdiff=DMR.methdiff, granges=TRUE, plot=plot, main=main, direction="both", ACF=ACF, fuzzypval=fuzzypval)
  } else if (mode==2){
    hyper.myDMR=eDMR.sub(myDiff, step=step, dist=dist, DMC.qvalue=DMC.qvalue, DMC.methdiff=DMC.methdiff, num.DMCs=num.DMCs, num.CpGs=num.CpGs,
                         DMR.methdiff=DMR.methdiff, granges=TRUE, plot=plot, main=main, direction="hyper", ACF=ACF)
    hypo.myDMR=eDMR.sub(myDiff, step=step, dist=dist, DMC.qvalue=DMC.qvalue, DMC.methdiff=DMC.methdiff, num.DMCs=num.DMCs, num.CpGs=num.CpGs,
                        DMR.methdiff=DMR.methdiff, granges=TRUE, plot=plot, main=main, direction="hypo", ACF=ACF, fuzzypval=fuzzypval)
    myDMR=c(hyper.myDMR, hypo.myDMR)
  }
  else {
    stop ("mode = 1 or 2")
  }
}
