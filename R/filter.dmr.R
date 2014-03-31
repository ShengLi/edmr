#' significant DMRs
#' @export
#' @param myDMR a \code{GRanges} object for DMR regions called.
#' @param DMR.qvalue qvalue cutoff for DMC definition, default: 0.001
#' @param mean.meth.diff cutoff of the DMR mean methylation difference, default=20
#' @param num.CpGs cutoff of the number DMCs in each region to call DMR, default: 5
#' @param num.DMCs cutoff of the number of CpGs, default: 3
#' @examples
#' library(methylKit)
#' library(GenomicRanges)
#' library(mixtools)
#' library(data.table)
#' data(edmr)
#' mydmr=edmr(chr22.myDiff, mode=1, ACF=TRUE)
#' mysigdmr=filter.dmr(mydmr)

filter.dmr=function(myDMR, DMR.qvalue=0.001, mean.meth.diff=20, num.CpGs=5, num.DMCs=3){
  x=myDMR; 
  if(class(x)=="GRanges"){
    idx=which(values(myDMR)[,"DMR.qvalue"]<=DMR.qvalue & abs(values(myDMR)[,"mean.meth.diff"])>=mean.meth.diff & values(myDMR)[,"num.CpGs"]>=num.CpGs & values(myDMR)[,"num.DMCs"] >=num.DMCs)
    return(x[idx])
  } else if(x=="data.frame"){
    idx=which(myDMR$DMR.qvalue<=DMR.qvalue & abs(myDMR$mean.meth.diff) >= mean.meth.diff & num.CpGs>=num.CpGs & num.DMCs >=num.DMCs)  
  }
  return(x[idx,])
}
