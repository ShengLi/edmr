#' edmr sub function
#' @importFrom GenomicRanges GRanges
eDMR.sub=function(myDiff, step=100, dist="none", DMC.qvalue=0.01, DMC.methdiff=25, num.DMCs=1, num.CpGs=3, DMR.methdiff=20, granges=TRUE, plot=FALSE, main="", direction="both", ACF=TRUE, fuzzypval=0.1){
  if(direction=="both"){
    print("DMR analysis for all detected CpGs...")
  } else if (direction=="hyper") {
    print("DMR analysis for hyper methylated CpGs...")
    idx=which(myDiff$meth.diff>0)
    myDiff=myDiff[idx,]
  } else if (direction=="hypo") {
    print("DMR analysis for hypo methylated CpGs...")
    idx=which(myDiff$meth.diff<0)
    myDiff=myDiff[idx,]
  } else {
    print ("parameter direction has too be both, hyper or hyper")
  }
  if(dist=="none"){
    mixmdl=myDiff.to.mixmdl(myDiff, plot=plot, main=main)
    dist=get.dist.cutoff(mixmdl)    
  }
  DMR=myDiffToDMR(myDiff, dist=dist, step=step, DMC.qvalue=DMC.qvalue, DMC.methdiff=DMC.methdiff, num.DMCs=num.DMCs, num.CpGs=num.CpGs, DMR.methdiff=DMR.methdiff, ACF=ACF, fuzzypval=fuzzypval)
  myDMR=DMR[,c(1:3,5:7,10,12)]
  if(granges) {
    myDMR.gr=with(myDMR, GRanges(chr, IRanges(as.numeric(start), as.numeric(end)), mean.meth.diff=mean.meth.diff, num.CpGs=num.CpGs, num.DMCs=num.DMCs, DMR.pvalue=DMR.q.pvalue, DMR.qvalue=DMR.q.qvalue))
    return(myDMR.gr)
  } else {
    return(myDMR)
  }
}
