#' calculate differentially methylated regions from \code{myDiff} object.
myDiffToDMR=function(myDiff, dist=100, step=100, DMC.qvalue=0.01, DMC.methdiff=25, num.DMCs=1, num.CpGs=3, DMR.methdiff=20, ACF=TRUE, fuzzypval=0.1){
  if (class(myDiff)=="methylDiff") input=data.frame(myDiff)
  else if(class(myDiff)=="data.frame") input=myDiff
  else print("Input object myDiff has too be methylDiff class or data.frame class")
  # prepare myDiff for DMR calling
  idx=with(input,order(chr, start))
  #myDiff.new=cbind(input[idx,2], input[idx,3]-1, input[idx,4:8])
  myDiff.new=cbind(input[idx,"chr"], input[idx,"start"]-1, input[idx,c("end","strand","pvalue", "qvalue", "meth.diff")])
  colnames(myDiff.new)=c("pchr","pstart","pend","pstrand","ppvalue","pqvalue","pmethdiff")
  
  peaks=getPeaks(myDiff.new, pcutoff=fuzzypval, dist=dist)
  #if(ACF==TRUE){
    res=getDMR(peaks, myDiff.new, pcutoff=fuzzypval, step=step, DMC.qvalue=DMC.qvalue, DMC.methdiff=DMC.methdiff, num.DMCs=num.DMCs, num.CpGs=num.CpGs, DMR.methdiff=DMR.methdiff, a=ACF)  
  #} else if (ACF==FALSE){
  #  res=getDMR2(peaks, myDiff.new, pcutoff=1, step=step, DMC.qvalue=DMC.qvalue, DMC.methdiff=DMC.methdiff, num.DMCs=num.DMCs, num.CpGs=num.CpGs, DMR.methdiff=DMR.methdiff)  
  #}
  return(res)
}
