#' Calculate emprical raw regions from \code{myDiff} object.
getPeaks=function(allMyDiff, pcutoff=0.1, dist=100){
  print("raw myDiff:")
  print(dim(allMyDiff))
  
  myDiff=allMyDiff[allMyDiff$ppvalue<=pcutoff,]
  probes.dist=diff(myDiff$pend)
  print(paste("max cpgs dist for region definition:", dist))
  bpoints=which(probes.dist>=dist | probes.dist<0)
  first.peak=c(as.character(myDiff[1,1]), myDiff[1,2],end=myDiff[bpoints[1],3])
  mid.peaks=cbind(myDiff[(bpoints[-length(bpoints)]+1),1:2],end=myDiff[bpoints[-1],3])
  last.idx=bpoints[length(bpoints)]+1
  last.peak=c(as.character(myDiff[last.idx,1]), myDiff[last.idx,2], end=myDiff[nrow(myDiff),3])
  peaks=as.data.frame(rbind(first.peak, mid.peaks, last.peak))
  colnames(peaks)=c("rchr","rstart","rend")
  peaks
}
