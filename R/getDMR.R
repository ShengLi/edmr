# Stouffer test
# @param x p-values.
stouffer=function(x)pnorm(sum(qnorm(x))/sqrt(length(x))) 

# Stoffer Liptak test
# @param pvals p-values.
# @param end myDiff \code{end} column.
# @param acf output object from ACF function.
# @param step base pairs in each step of auto-correlation calculation.
stouffer.liptak=function(pvals,end,acf=acf, step=100){
  n=length(pvals)
  if(n>1){
    comb=combn(n,2)
    a=matrix(1, ncol=n, nrow=n)
    for(cc in 1:ncol(comb)){
      i=comb[1,cc]
      j=comb[2,cc]
      dist=end[j] - end[i]
      if(dist<1) dist=1
      corr=acf$V3[ceiling(dist/step)]
      a[i,j]=a[j,i]=corr
    }
    sigma=a
    pvals[pvals>=1]=1.0-9e-16
    pnorm(sum(solve(chol(sigma, pivot=T))%*%qnorm(pvals))/sqrt(length(pvals)))   
  } else {
    pvals
  }
}
# Calculate the p-value for the differentially methylated regions (ACF=TRUE)
getDMR=function(peaks, allMyDiff, pcutoff=0.1,step=100, DMC.qvalue=0.01, DMC.methdiff=25, num.DMCs=1, num.CpGs=3, DMR.methdiff=20, a=TRUE){
#  library(GenomicRanges,quietly =TRUE)
#  library(data.table,quietly =TRUE)
  pmethdiff <- ppvalue <- pqvalue <- rchr <- rstart <- rend <- pend <- NULL
  myDiff=allMyDiff[allMyDiff$ppvalue<=pcutoff,]
  # in cases where there are many pvalues equal to zero 
  min.pval=min(9e-16,myDiff$ppvalue[myDiff$ppvalue!=0])
  min.qval=min(9e-16,myDiff$pqvalue[myDiff$pqvalue!=0])
  myDiff$ppvalue[myDiff$ppvalue==0]=min.pval
  myDiff$pqvalue[myDiff$pqvalue==0]=min.qval  
  
  myDiff.gr=GRanges(seqnames=Rle(myDiff$pchr), IRanges(start=as.integer(myDiff$pstart), end=as.integer(myDiff$pend)), strand=myDiff$pstrand, pval=myDiff$ppvalue, padj=myDiff$pqvalue, methDiff=myDiff$pmethdiff)
  peaks.gr=GRanges(seqnames=Rle(peaks$rchr), IRanges(start=as.integer(peaks$rstart), end=as.integer(peaks$rend)))
  overlap.idx=findOverlaps(myDiff.gr, peaks.gr)
  #dt.pk.myD=data.table(cbind(peaks[overlap.idx@to,], myDiff[overlap.idx@from,]))
  dt.pk.myD=data.table(cbind(peaks[overlap.idx@subjectHits,], myDiff[overlap.idx@queryHits,]))
  #print(head(dt.pk.myD))
  refine.pk.myD=dt.pk.myD[, list(medianmethdiff=median(pmethdiff),
                                meanmethdiff=mean(pmethdiff),
                                num_probes=length(ppvalue),
                                num_DMCs=length(which(pqvalue<=DMC.qvalue & abs(pmethdiff)>=DMC.methdiff))),
                         by=list(rchr,rstart,rend)]
#   refine.pk.myD=ddply(dt.pk.myD, .(rchr,rstart,rend), summarize, medianmethdiff=median(pmethdiff),
#                                                      meanmethdiff=mean(pmethdiff),
#                                                      num_probes=length(ppvalue),
#                                                      num_DMCs=length(which(pqvalue<=DMC.qvalue & abs(pmethdiff)>=DMC.methdiff)))
  refine.idx=which(refine.pk.myD$num_DMCs >= num.DMCs & refine.pk.myD$num_probes >= num.CpGs & abs(refine.pk.myD$meanmethdiff)>=DMR.methdiff)
  refine.pk.gr=refine.pk.myD[refine.idx,]
  refine.pk.idx=which(paste(dt.pk.myD$rchr,dt.pk.myD$rstart, dt.pk.myD$rend, sep="_") %in% paste(refine.pk.gr$rchr,refine.pk.gr$rstart, refine.pk.gr$rend, sep="_"))
  if(a){
    maxdist=max(dt.pk.myD[refine.pk.idx, as.integer(rend)-as.integer(rstart)])
    print(paste("auto correlation calculation. step:", step, "dist:", maxdist))
    acf=ACF(maxdist,step, allMyDiff)
#     res.pk.myD=ddply(dt.pk.myD[refine.pk.idx,], .(rchr,rstart,rend), summarize, 
#                      medianmethdiff=median(pmethdiff),
#                                              meanmethdiff=mean(pmethdiff),
#                                              num_probes=length(ppvalue),
#                                              num_DMCs=length(which(pqvalue <= DMC.qvalue & abs(pmethdiff) >= DMC.methdiff)),
#                                              methdiffs=paste(pmethdiff, collapse=","),
#                                              pvalues=paste(ppvalue,collapse=","),
#                                              r.sl.qval=stouffer.liptak(pqvalue[order(pend)], pend[order(pend)], acf),
#                                              r.sl.pval=stouffer.liptak(ppvalue[order(pend)], pend[order(pend)], acf))
    res.pk.myD=dt.pk.myD[refine.pk.idx, list(medianmethdiff=median(pmethdiff),
                                             meanmethdiff=mean(pmethdiff),
                                             num_probes=length(ppvalue),
                                             num_DMCs=length(which(pqvalue <= DMC.qvalue & abs(pmethdiff) >= DMC.methdiff)),
                                             methdiffs=paste(pmethdiff, collapse=","),
                                             pvalues=paste(ppvalue,collapse=","),
                                             r.sl.qval=stouffer.liptak(pqvalue, pend, acf),
                                             r.sl.pval=stouffer.liptak(ppvalue, pend, acf)),
                         by=list(rchr,rstart,rend)]    
  } else {
    res.pk.myD=dt.pk.myD[refine.pk.idx, list(medianmethdiff=median(pmethdiff), 
                                             meanmethdiff=mean(pmethdiff), 
                                             num_probes=length(ppvalue), 
                                             num_DMCs=length(which(pqvalue <= DMC.qvalue & abs(pmethdiff) >= DMC.methdiff)),
                                             methdiffs=paste(pmethdiff, collapse=","),
                                             pvalues=paste(ppvalue,collapse=","),
                                             r.sl.qval=stouffer(pqvalue),
                                             r.sl.pval=stouffer(ppvalue)),
                         by=list(rchr,rstart,rend)]
  }
  r.sl.padj=p.adjust(res.pk.myD$r.sl.pval,'fdr')
  r.sl.qadj=p.adjust(res.pk.myD$r.sl.qval,'fdr')
  DMR=cbind(res.pk.myD, r.sl.qadj, r.sl.padj)
  res=as.data.frame(DMR)
  colnames(res)=c("chr","start","end","median.meth.diff","mean.meth.diff","num.CpGs","num.DMCs","meth.diffs","pvalues","DMR.q.pvalue","DMR.p.pvalue","DMR.q.qvalue","DMR.p.qvalue")
  res
}

