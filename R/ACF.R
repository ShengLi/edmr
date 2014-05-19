
getCorr=function(lag.idx,myDiff,i){
  res=c()
  lag.idx.x=c()
  lag.idx.y=c()
  for(k in 1:i){
    idx=which(lag.idx==k)
    lag.idx.x=c(lag.idx.x,myDiff$ppvalue[idx])
    lag.idx.y=c(lag.idx.y,myDiff$ppvalue[idx+1])
    corr=cor(lag.idx.x,lag.idx.y)
    if(corr==0) corr=9e-16
    res=c(res,abs(corr))
  }
  res
}
# Auto correlation function
# @importFrom data.table data.table
ACF=function(dist,step, myDiff){
  probes.dist=myDiff$pend[-1]-myDiff$pend[-nrow(myDiff)]
  
  dist.div=probes.dist/step
  lag.idx=ceiling(dist.div)
  dt=data.table(cbind(dist=1:dist,key=ceiling(1:dist/step)))
  key <- NULL
  acf.dist=dt[,list(start=min(dist), end=max(dist)), by=key]
  acf.dist$end[nrow(acf.dist)]=acf.dist$start[nrow(acf.dist)]+step-1
#   dt=data.frame(dist=1:dist,key=ceiling(1:dist/step))
#   acf.dist=ddply(dt, .(key), summarize, start=min(dist), end=max(dist))
#   acf.dist$end[nrow(acf.dist)]=acf.dist$start[nrow(acf.dist)]+step-1
  acf=getCorr(lag.idx,myDiff,nrow(acf.dist)); 
  acfTable=cbind(V1=acf.dist$start,V2=acf.dist$end,V3=acf)
  data.frame(acfTable)
}
