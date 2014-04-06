# inner function 1
get.dist.myDiff <- function(myDiff){
  x=myDiff
  chr.list=unique(as.character(x$chr))
  chr.list=chr.list[chr.list!="chrM"]
  #y=foreach(chr = chr.list,.combine=c) %dopar% {diff(sort(x$start[x$chr==chr]))}
  y=c()
  for(chr in chr.list) {y=c(y, diff(sort(x$start[x$chr==chr])) )}
  y
}

#' Distance to mixmdl
#' @importFrom mixtools normalmixEM
dist_to_mixmdl <- function(dist)
{
  log2.distance=log2(dist[dist!=1])
  mixmdl=normalmixEM(log2.distance)
}

#' obtain mixtools model from \code{myDiff} object
#' @export
#' @param myDiff a \code{data.frame} object created by calculateDiffMeth from methylKit package and converted into data.frame. Required.
#' @param plot to plot or not the nearest distance distribution.
#' @param main title of the plot
#' @examples
#' library(edmr)
#' library(mixtools)
#' data(edmr)
#' 
#' myMixmdl=myDiff.to.mixmdl(myDiff)
myDiff.to.mixmdl=function(myDiff,plot=F, main=""){
  dist=get.dist.myDiff(as.data.frame(myDiff))
  mixmdl=dist_to_mixmdl(dist)
  if(plot){
    plotMdl1(mixmdl,main)
  }
  print(2^get.break_point(mixmdl))
  mixmdl
}
