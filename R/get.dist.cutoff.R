# inner function 1
# get break point to minimize the cost function.
get.break_point=function(mixmdl, itv=range(mixmdl$mu)){
  if(mixmdl$mu[1] < mixmdl$mu[2]) {
    i=1;k=2
  } else{i=2; k=1}
  break_point = optimize(
    f = function(x){
      p1=pnorm(x, mixmdl$mu[i],mixmdl$sigma[i], lower.tail = F)
      p2=pnorm(x, mixmdl$mu[k],mixmdl$sigma[k], lower.tail = T)
      cost=sum(mixmdl$lambda[i]*p1,mixmdl$lambda[k]*p2)
      return(cost)
    }
    , interval = itv
    , maximum=F
  )$minimum
  break_point
}
# get distance cutoff from \code{mixtools} object
get.dist.cutoff <- function(mixmdl){
  dist=round(2^get.break_point(mixmdl))
  dist
}
