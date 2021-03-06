\name{myDiff.to.mixmdl}
\alias{myDiff.to.mixmdl}
\title{obtain mixtools model from \code{myDiff} object}
\usage{
myDiff.to.mixmdl(myDiff, plot = F, main = "")
}
\arguments{
  \item{myDiff}{a \code{data.frame} object created by
  calculateDiffMeth from methylKit package and converted
  into data.frame. Required.}

  \item{plot}{to plot or not the nearest distance
  distribution.}

  \item{main}{title of the plot}
}
\description{
obtain mixtools model from \code{myDiff} object
}
\examples{
library(edmr)
library(mixtools)
data(edmr)

myMixmdl=myDiff.to.mixmdl(myDiff[1:5000,])
plotMdl1(myMixmdl, subtitle="example", cex.sub=1.2)
plotCost(myMixmdl, main="cost function")
}

