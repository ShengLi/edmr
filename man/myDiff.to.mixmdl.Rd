\name{myDiff.to.mixmdl}
\alias{myDiff.to.mixmdl}
\title{obtain mixtools model from \code{myDiff} object}
\usage{
myDiff.to.mixmdl(myDiff, plot = F, main = "")
}
\arguments{
  \item{myDiff}{a \code{myDiff} object from
  \code{methylKit} package. Required.}

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

myMixmdl=myDiff.to.mixmdl(myDiff)
}

