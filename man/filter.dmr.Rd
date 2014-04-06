\name{filter.dmr}
\alias{filter.dmr}
\title{significant DMRs}
\usage{
filter.dmr(myDMR, DMR.qvalue = 0.001, mean.meth.diff = 20, num.CpGs = 5,
  num.DMCs = 3)
}
\arguments{
  \item{myDMR}{a \code{GRanges} object for DMR regions
  called.}

  \item{DMR.qvalue}{qvalue cutoff for DMC definition,
  default: 0.001}

  \item{mean.meth.diff}{cutoff of the DMR mean methylation
  difference, default=20}

  \item{num.CpGs}{cutoff of the number DMCs in each region
  to call DMR, default: 5}

  \item{num.DMCs}{cutoff of the number of CpGs, default: 3}
}
\description{
significant DMRs
}
\examples{
library(methylKit)
library(GenomicRanges)
library(mixtools)
library(data.table)
data(edmr)
mydmr=edmr(myDiff, mode=1, ACF=TRUE)
mysigdmr=filter.dmr(mydmr)
}

