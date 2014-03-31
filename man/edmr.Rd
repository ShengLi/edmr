\name{edmr}
\alias{edmr}
\title{Empirical differentially methylated regions}
\usage{
edmr(myDiff, step = 100, dist = "none", DMC.qvalue = 0.01,
  DMC.methdiff = 25, num.DMCs = 1, num.CpGs = 3, DMR.methdiff = 20,
  plot = FALSE, main = "", mode = 1, ACF = TRUE, fuzzypval = 1)
}
\arguments{
  \item{myDiff}{a \code{myDiff} object from
  \code{methylKit} package. Required.}

  \item{step}{a numeric variable for calculating
  auto-correlation, default: 100.}

  \item{dist}{distance cutoff to call a gap for DMR,
  default: "none", which will be automatically determined
  by the bimodal normal distribution, default: 100.}

  \item{DMC.qvalue}{qvalue cutoff for DMC definition,
  default: 0.01}

  \item{DMC.methdiff}{methylation difference cutoff for DMC
  definition, default: 25.}

  \item{num.DMCs}{cutoff of the number DMCs in each region
  to call DMR, default: 1.}

  \item{num.CpGs}{cutoff of the number of CpGs, default:
  3.}

  \item{DMR.methdiff}{cutoff of the DMR mean methylation
  difference, default=20.}

  \item{plot}{plot the bimodal normal distribution fitting
  or not, default=FAlSE.}

  \item{main}{the title of the plot, if plot=TRUE.
  Default=FALSE.}

  \item{mode}{the mode of call DMRs. 1: using all CpGs
  together. 2: use unidirectional CpGs to call DMRs.
  default: 1.}

  \item{ACF}{p-value combination test with (TRUE, default)
  or without (FALSE) dependency adjustment.}

  \item{fuzzypval}{p-value cutoff for raw regions
  definition.}
}
\value{
\code{GRanges}
}
\description{
Comprehensive DMR analysis based on bimodal normal
distribution model and weighted cost function for regional
methylation analysis optimization. It captures the regional
methylation modification by taking the spatial distribution
of CpGs into account for the enrichment DNA methylation
sequencing data so as to optimize the definition of the
empirical regions. Combined with the dependent adjustment
for regional p-value combination.
}
\examples{
library(methylKit)
library(GenomicRanges)
library(mixtools)
library(data.table)
data(edmr)
mydmr=edmr(chr22.myDiff, mode=1, ACF=TRUE)
mysigdmr=filter.dmr(mydmr)
mydmr2=edmr(chr22.myDiff, mode=2, ACF=TRUE)
mydmr3=edmr(chr22.myDiff, mode=1, ACF=FALSE)
mydmr4=edmr(chr22.myDiff, mode=2, ACF=FALSE)
}

