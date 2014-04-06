\name{cpgi.anno}
\alias{cpgi.anno}
\title{generate CpG islands GRangesList object}
\usage{
cpgi.anno(file, shore.width = 2000, shelf.width = 2000)
}
\arguments{
  \item{bed}{file for CpG islands}
}
\description{
generate CpG islands GRangesList object
}
\examples{
library(GenomicRanges)
library(IRanges)
cpgi=cpgi.anno(file="http://edmr.googlecode.com/files/hg19_cpgisland_all.bed")
}

