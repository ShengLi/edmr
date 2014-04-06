\name{cpgi.anno}
\alias{cpgi.anno}
\title{generate CpG islands GRangesList object}
\usage{
cpgi.anno(file, shore.width = 2000, shelf.width = 2000)
}
\arguments{
  \item{file}{bed file for CpG islands}

  \item{shore.width}{width for CpG shores}

  \item{shelf.width}{width for CpG shelves}
}
\description{
generate CpG islands GRangesList object
}
\examples{
library(GenomicRanges)
library(IRanges)
cpgi.file=system.file("extdata", "chr22.hg19_cpgisland_all.bed.gz", package = "edmr")
cpgi=cpgi.anno(file=cpgi.file)
}

