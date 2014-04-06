\name{genebody.anno}
\alias{genebody.anno}
\title{generate genebody GRangesList object}
\usage{
genebody.anno(file)
}
\arguments{
  \item{file}{bed file for gene models}
}
\description{
generate genebody GRangesList object
}
\examples{
library(IRanges)
genebody.file=system.file("extdata", "chr22.hg19_refseq_all_types.bed.gz", package = "edmr")
genebody=genebody.anno(file=genebody.file)
}

