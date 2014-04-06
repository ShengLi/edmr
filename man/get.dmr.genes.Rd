\name{get.dmr.genes}
\alias{get.dmr.genes}
\title{get gene list based the genebody granges}
\usage{
get.dmr.genes(myDMR, subject, id.type = "gene.symbol")
}
\arguments{
  \item{myDMR}{DMRs predicted by \code{edmr}.}

  \item{subject}{GRanges used to annotate DMRs. For
  example, genebody$promoter will annotate the DMRs using
  promoters}

  \item{id.type}{the column names that will be used to
  annotate the DMR. default: "gene.symbol"}
}
\description{
get gene list based the genebody granges
}

