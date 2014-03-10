\name{TP53Genome}
\alias{TP53Genome}
\alias{TP53Which}
\title{
  Demo genome around TP53
}
\description{
  Returns a \code{\linkS4class{GmapGenome}} object consisting of the
  UCSC hg19 sequence centered on the region of the TP53 gene, with 1 Mb
  flanking sequence on each side. This is intended as a
  test/demonstration genome and can be used, e.g., in conjunction with
  the \code{LungCancerLines} data package.
}
\usage{
TP53Genome()
TP53Which()
}
\value{
  For \code{TP53Genome}, a \code{GmapGenome} object. If this is the
  first time the user has run this function, a side-effect will be the
  generation of an on-disk genome index, under the name
  \dQuote{TP53_demo_VERSION} in the default genome directory, where
  \code{VERSION} is the version of the TxDb package providing the bounds
  of the P53 gene.

  For \code{TP53Which}, a \code{GRanges} of the extents of the TP53
  gene, translated to the space of \code{TP53Genome}.
}
\author{
  Michael Lawrence, Cory Barr
}

\examples{
TP53Genome()
}
