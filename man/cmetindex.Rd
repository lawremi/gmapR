\name{cmetindex}
\alias{cmetindex}
\title{
  Call the cmetindex command
}
\description{
  Call the GMAP \command{cmetindex} command to build an index suitable
  for alignment of bisulfite-treated DNA, by allowing for C->T and G->A
  differences.

}
\usage{
cmetindex(db, use_snps = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{db}{
    The GmapGenome object
  }
  \item{use_snps}{
    A GmapSnps object for generating a SNP-tolerant index
  }
}
\author{
  Michael Lawrence
}

\examples{
\dontrun{
    library(BSgenome.Dmelanogaster.UCSC.dm3)
    flyGG <- GmapGenome(Dmelanogaster, create = TRUE)
    cmetindex(flyGG)
}
}
