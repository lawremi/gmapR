\name{GmapSnpDirectory-class}
\Rdversion{1.1}
\docType{class}
\alias{GmapSnpDirectory-class}
\alias{GmapSnpDirectory}
\alias{[[<-,GmapSnpDirectory,ANY,ANY-method}
\alias{length,GmapSnpDirectory-method}
\alias{names,GmapSnpDirectory-method}
\alias{path,GmapSnpDirectory-method}

\title{Class \code{"GmapSnpDirectory"}}
\description{
  This class represents a directory containig one or more sets of SNPs,
  each corresponding to a genome. These SNP databases enable
  SNP-tolerant alignment with GMAP and GSNAP. If the underlying files
  have not been created, this class provides a means to do so.
}
\section{Methods}{
  \describe{
    \item{[[<-}{\code{signature(x = "GmapSnpDirectory", i = "ANY", j = "ANY")}: ... }
    \item{length}{\code{signature(x = "GmapSnpDirectory")}: ... }
    \item{names}{\code{signature(x = "GmapSnpDirectory")}: ... }
    \item{path}{\code{signature(object = "GmapSnpDirectory")}: ... }
	 }
}
\author{
  Michael Lawrence
}
\keyword{classes}
