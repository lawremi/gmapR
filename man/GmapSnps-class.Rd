\name{GmapSnps-class}
\Rdversion{1.1}
\docType{class}
\alias{GmapSnps-class}
\alias{directory,GmapSnps-method}
\alias{GmapSnps}

\title{Class \code{"GmapSnps"}}
\description{
  This class represents a set of SNPs (single nucleotide polymorphisms)
  for use with GMAP and GSNAP (typically for SNP-tolerant alignment.)
}
\section{Objects from the Class}{
  ##TODO: doc these args
Objects can be created by calls of the form \code{GmapSnps(snps,
  directory, name, create)}.
}

\usage{GmapSnps(snps, directory, name = snps, create = FALSE, ...)}

\arguments{
  \item{snps}{A path to a VCF file}
  \item{directory}{The directory to create the IIT files used by GMAP
  and GSNAP}
  \item{name}{If provided, the name to give the database of SNPs. If not
  provided, defauts to the \code{snps} argument.}
  \item{create}{If the directory provided in the \code{directory}
  argument does not exist, create it.}
  }
  
\section{Accessors}{
    \describe{
      \item{}{\code{name(x)}: returns the name of the \code{GmapSnps} object}
      \item{}{\code{directory(x)}: returns the \code{GmapGenomeDirectory}
	that is the parent of the directory containing the index files for
	this \code{GmapSnps} object.}
    }
}

  
\section{Methods}{
  \describe{
    \item{directory}{\code{signature(x = "GmapSnps")}: ... }
  }
}
\author{
  Michael Lawrence
}
\keyword{classes}
