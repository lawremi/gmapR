\name{GmapGenomeDirectory-class}
\Rdversion{1.1}
\docType{class}
\alias{GmapGenomeDirectory-class}
\alias{path,GmapGenomeDirectory-method}
\alias{path,NULL-method}
\alias{genome,GmapGenomeDirectory-method}

\alias{GmapGenomeDirectory}

\title{Class \code{"GmapGenomeDirectory"}}
\description{
  The GmapGenomeDirectory class stores a path to a directory containing
  a one or more genome-specific subdirectories, each represented by a
  \code{\linkS4class{GmapGenome}}. Inside those directories are the
  files that the GMAP suite of tools uses for alignment, tallying, and
  other operations. This class is typically used to create a
  \code{GmapGenome} object. The default directory is
  \code{~/.local/share/gmap}, following the freedesktop.org XDG
  standard.
}
\section{Constructor}{
  \describe{
    \item{}{\code{GmapGenomeDirectory(path = getDefaultGmapGenomePath(),
        create = FALSE)}:
      Creates an object pointing to the directory at \code{path},
      creating it if it does not yet exist and \code{create} is \code{TRUE}.
    }
  }
}
\section{Methods}{
  \describe{
    \item{}{\code{path(object)}: gets the path to the genome directory.}
    \item{}{\code{genome(x)}: gets the names of the genomes in the
      directory.}
  }
}
\author{
  Michael Lawrence
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{GmapGenome-class}}
}

\examples{
gmapGenomePath <- file.path(getwd(), "newGmapGenomeDirectory")
gmapGenomeDirectory <- GmapGenomeDirectory(gmapGenomePath, create = TRUE)
}

\keyword{classes}
