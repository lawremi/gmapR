\name{GsnapOutput-class}
\Rdversion{1.1}
\docType{class}
\alias{GsnapOutput-class}
\alias{bamPaths,GsnapOutput-method}
\alias{path,GsnapOutput-method}
\alias{GsnapOutput}

\title{Class \code{"GsnapOutput"}}
\description{
  A \code{GsnapOutput} object stores locations of data output by the
  GSNAP alignment algorithm.
}

\section{Objects from the Class}{
  \code{GsnapOutput} objects are created from the \link{gsnap} function,
  though the function \code{GsnapOutput} can also be used as a constructor.
}

\section{Coercion}{
  In the code snippets below, \code{x} is a GsnapOutput object.
  
   \describe{
     \item{}{
       \code{as(x, BamFile)},
       \code{as(x, BamFileList)}:

       Returns either a BamFile or BamFileList object containing paths
     to the output of GSNAP.
   }
   \item{}{
     \code{asBam(x)}:

     converts all gsnap SAM files to BAM files and creates the .bai
     index files.
   }
 }
}

\author{
  Michael Lawrence
}
\seealso{
\code{\link{gsnap}}
}
\keyword{classes}
