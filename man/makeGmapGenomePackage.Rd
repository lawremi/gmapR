\name{makeGmapGenomePackage}
\alias{makeGmapGenomePackage}

\title{
  Function to create a GmapGenome package from a GmapGenome object
}
\description{
  A \code{GmapGenome} object is required to align reads using the GSNAP or GMAP
  algorithms. The \code{makeGmapGenomePackage} function allows users to
  save a particular \code{GmapGenome} object in an R package.
}
\usage{
makeGmapGenomePackage(gmapGenome, version, maintainer, author,
destDir = ".", license = "Artistic-2.0", pkgName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gmapGenome}{A \link{GmapGenome} object.}
  \item{version}{The version number of this package.}
  \item{maintainer}{The maintainer of the package. The string must
  contain a valid email address.}
  \item{author}{The author of the package}
  \item{destDir}{The path that the new \code{GmapGenome} package should
  be created at.}
  \item{license}{The package's license (and its version)}
  \item{pkgName}{The name the package should have. Though free form,
  names of the form GmapGenome.Organism.Source.Build are
  recommended. E.g., GmapGenome.Hsapiens.UCSC.hg19}
}

\author{
  Cory Barr
}

\seealso{
\code{\link{GmapGenome}}
}

\examples{
\dontrun{
library(gmapR)

if (!require(BSgenome.Dmelanogaster.UCSC.dm3)) {
  library(BiocManager)
  BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm3")
  library(BSgenome.Dmelanogaster.UCSC.dm3)
}

gmapGenomePath <- file.path(getwd(), "flyGenome")
if (file.exists(gmapGenomePath)) unlink(gmapGenomePath, recursive=TRUE)
ggd <- GmapGenomeDirectory(gmapGenomePath, create = TRUE)
gmapGenome <- GmapGenome(genome=Dmelanogaster,
                         directory = ggd,
                         name = "dm3",
                         create = TRUE)

makeGmapGenomePackage(gmapGenome=gmapGenome,
                      version="0.1.0",
                      maintainer="<your.name@somewhere.com>",
                      author="Your Name",
                      destDir=".",
                      license="Artistic-2.0",
                      pkgName="GmapGenome.Dmelanogaster.UCSC.dm3")

}
}