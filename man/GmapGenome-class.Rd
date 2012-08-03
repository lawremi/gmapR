\name{GmapGenome-class}
\Rdversion{1.1}
\docType{class}
\alias{GmapGenome-class}
\alias{path,GmapGenome-method}
\alias{genome,GmapGenome-method}
\alias{seqinfo,GmapGenome-method}
\alias{GmapGenome}

\title{Class \code{"GmapGenome"}}
\description{
  The GmapGenome class represents a genome that has been indexed for use
  with the GMAP suite of tools. It is typically used as a parameter to
  the functions \code{gsnap} and \code{\link{bam_tally}}. This class
  also provides the means to index new genomes, from either a FASTA file
  or a \code{BSgenome} object. Genome indexes are typically stored in a
  centralized directory on the file system and are identified by a
  string key.
}

\section{Constructor}{
  \describe{
    \item{}{\code{GmapGenome(genome, directory =
        GmapGenomeDirectory(create = create), name = genomeName(genome),
        create = FALSE, ...)}:
      
      Creates a \code{GmapGenome} corresponding to the \code{genome}
      argument, which may be either a string identifier of the genome
      within \code{directory}, a \code{\link[rtracklayer]{FastaFile}} or
      \code{\link[Biostrings]{DNAStringSet}} of the genome sequence, or
      a \code{\link[BSgenome:BSgenome-class]{BSgenome}} object.

      The genome index is stored in \code{directory} argument, which may
      be either a \code{\linkS4class{GmapGenomeDirectory}} object, or a
      string path.

      The \code{name} argument is the actual key used for storing the
      genome index within \code{directory}. If \code{genome} is a
      string, it is taken as the key. If a \code{FastaFile}, it is the
      basename of the file without the extension. If a \code{BSgenome},
      it is the \code{providerVersion}. Otherwise, the \code{name} must
      be specified. If \code{create} is \code{TRUE}, the genome index is
      created if one with that name does not already exist. This
      obviously only works if \code{genome} actually contains the genome
      sequence.

      The first example below gives the typical and recommended usage
      when implementing a reproducible analysis.
    }
  }
}

\section{Accessors}{
  \describe{
    \item{}{\code{path(object)}: returns the path to the directory
      containing the genome index files.
    }
    \item{}{\code{directory(x)}: returns the \code{GmapGenomeDirectory}
      that is the parent of the directory containing the index files for
      this genome.
    }
    \item{}{\code{genome(x)}: gets the name of this genome.}
    \item{}{\code{seqinfo(x)}: gets the \code{\link[GenomicRanges]{Seqinfo}}
      for this genome; only sequence names and lengths are available.
    }
  }
}

\author{
  Michael Lawrence
}
\examples{
\dontrun{
library(BSgenome.Dmelanogaster.UCSC.dm3)
flyGG <- GmapGenome(Dmelanogaster, create = TRUE)

## access system-wide genome using a key
flyGG <- GmapGenome(genome = "dm3",
                    directory = path.expand("/usr/share/gmap"))
}
}
\keyword{classes}