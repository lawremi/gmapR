\name{gmap_build-methods}
\docType{methods}
\alias{gmap_build-methods}
\alias{gmap_build,ANY,GmapGenome-method}
\alias{gmap_build,character,GmapGenome-method}
\alias{gmap_build,DNAStringSet,GmapGenome-method}
\title{Build Gmap/Gsnap Genome}
\description{
  Construct the IIT (interval index tree) needed from the GMAP suite of tools to run from a
  genome file. IIT files are an oligomer index and what allow GMAP and GSNAP to
  efficiently lookup interval information for fast genomic
  mapping. \href{http://bioinformatics.oxfordjournals.org/content/26/7/873.abstract}{Fast
    and SNP-tolerant detection of complex variants and splicing in short
    reads} offers an depth explication of IIT files and their use in GMAP and GSNAP. 
}
\section{Methods}{
\describe{

\item{\code{signature(x = "ANY", genome = "GmapGenome")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(x = "character", genome = "GmapGenome")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(x = "DNAStringSet", genome = "GmapGenome")}}{
%%  ~~describe this method here~~
}
}}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }

\arguments{
\item{d}{genome name}
\item{D}{destination directory for installation (defaults to gmapdb directory specified at configure time}
\item{k}{k-mer value for genomic index (allowed: 12..15, default 14)}
\item{S}{do not order chromosomes in numeric/alphabetic order, but use order in FASTA file(s)}
\item{g}{files are gzipped, so need to gunzip each file first}
}

\examples{
\dontrun{flyGG <- GmapGenome(genome = "dm3",
                    directory = ggd)
gmap_build(x=Dmelanogaster, genome=flyGG)
}
}

