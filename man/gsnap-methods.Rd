\name{gsnap-methods}
\docType{methods}
\alias{gsnap-methods}
\alias{gsnap,character,characterORNULL,GsnapParam-method}
\alias{gsnap}

\title{Align a Set of Reads Using the GSNAP Aligner}

\description{
  Given a set of alignments, align them to a genome using the GSNAP
  algorithm. The GSNAP algorithm contains a number of features making it
  a very high quality algorithm for dealing with short reads and those
  from RNA-seq data in particular. Via the \code{GsnapParam} class
  and the \code{gsnap} function, R users are given complete control over GSNAP.
}

\usage{
\S4method{gsnap}{character,characterORNULL,GsnapParam}(input_a, input_b, params,
                   output = file_path_sans_ext(input_a, TRUE),
                   consolidate = TRUE, ...)
}

\arguments{
  \item{input_a}{A path to the FASTA file containing reads to align against a
  \code{GmapGenome} object. If the sequencing data is single-end, this
  is the only FASTA file used as input.}
  \item{input_b}{If provided, a path to the FASTA file containing the
  second set of reads from paired-end sequencing data.}
  \item{params}{A \code{GsnapParam} object to configure the behavior of GSNAP.}
  \item{output}{The output directory for the GSNAP alignments. If not
  provided, the directory is the same directory that contains input_a.}
  \item{consolidate}{If GSNAP is run with multiple worker threads, each
  thread will output its own set of files. If consolidate is set to
  TRUE, these files will be merged. The default is TRUE.}
  \item{...}{Additional arguments to pass to GSNAP not specifically
  supported by the \code{gmapR} package.}
}

\value{
  A \code{GsnapOutput} class.
}

\author{Michael Lawrence}
