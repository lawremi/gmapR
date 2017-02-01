\name{gsnap-methods}
\docType{methods}
\alias{gsnap-methods}
\alias{gsnap,character,character_OR_NULL,GsnapParam-method}
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
\S4method{gsnap}{character,character_OR_NULL,GsnapParam}(input_a, input_b, params,
                   output = file.path(getwd(),
                                      file_path_sans_ext(basename(input_a),
                                                         TRUE)),
                   consolidate = TRUE, ...)
}

\arguments{
  \item{input_a}{A path to the FASTA file containing reads to align against a
  \code{GmapGenome} object. If the sequencing data is single-end, this
  is the only FASTA file used as input.}
  \item{input_b}{If provided, a path to the FASTA file containing the
  second set of reads from paired-end sequencing data.}
  \item{params}{A \code{GsnapParam} object to configure the behavior of
    GSNAP.}
  \item{output}{The output path for the GSNAP alignments. The results
    will be saved in \code{dirname(output)}. If \code{split_output} in
    \code{params} is \code{TRUE}, \code{basename(output)} is used as the
    common stem for the multiple output files. Otherwise, the results
    are saved to a single SAM file, its path formed by adding the
    \dQuote{sam} extention to \code{output}.
  }
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
