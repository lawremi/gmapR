\name{bam_tally-methods}
\docType{methods}
\alias{bam_tally-methods}
\alias{bam_tally,BamFile-method}
\alias{bam_tally,character-method}
\alias{bam_tally,GmapBamReader-method}
\alias{bam_tally}

\title{Per-position Alignment Summaries}

\description{
  Given a set of alignments, for each position in the genome output
  counts for the reference allele and all alternate alleles. Often used
  as a precursor to detecting variants. Indels will be supported soon.
}

\usage{
\S4method{bam_tally}{BamFile}(x, param, ...)
\S4method{bam_tally}{character}(x, param, ...)
}

\arguments{
  \item{x}{a \code{BamFile} object or string path to a BAM file to read}
  \item{param}{The \code{\linkS4class{BamTallyParam}} object with
    parameters for the tally operation. }
  \item{...}{Arguments that override settings in \code{param}.}
}

\value{
  A \code{\link[GenomicRanges]{GRanges}}, with a range for each position
  that passed the filters, and with the following \code{elementMetadata}
  columns:
  \item{location}{A string representation of the location, of the form
    \dQuote{chr:pos}. This makes it easy, e.g., to check for the
    presence of a variant in another result object.}
  \item{ref}{The reference base at that position.}
  \item{alt}{The base for the alternate allele, \code{NA} for the
    reference allele row.}
  \item{ncycles}{The number of unique cycles at which the alternate allele was
    observed, \code{NA} for the reference allele row.}
  \item{ncycles.ref}{The number of unique cycles at which the reference
    allele was observed.}
  \item{count}{The number of reads with the alternate allele,
    \code{NA} for the reference allele row.}
  \item{count.ref}{The number of reads with the reference allele.}
  \item{count.total}{The total number of reads at that position,
    including reference and all alternates.}
  \item{high.quality}{The number of reads for the alternate allele that were
    above \code{high_quality_cutoff}, \code{NA} for the reference allele
    row.}
  \item{high.quality.ref}{The number of reads for the reference allele that were
    above \code{high_quality_cutoff}.}
  \item{mean.quality}{The mean mapping quality for the alternate allele,
    \code{NA} for the reference allele row.}
  \item{mean.quality.ref}{The mean mapping quality for the reference
    allele.}
  \item{count.pos}{The number of positive strand reads for the alternate
    allele, \code{NA} for the reference allele row.}
  \item{count.pos.ref}{The number of positive strand reads for the reference
    allele.}
  \item{count.neg}{The number of negative strand reads for the alternate
    allele, \code{NA} for the reference allele row.}
  \item{count.neg.ref}{The number of negative strand reads for the reference
    allele.}

  An additional column is present for each bin formed by
  the \code{cycle_breaks} parameter, with the read count for that bin.
}

\author{Michael Lawrence}
