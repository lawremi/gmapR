\name{bam_tally-methods}
\docType{methods}
\alias{bam_tally-methods}
\alias{bam_tally,BamFile-method}
\alias{bam_tally,character-method}
\alias{bam_tally,GmapBamReader-method}
\alias{genome,TallyIIT-method}
\alias{bam_tally}
\alias{variantSummary}

\title{Per-position Alignment Summaries}

\description{
  Given a set of alignments, for each position in the genome output
  counts for the reference allele and all alternate alleles. Often used
  as a precursor to detecting variants. Indels will be supported soon.
}

\usage{
\S4method{bam_tally}{BamFile}(x, param, ...)
\S4method{bam_tally}{character}(x, param, ...)
variantSummary(x, read_pos_breaks = NULL, high_base_quality = 0L,
               keep_ref_rows = FALSE)
}

\arguments{
  \item{x}{a \code{BamFile} object or string path to a BAM file to read}
  \item{param}{The \code{\linkS4class{BamTallyParam}} object with
    parameters for the tally operation. }
  \item{read_pos_breaks}{The breaks, like those passed to \code{\link{cut}}
    for aggregating the per-read position counts. If \code{NULL}, no per-cycle
    counts are returned.}
  \item{high_base_quality}{The minimum mapping quality for a
    read to be counted as high quality.}
  \item{keep_ref_rows}{Whether to keep the rows describing only the
    reference calls, i.e., where ref and alt are the same. These are
    useful when one needs the reference counts even when there are no
    alts at that position.}
  \item{...}{Arguments that override settings in \code{param}.}
}

\value{
  The \code{bam_tally} function returns an opaque pointer to a C-level
  data structure with the class \dQuote{TallyIIT}. Currently, the only
  operation applicable to this object is \code{variantSummary}.
  
  The \code{variantSummary} function returns
  a \code{\link[VariantAnnotation]{VRanges}}, with a range for each position
  that passed the filters. The depth columns correspond to the counts
  after quality filtering (except for indels, for which there is no
  quality filtering). The following \code{elementMetadata}
  columns are also present:
  \item{raw.count}{The number of reads with the alternate allele,
    \code{NA} for the reference allele row.}
  \item{raw.count.ref}{The number of reads with the reference allele.}
  \item{raw.count.total}{The total number of reads at that position,
    including reference and all alternates.}
  \item{mean.quality}{The mean base quality for the alt allele,
    truncated at \code{high_base_quality}.}
  \item{mean.quality.ref}{The mean base quality for the ref allele,
    truncated at \code{high_base_quality}.}
  \item{count.plus}{The number of positive strand reads for the alternate
    allele, \code{NA} for the reference allele row.}
  \item{count.plus.ref}{The number of positive strand reads for the reference
    allele.}
  \item{count.minus}{The number of negative strand reads for the alternate
    allele, \code{NA} for the reference allele row.}
  \item{count.minus.ref}{The number of negative strand reads for the reference
    allele.}
  \item{read.pos.mean}{Mean read position for the alt allele.}
  \item{read.pos.mean.ref}{Mean read position for the ref allele.}
  \item{read.pos.var}{Variance in the read positions for the alt allele.}
  \item{read.pos.var.ref}{Variance in the read positions for the ref allele.}
  \item{mdfne}{Median distance from nearest end for the alt allele.}
  \item{mdfne.ref}{Median distance from nearest end for the ref allele.}
  
  An additional column is present for each bin formed by
  the \code{read_pos_breaks} parameter, with the read count for that bin.
}

\seealso{\code{tallyVariants} in the VariantTools package provides a
  high-level wrapper for this functionality.}

\author{Michael Lawrence}
