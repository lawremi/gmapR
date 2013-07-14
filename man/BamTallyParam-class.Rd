\name{BamTallyParam-class}
\Rdversion{1.1}
\docType{class}
\alias{BamTallyParam-class}
\alias{coerce,BamTallyParam,list-method}
\alias{as.list,BamTallyParam-method}
\alias{BamTallyParam}

\title{Class \code{"BamTallyParam"}}
\description{
  A \code{BamTallyParam} object stores parameters for
  \code{\link{bam_tally}}. The function of the same name serves as its
  constructor. 
}

\usage{
  BamTallyParam(genome, which = RangesList(), read_pos_breaks = NULL,
                high_base_quality = 0L,
                minimum_mapq = 0L,
                concordant_only = FALSE, unique_only = FALSE,
                primary_only = FALSE, ignore_duplicates = FALSE,
                min_depth = 0L, variant_strand = 0L,
                ignore_query_Ns = FALSE,
                indels = FALSE)
}
\arguments{
  \item{genome}{A \code{GmapGenome} object, or something coercible to one.}
  \item{which}{A \code{RangesList} or something coercible to
    one that limits the tally to that range or set of ranges. By
    default, the entire genome is processed.
  }
  \item{read_pos_breaks}{The breaks, like those passed to \code{\link{cut}}
    for aggregating the per-read position counts. If \code{NULL}, no per-cycle
    counts are returned.}
  \item{high_base_quality}{The minimum mapping quality for a
    read to be counted as high quality.}
  \item{minimum_mapq}{Minimum mapping quality for a read to be counted
    at all.}
  \item{concordant_only}{Consider only what gnsap
    calls \dQuote{concordant} alignments.
  }
  \item{unique_only}{Consider only the uniquly mapped reads.}
  \item{primary_only}{Consider only primary pairs.}
  \item{ignore_duplicates}{Whether to ignore the reads flagged as
    PCR/optical duplicates.
  }
  \item{min_depth}{The minimum number of reads overlapping a position for
    it to be counted.}
  \item{variant_strand}{The number of strands on which a variant must be
    seen for it to be counted. This means that a value of 0 will report
    reference alleles in addition to variants. A value of 1 will report
    only positions where a variant was seen on at least one strand, and
    2 requires the variant be seen on both strands. Setting this to 1
    is a good way to save resources.}
  \item{ignore_query_Ns}{Whether to ignore the N base pairs when
    counting. Can save a lot of resources when processing low quality data.}
  \item{indels}{Whether to return indel counts; not supported yet.}
}

\seealso{
  \code{\link{bam_tally}}
}
