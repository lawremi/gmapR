\name{GsnapParam-class}
\Rdversion{1.1}
\docType{class}
\alias{GsnapParam-class}
\alias{GsnapParam}

\title{Class \code{"GsnapParam"}}
\description{
  A \code{GsnapParam} object stores parameters for
  \code{\link{gsnap}}. The function of the same name serves as its
  constructor. 
}

\usage{
GsnapParam(genome, unique_only = FALSE, molecule = c("RNA", "DNA"),
           max_mismatches = NULL,
           suboptimal_levels = 0L, mode = "standard",
           snps = NULL,
           npaths = if (unique_only) 1L else 100L,
           quiet_if_excessive = unique_only, nofails = unique_only,
           split_output = !unique_only,
           novelsplicing = FALSE, splicing = NULL, 
           nthreads = 1L, part = NULL, batch = "2",
           terminal_threshold = if (molecule == "DNA") 1000L else 2L,
           gmap_mode = if (molecule == "DNA") "none" else
                       "pairsearch,terminal,improve",
           clip_overlap = FALSE, ...)
}

\arguments{
  \item{genome}{A GmapGenome object to align against}
  \item{unique_only}{Whether only alignments with a unique match should be
    output. The default is FALSE.}
  \item{molecule}{The type of molecule sequenced; used to determine
    appropriate parameter defaults.}
  \item{max_mismatches}{The maximum number of mismatches to allow per
    alignment. If NULL, then the value defaults to ((readlength + 2) / 12 - 2))}
  \item{suboptimal_levels}{Report suboptimal hits beyond best hit. The
    default is 0L.}
  \item{mode}{The alignment mode. It can be "standard", "cmet-stranded",
    "cmet-nonstranded", "atoi-stranded", or "atoi-nonstranded". The default 
    is "standard".}
  \item{snps}{If not NULL, then a GmapSnps object. Provided SNPs will not
    count as mismatches.}
  \item{npaths}{The maximum number of paths to print.}
  \item{quiet_if_excessive}{If more than maximum number of paths are
    found, then no alignment from the read will be in the output.}
  \item{nofails}{Exclude failed alignments from output}
  \item{split_output}{Basename for multiple-file output, separately for
    nomapping, halfmapping_uniq, halfmapping_mult, unpaired_uniq,
    unpaired_mult, paired_uniq, paired_mult, concordant_uniq, and
    concordant_mult results (up to 9 files, or 10 if --fails-as-input is
    selected, or 3 for single-end reads)}
  \item{novelsplicing}{Logical indicating whether to look for novel
    splicing events. FALSE is the default.}
  \item{splicing}{If not NULL, a GmapSplices object. NULL is the default.}
  \item{nthreads}{The number of worker threads gsnap should use to align.}
  \item{part}{If not NULL, then process only the i-th out of every n
    sequences e.g., 0/100 or 99/100 (useful for distributing jobs to a
    computer farm). If NULL, then all sequences are processed. NULL is
    the default.}
  \item{batch}{This argument allows control over gsnap's memory mapping
    and allocation. The default is mode 2.
    Mode 0: \{offsets=allocate, positions=mmap, genome=mmap\},
    Mode 1: \{offsets=allocate, positions=mmap & preload,genome=mmap\},
    Mode 2: \{offsets=allocate, positions=mmap & preload,genome=mmap &
    preload\}, 
    Mode 3: \{offsets=allocate, positions=allocate,genome=mmap & preload\},
    Mode 4: \{offsets=allocate, positions=allocate,genome=allocate\}}
  \item{...}{Additional parameters for gsnap. See gsnap's full
    documentation for those available.}
  \item{terminal_threshold}{If this number of mismatches is exceeded,
    GSNAP will attempt to align from one of the sequence, eventually
    giving up and discarding the rest of the sequence. This is called a
    \dQuote{terminal alignment}. By setting this to a high value, we
    have effectively disabled it for DNA, since terminal alignments were
    motivated by splicing alignment problems and other special cases.
  }
  \item{gmap_mode}{Specifies the GMAP pipeline executed when GSNAP
    delegates to GMAP (a Smith-Waterman aligner) in difficult cases. We
    have disabled this for DNA, since such difficult cases are only
    anticipated in the context of splicing or complex rearrangements.}
  \item{clip_overlap}{Whether to equally clip paired ends that
    overlap each other (due to the fragment length being shorter than 2X
    the read length). This can be important for getting accurate counts
    from \code{bam_tally}.
  }
}

\seealso{
  \code{\link{gsnap}}
}
