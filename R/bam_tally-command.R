### =========================================================================
### bam_tally program
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Raw tally result
###

setClass("TallyIIT", representation(genome = "GmapGenome",
                                    bam = "BamFile",
                                    xs = "logical",
                                    read_pos = "logical",
                                    nm = "logical",
                                    codon = "logical"),
         contains = "IIT")

TallyIIT <- function(ptr, genome, bam, xs, read_pos, nm, codon) {
  new("TallyIIT", ptr = ptr, genome = genome, bam = bam, xs = xs,
      read_pos = read_pos, nm = nm, codon = codon)
}

setMethod("genome", "TallyIIT", function(x) x@genome)

bamFile <- function(x) x@bam

setMethod("show", "TallyIIT", function(object) {
  cat("Tally IIT object for '", path(bamFile(object)), "' on '",
      genome(genome(object)), "'\n", sep = "")
})

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### High-level wrapper
###

setGeneric("bam_tally",
           function(x, param, ...)
           standardGeneric("bam_tally"),
           signature = "x")

setMethod("bam_tally", "BamFile",
          function(x, param, ...)
          {
            x <- GmapBamReader(x)
            callGeneric()
          })

setMethod("bam_tally", "character",
          function(x, param, ...)
          {
            x <- BamFile(x)
            callGeneric()
          })

setMethod("bam_tally", "GmapBamReader",
          function(x, param, ...)
          {
            param_list <- as.list(param)
            args <- list(...)
            param_list[names(args)] <- args
            genome <- param_list$genome

            ##verify genome has been created

            param_list$db <- genome(genome)
            param_list$genome_dir <- path(directory(genome))
            if (!.gmapGenomeCreated(genome)) {
              stop("The GmapGenome object has not yet been created. ",
                   "One solution is to run the GmapGenome constructor ",
                   "with create=TRUE")
            }

            param_list$genome <- NULL
            
            TallyIIT(do.call(.bam_tally_C, c(list(x), param_list)), genome,
                     as(x, "BamFile"), xs=param_list$xs,
                     read_pos=param_list$read_pos, nm=param_list$nm,
                     codon=!is.null(param_list$exon_iit))
          })

variantSummary <- function(x, read_pos_breaks = NULL,
                           keep_ref_rows = FALSE, read_length = NA_integer_,
                           high_nm_score = NA_integer_)
{
  read_length <- as.integer(read_length)
  if (length(read_length) != 1L) {
    stop("'read_length' must be a single integer")
  }
  high_nm_score <- as.integer(high_nm_score)
  stopifnot(length(high_nm_score) == 1L)
  if (!is.na(high_nm_score) && !x@nm) {
      stop("'high_nm_score is not NA but NM scores were not tallied")
  }
  if (length(read_pos_breaks) > 0L) {
    if (!x@read_pos) {
      stop("'read_pos_breaks' non-empty, but read positions were not tallied")
    }
    read_pos_breaks <- as.integer(read_pos_breaks)
    if (any(is.na(read_pos_breaks)))
        stop("'read_pos_breaks' should not contain missing values")
    if (length(read_pos_breaks) < 2)
        stop("'read_pos_breaks' needs at least two elements to define a bin")
    if (is.unsorted(read_pos_breaks))
        stop("'read_pos_breaks' must be sorted")
  }

  tally <- .Call(R_tally_iit_parse, x@ptr,
                 read_pos_breaks,
                 NULL, read_length, x@xs, high_nm_score)
  
  tally_names <- c("seqnames", "pos", "ref", "alt",
                   "n.read.pos", "n.read.pos.ref",
                   "count", "count.ref", "raw.count.total", "count.total",
                   "count.plus", "count.plus.ref",
                   "count.minus", "count.minus.ref",
                   "count.del.plus", "count.del.minus",
                   "read.pos.mean", "read.pos.mean.ref",
                   "read.pos.var", "read.pos.var.ref",
                   "mdfne", "mdfne.ref", "codon.strand",
                   "count.xs.plus", "count.xs.plus.ref",
                   "count.xs.minus", "count.xs.minus.ref",
                   "count.high.nm", "count.high.nm.ref")

  break_names <- character()
  if (length(read_pos_breaks) > 0L) {
    read_pos_breaks <- as.integer(read_pos_breaks)
    break_names <- paste("readPosCount", head(read_pos_breaks, -1),
                         tail(read_pos_breaks, -1), sep = ".")
    tally_names <- c(tally_names, break_names)
  }
  names(tally) <- tally_names
  tally$codon.strand <- if (x@codon) {
      structure(tally$codon.strand, class="factor", levels=levels(strand()))
  }
  tally <- Filter(Negate(is.null), tally)

  if (!keep_ref_rows) {
    variant_rows <- !is.na(tally$alt)
    if (!all(variant_rows))
      tally <- lapply(tally, `[`, variant_rows)
  }
  
  meta_names <- setdiff(names(tally),
                        c("seqnames", "pos", "ref", "alt", "count",
                          "count.ref", "count.total"))
  genome <- genome(x)
  indel <- nchar(tally$ref) == 0L | nchar(tally$alt) == 0L
  metacols <- DataFrame(tally[meta_names])
  mcols(metacols) <- variantSummaryColumnDescriptions(read_pos_breaks, x@xs,
                                                      high_nm_score, x@codon)

  gr <- with(tally,
             VRanges(seqnames,
                     IRanges(pos,
                             width = ifelse(nchar(alt) == 0 & !is.na(alt),
                                            nchar(ref), 1L)),
                     ref, alt,
                     count.total, count.ref, count,
                     seqlengths = seqlengths(genome)))
  mcols(gr) <- metacols
  checkTallyConsistency(gr)
  ## important to preserve seqlevel ordering compatible with 'genome'
  seqinfo(gr) <- merge(seqinfo(genome), seqinfo(bamFile(x)))
  gr <- keepSeqlevels(gr, intersect(seqlevels(gr), seqlevels(bamFile(x))))
  gr <- normalizeIndelAlleles(gr, genome)
  gr
}

checkTallyConsistency <- function(x) {
  with(mcols(x), {
    stopifnot(all(count.plus + count.minus == altDepth(x), na.rm=TRUE))
    stopifnot(all(count.plus.ref + count.minus.ref == refDepth(x)))
  })
}

normalizeIndelAlleles <- function(x, genome) {
  is.indel <- nchar(ref(x)) == 0L | (nchar(alt(x)) == 0L & !is.na(alt(x)))
  if (any(is.indel)) {
    indels <- x[is.indel]
    flanks <- flank(indels, 1)
    anchor <- getSeq(genome, flanks)
    ref(x)[is.indel] <- paste0(anchor, ref(indels))
    alt(x)[is.indel] <- paste0(anchor, alt(indels))
    ranges(x)[is.indel] <- resize(ranges(flanks), nchar(ref(x)[is.indel]))
  }
  x
}

guessReadLengthFromBam <- function(x, n=100L) {
  ga <- readGAlignments(BamFile(x, yieldSize=n))
  readlen <- unique(qwidth(ga))
  if (length(readlen) != 1L)
    NA
  else readlen
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level wrappers
###

normArgSingleInteger <- function(x) {
  name <- deparse(substitute(x))
  x <- as.integer(x)
  if (!isSingleInteger(x))
    stop("'", name, "' should be a single, non-NA integer")
  x
}
normArgSingleNumber <- function(x) {
    name <- deparse(substitute(x))
    x <- as.numeric(x)
    if (!isSingleNumber(x))
        stop("'", name, "' should be a single, non-NA number")
    x
}
normArgTRUEorFALSE <- function(x) {
  name <- deparse(substitute(x))
  if (!isTRUEorFALSE(x))
    stop("'", name, "' should be TRUE or FALSE")
  x
}

normArgSingleCharacterOrNULL <- function(x) {
  name <- deparse(substitute(x))
  if (!is.null(x) && (!is(x, "character") || length(x) != 1))
    stop("'", name, "' should be a single character value")
  x
}

.bam_tally_C <- function(bamreader, genome_dir = NULL, db = NULL,
                         which = NULL,
                         high_base_quality = 0L, desired_read_group = NULL,
                         alloclength = 200000L,
                         minimum_mapq = 0L, good_unique_mapq = 35L,
                         maximum_nhits = 1000000L,
                         concordant_only = FALSE, unique_only = FALSE,
                         primary_only = FALSE, ignore_duplicates = FALSE,
                         min_depth = 0L, variant_strand = 0L, variant_pct = 0,
                         ignore_query_Ns = FALSE,
                         indels = FALSE,
                         blocksize = 1000L, verbosep = FALSE,
                         min_softclip = 0L, max_softclip = 0L,
                         exon_iit = NULL, xs = FALSE, read_pos = FALSE,
                         min_base_quality = 0L, noncovered = FALSE, nm = FALSE)
{
  if (!is(bamreader, "GmapBamReader"))
    stop("'bamreader' must be a GmapBamReader")
  if (length(which) > 0L) {
    which <- list(as.character(seqnames(which)), start(which), end(which))
  } else which <- NULL
  if (!is.null(genome_dir) && !isSingleString(genome_dir))
    stop("'genome_dir' must be NULL or a single, non-NA string")
  if (!is.null(db) && !isSingleString(db))
    stop("'db' must be NULL or a single, non-NA string")
  if (!is.null(desired_read_group) && !isSingleString(desired_read_group))
    stop("'desired_read_group' must be NULL or a single, non-NA string")
  .Call(R_Bamtally_iit, bamreader@.extptr, genome_dir, db, which,
        desired_read_group,
        normArgSingleInteger(alloclength),
        normArgSingleInteger(minimum_mapq),
        normArgSingleInteger(good_unique_mapq),
        normArgSingleInteger(maximum_nhits),
        normArgTRUEorFALSE(concordant_only),
        normArgTRUEorFALSE(unique_only),
        normArgTRUEorFALSE(primary_only),
        normArgTRUEorFALSE(ignore_duplicates),
        normArgSingleInteger(min_depth),
        normArgSingleInteger(variant_strand),
        normArgSingleNumber(variant_pct),
        normArgTRUEorFALSE(ignore_query_Ns),
        normArgTRUEorFALSE(indels),
        normArgSingleInteger(blocksize),
        normArgTRUEorFALSE(verbosep),
        normArgSingleInteger(min_softclip),
        normArgSingleInteger(max_softclip),
        normArgSingleCharacterOrNULL(exon_iit),
        normArgTRUEorFALSE(xs),
        normArgTRUEorFALSE(read_pos),
        normArgSingleInteger(min_base_quality),
        normArgTRUEorFALSE(noncovered),
        normArgTRUEorFALSE(nm))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Column metadata
###

variantSummaryColumnDescriptions <-
    function(read_pos_breaks, xs, high_nm, codon)
{
  desc <- c(
    n.read.pos = "Number of unique read positions for the ALT",
    n.read.pos.ref = "Number of unique read positions for the REF",
    raw.count.total = "Total depth, before quality filtering",  
    count.plus = "Positive strand ALT count",
    count.plus.ref = "Positive strand REF count",
    count.minus = "Negative strand ALT count",
    count.minus.ref = "Negative strand REF count",
    count.del.plus = "Plus strand deletion count",
    count.del.minus = "Count strand deletion count",
    read.pos.mean = "Average read position for the ALT",
    read.pos.mean.ref = "Average read position for the ALT",
    read.pos.var = "Variance in read position for the ALT",
    read.pos.var.ref = "Variance in read position for the REF",
    mdfne = "Median distance from nearest end of read for the ALT",
    mdfne.ref = "Median distance from nearest end of read for the REF",
      if (codon) {
          c(codon.strand = "Strand of transcription for the codon")
      },
      if (xs) {
          c(count.xs.plus = "Plus strand XS counts",
            count.xs.plus.ref = "Plus strand reference XS counts",
            count.xs.minus = "Minus strand XS counts",
            count.xs.minus.ref = "Minus strand reference XS counts")
      },
      count.high.nm = paste("Count of reads with >=", high_nm,
          "NM score for the ALT"),
      count.high.nm.ref = paste("Count of reads with >=", high_nm,
          "NM score for the REF"))
  if (length(read_pos_breaks) > 0L) {
    break_desc <- paste0("Raw ALT count in read position range [",
                         head(read_pos_breaks, -1), ",",
                         tail(read_pos_breaks, -1), ")")
    desc <- c(desc, break_desc)
  }
  DataFrame(Description = desc)
}
