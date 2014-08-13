### =========================================================================
### bam_tally program
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Raw tally result
###

setClass("TallyIIT", representation(ptr = "externalptr",
                                    genome = "GmapGenome",
                                    bam = "BamFile"))

TallyIIT <- function(ptr, genome, bam) {
  new("TallyIIT", ptr = ptr, genome = genome, bam = bam)
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
                     as(x, "BamFile"))
          })

variantSummary <- function(x, read_pos_breaks = NULL, high_base_quality = 0L,
                           keep_ref_rows = FALSE, read_length = NA_integer_)
{
  read_length <- as.integer(read_length)
  if (length(read_length) != 1L) {
    stop("'read_length' must be a single integer")
  }
  tally <- .Call(R_tally_iit_parse, x@ptr,
                 read_pos_breaks,
                 normArgSingleInteger(high_base_quality),
                 NULL, read_length)
  
  tally_names <- c("seqnames", "pos", "ref", "alt",
                   "n.read.pos", "n.read.pos.ref",
                   "raw.count", "raw.count.ref",
                   "raw.count.total",
                   "high.quality", "high.quality.ref",
                   "high.quality.total", "mean.quality",
                   "mean.quality.ref",
                   "count.plus", "count.plus.ref",
                   "count.minus", "count.minus.ref",
                   "read.pos.mean", "read.pos.mean.ref",
                   "read.pos.var", "read.pos.var.ref",
                   "mdfne", "mdfne.ref", "codon.dir")
  break_names <- character()
  if (length(read_pos_breaks) > 0L) {
    read_pos_breaks <- as.integer(read_pos_breaks)
    break_names <- paste("readPosCount", head(read_pos_breaks, -1),
                         tail(read_pos_breaks, -1), sep = ".")
    tally_names <- c(tally_names, break_names)
  }
  names(tally) <- tally_names

  if (!keep_ref_rows) {
    variant_rows <- !is.na(tally$alt)
    if (!all(variant_rows))
      tally <- lapply(tally, `[`, variant_rows)
  }
  
  meta_names <- setdiff(tally_names,
                        c("seqnames", "pos", "ref", "alt", "high.quality",
                          "high.quality.ref", "high.quality.total"))
  genome <- genome(x)
  indel <- nchar(tally$ref) == 0L | nchar(tally$alt) == 0L
  metacols <- DataFrame(tally[meta_names])
  mcols(metacols) <- variantSummaryColumnDescriptions(read_pos_breaks)

  samecols = tally[["ref"]] != tally[["alt"]]
  gr <- with(tally,
             VRanges(seqnames,
                     IRanges(pos,
                             width = ifelse(nchar(alt) == 0, nchar(ref), 1L)),
                     ref, alt,
                     ifelse(indel, raw.count.total, high.quality.total),
                     ifelse(indel, raw.count.ref, high.quality.ref),
                     ifelse(indel, raw.count, high.quality),
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
    stopifnot(all(raw.count + raw.count.ref <= raw.count.total, na.rm=TRUE))
    stopifnot(all(altDepth(x) <= raw.count, na.rm=TRUE))
    stopifnot(all(refDepth(x) <= raw.count.ref, na.rm=TRUE))
    stopifnot(all(count.plus + count.minus == raw.count, na.rm=TRUE))
    stopifnot(all(count.plus.ref + count.minus.ref == raw.count.ref))
  })
}

normalizeIndelAlleles <- function(x, genome) {
  is.indel <- nchar(ref(x)) == 0L | nchar(alt(x)) == 0L
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
                         which = NULL, read_pos_breaks = NULL,
                         high_base_quality = 0L, desired_read_group = NULL,
                         alloclength = 200000L,
                         minimum_mapq = 0L, good_unique_mapq = 35L,
                         maximum_nhits = 1000000L,
                         concordant_only = FALSE, unique_only = FALSE,
                         primary_only = FALSE, ignore_duplicates = FALSE,
                         min_depth = 0L, variant_strand = 0L,
                         ignore_query_Ns = FALSE,
                         indels = FALSE,
                         blocksize = 1000L, verbosep = FALSE,
                         include_soft_clips = 0L,
                         exon_iit = NULL)
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
  if (!is.null(read_pos_breaks)) {
    read_pos_breaks <- as.integer(read_pos_breaks)
    if (any(is.na(read_pos_breaks)))
      stop("'read_pos_breaks' should not contain missing values")
    if (length(read_pos_breaks) < 2)
      stop("'read_pos_breaks' needs at least two elements to define a bin")
    if (is.unsorted(read_pos_breaks))
      stop("'read_pos_breaks' must be sorted")
  }
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
        normArgTRUEorFALSE(ignore_query_Ns),
        normArgTRUEorFALSE(indels),
        normArgSingleInteger(blocksize),
        normArgTRUEorFALSE(verbosep),
        normArgSingleInteger(include_soft_clips),
        normArgSingleCharacterOrNULL(exon_iit))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Column metadata
###

variantSummaryColumnDescriptions <- function(read_pos_breaks) {
  desc <- c(
    n.read.pos = "Number of unique read positions for the ALT",
    n.read.pos.ref = "Number of unique read positions for the REF",
    raw.count = "Raw ALT count",
    raw.count.ref = "Raw REF count",
    raw.count.total = "Raw total count",
    mean.quality = "Average ALT base quality",
    mean.quality.ref = "Average REF base quality",
    count.plus = "Raw positive strand ALT count",
    count.plus.ref = "Raw positive strand REF count",
    count.minus = "Raw negative strand ALT count",
    count.minus.ref = "Raw negative strand REF count",
    read.pos.mean = "Average read position for the ALT",
    read.pos.mean.ref = "Average read position for the ALT",
    read.pos.var = "Variance in read position for the ALT",
    read.pos.var.ref = "Variance in read position for the REF",
    mdfne = "Median distance from nearest end of read for the ALT",
    mdfne.ref = "Median distance from nearest end of read for the REF",
    codon.dir = "Direction of transcription for the codon. 0=positive, 1=negative, 2=NA (not a codon)" )
  if (length(read_pos_breaks) > 0L) {
    break_desc <- paste0("Raw ALT count in read position range [",
                         head(read_pos_breaks, -1), ",",
                         tail(read_pos_breaks, -1), ")")
    desc <- c(desc, break_desc)
  }
  DataFrame(Description = desc)
}
