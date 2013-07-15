### =========================================================================
### bam_tally program
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Raw tally result
###

setClass("TallyIIT", representation(ptr = "externalptr", genome = "GmapGenome"))

TallyIIT <- function(ptr, genome) {
  new("TallyIIT", ptr = ptr, genome = genome)
}

setMethod("genome", "TallyIIT", function(x) x@genome)

setMethod("show", "TallyIIT", function(object) {
  cat("Tally IIT object for genome '", genome(genome(object)), "'\n", sep = "")
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
            TallyIIT(do.call(.bam_tally_C, c(list(x), param_list)), genome)
          })

variantSummary <- function(x, read_pos_breaks = NULL, high_base_quality = 0L)
{
  tally <- .Call(R_tally_iit_parse, x@ptr,
                 read_pos_breaks,
                 normArgSingleInteger(high_base_quality),
                 NULL)

  tally_names <- c("seqnames", "pos", "ref", "alt", "n.read.pos",
                   "n.read.pos.ref", "raw.count", "raw.count.ref",
                   "raw.count.total",
                   "high.quality", "high.quality.ref",
                   "high.quality.total", "mean.quality",
                   "mean.quality.ref",
                   "count.pos", "count.pos.ref",
                   "count.neg", "count.neg.ref",
                   "read.pos.mean", "read.pos.mean.ref",
                   "read.pos.var", "read.pos.var.ref")
  if (length(read_pos_breaks) > 0L) {
    read_pos_breaks <- as.integer(read_pos_breaks)
    break_names <- paste("readPosCount", head(read_pos_breaks, -1),
                         tail(read_pos_breaks, -1), sep = ".")
    tally_names <- c(tally_names, break_names)
  }
  names(tally) <- tally_names

  variant_rows <- !is.na(tally$alt)
  if (!all(variant_rows))
    tally <- lapply(tally, `[`, variant_rows)
  
  meta_names <- setdiff(tally_names,
                        c("seqnames", "pos", "ref", "alt", "high.quality",
                          "high.quality.total"))
  genome <- genome(x)
  indel <- nchar(tally$ref) == 0L | nchar(tally$alt) == 0L
  gr <- with(tally,
             VRanges(seqnames,
                     IRanges(pos,
                             width = ifelse(nchar(alt) == 0, nchar(ref), 1L)),
                     ref, alt,
                     ifelse(indel, raw.count.total, high.quality.total),
                     ifelse(indel, raw.count.ref, high.quality.ref),
                     ifelse(indel, raw.count, high.quality),
                     DataFrame(tally[meta_names]),
                     seqlengths = seqlengths(genome)))
  seqinfo(gr) <- seqinfo(genome)
  gr <- normalizeIndelAlleles(gr, genome)
  gr
}

normalizeIndelAlleles <- function(x, genome) {
  is.indel <- nchar(ref(x)) == 0L | nchar(alt(x)) == 0L
  if (any(is.indel)) {
    indels <- x[is.indel]
    indels <- shift(indels, -1)
    anchor <- getSeq(genome, indels)
    ref(indels) <- paste0(anchor, ref(indels))
    alt(indels) <- paste0(anchor, alt(indels))
    x[is.indel] <- indels
  }
  x
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level wrappers
###

normArgSingleInteger <- function(x) {
  name <- deparse(substitute(x))
  x <- as.integer(x)
  if (!IRanges:::isSingleInteger(x))
    stop("'", name, "' should be a single, non-NA integer")
  x
}
normArgTRUEorFALSE <- function(x) {
  name <- deparse(substitute(x))    
  if (!isTRUEorFALSE(x))
    stop("'", name, "' should be TRUE or FALSE")
  x
}

.bam_tally_C <- function(bamreader, genome_dir = NULL, db = NULL,
                         which = NULL, read_pos_breaks = NULL,
                         high_base_quality = 0L, alloclength = 200000L,
                         minimum_mapq = 0L, good_unique_mapq = 35L,
                         maximum_nhits = 1000000L,
                         concordant_only = FALSE, unique_only = FALSE,
                         primary_only = FALSE, ignore_duplicates = FALSE,
                         min_depth = 0L, variant_strand = 0L,
                         ignore_query_Ns = FALSE,
                         indels = FALSE,
                         blocksize = 1000L, verbosep = FALSE)
{
  if (!is(bamreader, "GmapBamReader"))
    stop("'bamreader' must be a GmapBamReader")
  if (length(which) > 0L) {
    which <- list(as.character(seqnames(which)), start(which), end(which))
  } else which <- NULL
  if (!is.null(genome_dir) && !IRanges:::isSingleString(genome_dir))
    stop("'genome_dir' must be NULL or a single, non-NA string")
  if (!is.null(db) && !IRanges:::isSingleString(db))
    stop("'db' must be NULL or a single, non-NA string")
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
        normArgTRUEorFALSE(verbosep))
}
