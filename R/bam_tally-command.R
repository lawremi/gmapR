### =========================================================================
### bam_tally program
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### High-level wrapper
###

setGeneric("bam_tally",
           function(x, genome, param = BamTallyParam(), ...)
           standardGeneric("bam_tally"),
           signature = "x")

setMethod("bam_tally", "BamFile",
          function(x, genome, param, ...)
          {
            x <- GmapBamReader(x)
            callGeneric()
          })

setMethod("bam_tally", "character",
          function(x, genome, param, ...)
          {
            x <- BamFile(x)
            callGeneric()
          })

setMethod("bam_tally", "GmapBamReader",
          function(x, genome, param, ...)
          {
            param_list <- as.list(param)
            args <- list(...)
            param_list[names(args)] <- args
            param_list$db <- genome(genome)
            param_list$genome_dir <- path(directory(genome))
            tally <- do.call(.bam_tally_C, c(list(x), param_list))
            tally_names <- c("seqnames", "pos", "ref", "alt", "ncycles",
                             "ncycles.ref", "count", "count.ref",
                             "count.total", "high.quality", "high.quality.ref",
                             "mean.quality", "mean.quality.ref",
                             "count.pos", "count.pos.ref",
                             "count.neg", "count.neg.ref")
            cycle_breaks <- param_list$cycle_breaks
            if (!is.null(cycle_breaks)) {
              cycle_breaks <- as.integer(cycle_breaks)
              break_names <- paste("cycleCount", head(cycle_breaks, -1),
                                    tail(cycle_breaks, -1), sep = ".")
              tally_names <- c(tally_names, break_names)
            }
            names(tally) <- tally_names
            if (param_list$variant_strand > 0) {
              variant_rows <- !is.na(tally$alt)
              tally <- lapply(tally, `[`, variant_rows)
            }
            gr <- GRanges(tally$seqnames,
                          IRanges(tally$pos,
                                  width = rep.int(1L, length(tally$pos))),
                          strand = Rle("+", length(tally$pos)),
                          location = paste(tally$seqnames, tally$pos,
                            sep = ":"),
                          DataFrame(tail(tally, -2)),
                          seqlengths = seqlengths(genome))
            seqinfo(gr) <- seqinfo(genome)
            gr
          })

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
                         which = NULL, cycle_breaks = NULL,
                         high_quality_cutoff = 0L, alloclength = 200000L,
                         minimum_mapq = 0L, good_unique_mapq = 35L,
                         maximum_nhits = 1000000L,
                         concordant_only = FALSE, unique_only = FALSE,
                         primary_only = FALSE,
                         min_depth = 0L, variant_strand = 0L,
                         ignore_query_Ns = FALSE,
                         indels = FALSE,
                         blocksize = 1000L, verbosep = FALSE)
{
  if (!is(bamreader, "GmapBamReader"))
    stop("'bamreader' must be a GmapBamReader")
  if (!is.null(which)) {
    which <- as(which, "RangesList")
    spaceIsNULL <- ifelse(is.null(space(which)), TRUE, FALSE)
    which <- list(as.character(space(which)),
                  unlist(start(which), use.names = FALSE),
                  unlist(end(which), use.names = FALSE))
    if (spaceIsNULL) which <- NULL
  }
  if (!is.null(genome_dir) && !IRanges:::isSingleString(genome_dir))
    stop("'genome_dir' must be NULL or a single, non-NA string")
  if (!is.null(db) && !IRanges:::isSingleString(db))
    stop("'db' must be NULL or a single, non-NA string")
  if (!is.null(cycle_breaks)) {
    cycle_breaks <- as.integer(cycle_breaks)
    if (any(is.na(cycle_breaks)))
      stop("'cycle_breaks' should not contain missing values")
    if (length(cycle_breaks) < 2)
      stop("'cycle_breaks' needs at least two elements to define a bin")
    if (is.unsorted(cycle_breaks))
      stop("'cycle_breaks' must be sorted")
  }
  .Call(R_Bamtally_iit, bamreader@.extptr, genome_dir, db, which,
        cycle_breaks,
        normArgSingleInteger(high_quality_cutoff),
        normArgSingleInteger(alloclength),
        normArgSingleInteger(minimum_mapq),
        normArgSingleInteger(good_unique_mapq),
        normArgSingleInteger(maximum_nhits),
        normArgTRUEorFALSE(concordant_only),
        normArgTRUEorFALSE(unique_only),
        normArgTRUEorFALSE(primary_only),
        normArgSingleInteger(min_depth),
        normArgSingleInteger(variant_strand),
        normArgTRUEorFALSE(ignore_query_Ns),
        normArgTRUEorFALSE(indels),
        normArgSingleInteger(blocksize),
        normArgTRUEorFALSE(verbosep))
}
