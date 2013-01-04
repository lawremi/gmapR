### =========================================================================
### BamTallyParam class
### -------------------------------------------------------------------------
###
### Collects parameters for the bam_tally command
###

setClass("BamTallyParam",
         representation(genome = "GmapGenome",
                        which = "RangesList",
                        cycle_breaks = "integerORNULL",
                        high_base_quality = "integer",
                        minimum_mapq = "integer",
                        concordant_only = "logical",
                        unique_only = "logical",
                        primary_only = "logical",
                        min_depth = "integer",
                        variant_strand = "integer",
                        ignore_query_Ns = "logical",
                        indels = "logical"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

normArgWhich <- function(x) {
  if (is(x, "GenomicRanges"))
    x <- split(ranges(x), seqnames(x))
  else if (!is(x, "RangesList"))
    stop("'which' must be a GenomicRanges or RangesList")
  x
}

BamTallyParam <- function(genome, which = RangesList(),
                          cycle_breaks = NULL,
                          high_base_quality = 0L,
                          minimum_mapq = 0L,
                          concordant_only = FALSE, unique_only = FALSE,
                          primary_only = FALSE,
                          min_depth = 0L, variant_strand = 0L,
                          ignore_query_Ns = FALSE,
                          indels = FALSE)
{
  args <- names(formals(sys.function()))
  params <- mget(args, environment())
  params$genome <- as(genome, "GmapGenome")
  params$which <- normArgWhich(which)
  integer_params <- c("high_base_quality", "minimum_mapq", "min_depth",
                      "variant_strand")
  params[integer_params] <- lapply(params[integer_params], as.integer)
  do.call(new, c("BamTallyParam", params))  
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

setAs("BamTallyParam", "list", function(from) {
  sapply(slotNames(from), slot, object = from, simplify = FALSE)
})

setMethod("as.list", "BamTallyParam", function(x) as(x, "list"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

setMethod("show", "BamTallyParam", function(object) {
  cat("A BamTallyParam object\n")
  cat("genome:", genome(object@genome), "\n")
  which <- if (length(object@which) > 0) {
    paste0(space(object@which), ":", unlist(start(object@which)), "-",
           unlist(end(object@which)))
  } else "whole genome"
  cat(IRanges:::labeledLine("which", which))
  otherSlotNames <- setdiff(slotNames(object), c("genome", "which"))
  slots <- sapply(otherSlotNames, slot, object = object, simplify = FALSE)
  cat(mapply(IRanges:::labeledLine, names(slots), slots, count = FALSE),
      sep = "")
})
