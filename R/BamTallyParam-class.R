### =========================================================================
### BamTallyParam class
### -------------------------------------------------------------------------
###
### Collects parameters for the bam_tally command
###

setClass("BamTallyParam",
         representation(genome = "GmapGenome",
                        which = "GenomicRanges",
                        desired_read_group = "characterORNULL",
                        minimum_mapq = "integer",
                        concordant_only = "logical",
                        unique_only = "logical",
                        primary_only = "logical",
                        ignore_duplicates = "logical",
                        min_depth = "integer",
                        variant_strand = "integer",
                        ignore_query_Ns = "logical",
                        indels = "logical"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

normArgWhich <- function(x, genome) {
  if (!is(x, "GenomicRanges"))
    stop("'which' must be a GenomicRanges")
  si <- seqinfo(genome)
  seqinfo(x, new2old = match(seqlevels(si), seqlevels(x))) <-
    merge(si, seqinfo(x))
  x
}

BamTallyParam <- function(genome, which = GRanges(),
                          desired_read_group = NULL, minimum_mapq = 0L,
                          concordant_only = FALSE, unique_only = FALSE,
                          primary_only = FALSE, ignore_duplicates = FALSE,
                          min_depth = 0L, variant_strand = 0L,
                          ignore_query_Ns = FALSE,
                          indels = FALSE)
{
  if (!is.null(desired_read_group) && !isSingleString(desired_read_group))
    stop("'desired_read_group' must be NULL or a single, non-NA string")
  if (!isSingleNumber(minimum_mapq) || minimum_mapq < 0)
    stop("minimum_mapq must be a single, non-negative, non-NA number")
  if (!isTRUEorFALSE(concordant_only))
    stop("concordant_only must be TRUE or FALSE")
  if (!isTRUEorFALSE(unique_only))
    stop("unique_only must be TRUE or FALSE")
  if (!isTRUEorFALSE(primary_only))
    stop("primary_only must be TRUE or FALSE")
  if (!isTRUEorFALSE(ignore_duplicates))
    stop("ignore_duplicates must be TRUE or FALSE")
  if (!isSingleNumber(min_depth) || min_depth < 0)
    stop("min_depth must be a single, non-negative, non-NA number")
  if (!variant_strand %in% c(0, 1, 2))
    stop("variant_strand must be one of 0, 1, or 2")
  if (!isTRUEorFALSE(ignore_query_Ns))
    stop("ignore_query_Ns must be TRUE or FALSE")
  if (!isTRUEorFALSE(indels))
    stop("indels must be TRUE or FALSE")
  args <- names(formals(sys.function()))
  params <- mget(args, environment())
  params$genome <- as(genome, "GmapGenome")
  params$which <- normArgWhich(which, params$genome)
  integer_params <- c("minimum_mapq", "min_depth", "variant_strand")
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

showSlot <- function(name, value, ...) {
  IRanges:::labeledLine(name, showAsCell(value), ...)
}

showSlots <- function(object, exclude = character(), ...) {
  snames <- setdiff(slotNames(object), exclude)
  slots <- sapply(snames, slot, object = object, simplify = FALSE)
  mapply(showSlot, names(slots), slots, MoreArgs = list(...))
}

setMethod("show", "BamTallyParam", function(object) {
  cat("A", class(object), "object\n", sep = " ")
  cat(showSlots(object, count = FALSE), sep = "")
})
