### =========================================================================
### GsnapParam class
### -------------------------------------------------------------------------
###
### High-level interface to gsnap. As a complex operation, we need a
### formal parameter object. It should only formally represent the
### most commonly used parameters. The rest fall into the 'extra'
### list.
###

setClassUnion("integer_OR_NULL", c("integer", "NULL"))

setClass("GsnapParam",
         representation(max_mismatches = "integer_OR_NULL",
                        suboptimal_levels = "integer",
                        novelsplicing = "logical",
                        splicing = "character_OR_NULL",
                        terminal_threshold = "integer",
                        gmap_mode = "character",
                        clip_overlap = "logical"),
         contains="GmapAlignerParam")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

gsnap_part <- function(x) {
  x@part
}
gsnap_batch <- function(x) {
  x@batch
}
gsnap_max_mismatches <- function(x) {
  x@max_mismatches
}
gsnap_suboptimal_levels <- function(x) {
  x@suboptimal_levels
}
gsnap_use_snps <- function(x) {
  x@use_snps
}
gsnap_snpsdir <- function(x) {
  x@snpsdir
}
gsnap_mode <- function(x) {
  x@mode
}
gsnap_nthreads <- function(x) {
  x@nthreads
}
gsnap_novelsplicing <- function(x) {
  x@novelsplicing
}
gsnap_splicing <- function(x) {
  x@splicing
}
gsnap_npaths <- function(x) {
  x@npaths
}
gsnap_quiet_if_excessive <- function(x) {
  x@quiet_if_excessive
}
gsnap_nofails <- function(x) {
  x@nofails
}
gsnap_split_output <- function(x) {
  x@split_output
}
gsnap_extra <- function(x) {
  x@extra
}

`gsnap_part<-` <- function(x, value) {
  x@part <- value
  x
}
`gsnap_batch<-` <- function(x, value) {
  x@batch <- value
  x
}
`gsnap_max_mismatches<-` <- function(x, value) {
  x@max_mismatches <- value
  x
}
`gsnap_suboptimal_levels<-` <- function(x, value) {
  x@suboptimal_levels <- value
  x
}
`gsnap_use_snps<-` <- function(x, value) {
  x@use_snps <- value
  x
}
`gsnap_snpsdir<-` <- function(x, value) {
  x@snpsdir <- value
  x
}
`gsnap_mode<-` <- function(x, value) {
  x@mode <- value
  x
}
`gsnap_nthreads<-` <- function(x, value) {
  x@nthreads <- value
  x
}
`gsnap_novelsplicing<-` <- function(x, value) {
  x@novelsplicing <- value
  x
}
`gsnap_splicing<-` <- function(x, value) {
  x@splicing <- value
  x
}
`gsnap_npaths<-` <- function(x, value) {
  x@npaths <- value
  x
}
`gsnap_quiet_if_excessive<-` <- function(x, value) {
  x@quiet_if_excessive <- value
  x
}
`gsnap_nofails<-` <- function(x, value) {
  x@nofails <- value
  x
}
`gsnap_split_output<-` <- function(x, value) {
  x@split_output <- value
  x
}
`gsnap_extra<-` <- function(x, value) {
  x@extra <- value
  x
}

## etc..

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

newGmapAlignerParam <- function(Class, genome, snps) {
    if (missing(genome))
        stop("The 'genome' must be specified (and coercible to GmapGenome)")
    args <- formals(sys.function(sys.parent(1L)))
    params <- mget(names(args), parent.frame())
    params$unique_only <- NULL
    paramClasses <- getSlots(Class)
    paramClasses <- paramClasses[setdiff(names(paramClasses),
                                         c("extra", "snps"))]
    params <- mapply(as, params[names(paramClasses)], paramClasses,
                     SIMPLIFY = FALSE)
    if (!is.null(snps)) {
        if (!is(snps, "GmapSnps")) {
            snps <- GmapSnps(snps, genome)
        }
        params$snps <- snps
    }
    params$extra <- evalq(list(...), parent.frame())
    do.call(new, c(Class, params))
}

GsnapParam <- function(genome, unique_only = FALSE, molecule = c("RNA", "DNA"),
                       max_mismatches = NULL,
                       suboptimal_levels = 0L, mode = "standard",
                       snps = NULL,
                       npaths = if (unique_only) 1L else 100L,
                       quiet_if_excessive = unique_only, nofails = unique_only,
                       split_output = !unique_only,
                       novelsplicing = FALSE, splicing = NULL, 
                       nthreads = 1L, part = NULL, batch = "2",
                       terminal_threshold =
                           if (molecule == "DNA") 1000L else 2L,
                       gmap_mode = if (molecule == "DNA") "none"
                                   else "pairsearch,terminal,improve",
                       clip_overlap = FALSE, ...)
{
  molecule <- match.arg(molecule)
  newGmapAlignerParam("GsnapParam", genome, snps)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

GmapAlignerParam_asList <- function(from) {
  to <- lapply(slotNames(from), slot, object = from)
  names(to) <- slotNames(from)
  to$split_output <- if (to$split_output) tolower(sub("Param", "", class(from)))
  to$db <- genome(to$genome)
  to$dir <- path(directory(to$genome))
  to$genome <- NULL
  to$use_snps <- name(to$snps)
  to$snpsdir <- path(directory(to$snps))
  extras <- to$extra
  to <- c(to, extras)
  to$extra <- NULL  
  to
}

setAs("GsnapParam", "list", function(from) {
          to <- GmapAlignerParam_asList(from)
          to$novelsplicing <- as.integer(to$novelsplicing)
          to <- rename(to, splicing = "use_splicing")
          to
      })

as.list.GmapAlignerParam <- function(x, ...) as(x, "list")

setAs("ANY", "character_OR_NULL", function(from) {
  if (is.null(from))
    NULL
  else as.character(from)
})
setAs("ANY", "integer_OR_NULL", function(from) {
  if (is.null(from))
    NULL
  else as.integer(from)
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

setMethod("show", "GmapAlignerParam", function(object) {
  slots <- lapply(slotNames(object), slot, object = object)
  names(slots) <- slotNames(object)
  slots$genome <- paste0(slots$genome@name,
                         " (", path(directory(slots$genome)), ")")
  if (!is.null(slots$snps)) {
    slots$snps <- paste0(name(slots$snps),
                         " (", path(directory(slots$snps)), ")")
  }
  cat("A ", class(object), " object\n",
      paste0(names(slots), ": ", slots, collapse = "\n"), "\n", sep = "")
})
