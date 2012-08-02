### =========================================================================
### GsnapParam class
### -------------------------------------------------------------------------
###
### High-level interface to gsnap. As a complex operation, we need a
### formal parameter object. It should only formally represent the
### most commonly used parameters. The rest fall into the 'extra'
### list.
###

setClassUnion("integerORNULL", c("integer", "NULL"))
setClassUnion("GmapSnpsORNULL", c("GmapSnps", "NULL"))

setClass("GsnapParam",
         representation(genome = "GmapGenome",
                        part = "characterORNULL", # used by parallelized_gsnap
                        batch = "character", # weird "0", "1", ... 
                        max_mismatches = "integerORNULL",
                        suboptimal_levels = "integer",
                        snps = "GmapSnpsORNULL",
                        mode = "character",
                        nthreads = "integer",
                        novelsplicing = "logical",
                        splicing = "characterORNULL",
                        npaths = "integer",
                        quiet_if_excessive = "logical",
                        nofails = "logical", 
                        split_output = "logical",
                        extra = "list"))

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

GsnapParam <- function(genome, unique_only = FALSE,
                       max_mismatches = NULL,
                       suboptimal_levels = 0L, mode = "standard",
                       snps = NULL,
                       npaths = if (unique_only) 1L else 100L,
                       quiet_if_excessive = unique_only, nofails = unique_only,
                       split_output = !unique_only,
                       novelsplicing = 0L, splicing = NULL, 
                       nthreads = 1L, part = NULL, batch = "2", ...) {
  params <- formals(sys.function())
  mc <- as.list(match.call(expand.dots = FALSE))[-1L]
  params[names(mc)] <- mc
  params$unique_only <- NULL
  if (!is.null(snps)) {
    if (!is(snps, "GmapSnps")) {
      params$snps <- GmapSnps(snps, genome)
    }
  }
  extra <- params$...
  ##TODO: this breaks if ... has anything in it
  if (!missing(extra) && !is.null(names(extra))) {
    params$extra <- params$...
  }
  params$... <- NULL
  do.call(new, c("GsnapParam", params))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

setAs("GsnapParam", "list", function(from) {
  to <- lapply(slotNames(from), slot, object = from)
  names(to) <- slotNames(from)
  to$split_output <- if (to$split_output) "gsnap" else NULL
  to$db <- genome(to$genome)
  to$dir <- path(directory(to$genome))
  to$genome <- NULL
  if(!is.null(to$snps)) {
    to$use_snps <- name(to$snps)
  } else {
    to$use_snps <- NULL
  }
  to$snpsdir <- path(directory(to$snps))
  to$snps <- NULL
  to
})

setMethod("as.list", "GsnapParam", function(x) as(x, "list"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

setMethod("show", "GsnapParam", function(object) {
  slots <- lapply(slotNames(object), slot, object = object)
  names(slots) <- slotNames(object)
  slots$genome <- paste0(slots$genome@name,
                         " (", path(directory(slots$genome)), ")")
  if (!is.null(slots$snps)) {
    slots$snps <- paste0(name(slots$snps),
                         " (", path(directory(slots$snps)), ")")
  }
  cat("A GsnapParams object\n",
      paste0(names(slots), ": ", slots, collapse = "\n"), "\n", sep = "")
})
