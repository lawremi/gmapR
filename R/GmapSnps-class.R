### =========================================================================
### GmapSnps class
### -------------------------------------------------------------------------
###
### A set of SNPs in a GMAP SNP directory
###

setClass("GmapSnps",
         representation(name = "character",
                        directory = "GmapSnpDirectory"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

setGeneric("name", function(x) standardGeneric("name"))

setMethod("name", "GmapSnps", function(x) {
  x@name
})

setMethod("directory", "GmapSnps", function(x) {
  x@directory
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

GmapSnps <- function(snps, directory, name = snps, create = FALSE, ...)
{
  if (!isTRUEorFALSE(create))
    stop("'create' must be TRUE or FALSE")
  if (isSingleString(directory) || is(directory, "GmapGenome"))
    directory <- GmapSnpDirectory(directory, create = create)
  if (!is(directory, "GmapSnpsDirectory"))
    stop("'directory' must be a GmapSnpsDirectory or a path to one")
  if (!isSingleString(name))
    stop("'name' must be a single, non-NA string")
  db <- new("GmapSnps", name = name, directory = directory)
  if (create) {
    if (name %in% names(directory))
      message("NOTE: snps db '", name, "' already exists, not overwriting")
    else snps(db, ...) <- snps
  }
  db
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

setMethod("show", "GmapSnps", function(object) {
  cat("GmapSnps object\nname:", name(object), "\ndirectory:",
      path(directory(object)), "\n")
})
