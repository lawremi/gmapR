### =========================================================================
### GmapBamReader class
### -------------------------------------------------------------------------
###
### The C-level BAM file iterator used by BAM processors in gmap/gstruct.
###
### This is very similar to BamFile from Rsamtools, but it iterates
### line-by-line, not by chunk.
###

setClass("GmapBamReader",
         representation(.extptr = "externalptr", path = "character"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

setMethod("path", "GmapBamReader", function(object) object@path)

## setReplaceMethod("bamWhich", c("GmapBamReader", "ANY"),
##                  function(object, value) {
##                    bamWhich(object) <- as(value, "RangesList")
##                    object
##                  })

## setReplaceMethod("bamWhich", c("GmapBamReader", "RangesList"),
##                  function(object, value) {
##                    if (length(value) != 1L || length(value[[1]]) != 1L)
##                      stop("'value' must be a RangesList with a single, ",
##                           "length-one element")
##                    .Call(R_Bamread_limit_region, object, names(value),
##                          start(value[[1]]), end(value[[1]]))
##                  })

## setReplaceMethod("bamWhich", c("GmapBamReader", "NULL"),
##                  function(object, value) {
##                    .Call(R_Bamread_unlimit_region, object)
##                  })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

GmapBamReader <- function(path, which = NULL) {
  if (is(path, "BamFile"))
    path <- path(path)
  if (!IRanges:::isSingleString(path))
    stop("'path' must be a single, non-NA string")
  path <- path.expand(path)
  obj <- new("GmapBamReader", .extptr = .Call(R_Bamread_new, path), path = path)
### TODO:
  ##  bamWhich(obj) <- which
  obj
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

setAs("BamFile", "GmapBamReader", function(from) {
  GmapBamReader(from)
})

setAs("GmapBamReader", "BamFile", function(from) {
  BamFile(path(from))
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

setMethod("show", "GmapBamReader", function(object) {
  cat("GmapBamReader object\npath:", path(object), "\n")
})
