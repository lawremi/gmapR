### =========================================================================
### GsnapOutput class
### -------------------------------------------------------------------------
###
### Stores one or more paths to gSNAP output directories. This is
### based on the assumption that the alignments from some number of
### samples will each land in a unique directory.
###

setOldClass("numeric_version")

setClass("GsnapOutput",
         representation(path = "character",
                        ##param = "GsnapParam", TODO
                        version = "numeric_version"),
         prototype = list(
           version = numeric_version("0")
           ##param = GsnapParam()
           ))

setMethod("path", "GsnapOutput", function(object) object@path)

setMethod("bamPaths", "GsnapOutput", function(x) {
  paths <- list_files_with_exts(path(x), "bam", full.names = TRUE)
  names(paths) <- sub(".*\\.([^.]*)\\.bam$", "\\1", paths)
  paths
})

GsnapOutput <- function(path) {
  if (!is.character(path) || any(is.na(path)))
    stop("'path' must be a character vector without any NA's")
  new("GsnapOutput", path = path)
}
