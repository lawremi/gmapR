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

## should contain any category ever known to be used by gsnap
.ALIGNMENT_CATEGORIES <- c("concordant_uniq", "concordant_mult",
                           "concordant_transloc", "paired_mult",
                           "paired_uniq_inv", "paired_uniq_long",
                           "paired_uniq_scr", "unpaired_uniq", "unpaired_mult",
                           "unpaired_transloc",  "halfmapping_uniq",
                           "halfmapping_mult", "halfmapping_transloc",
                           "nomapping")

setMethod("bamPaths", "GsnapOutput", function(x) {
  paths <- list_files_with_exts(path(x), "bam", full.names = TRUE)
  names(paths) <- sub(".*\\.([^.]*)\\.bam$", "\\1", paths)
  paths <- paths[names(paths) %in% .ALIGNMENT_CATEGORIES]
  paths
})

##' Constructor for GsnapOutput class, which represents the output
##' (directory) of a gsnap run.
##'
##' 
##' @title GsnapOutput constructor
##' @param path Path to the gsnap output directory
##' @return New instance
##' @author Michael Lawrence
##' @export
GsnapOutput <- function(path) {
  if (!is.character(path) || any(is.na(path)))
    stop("'path' must be a character vector without any NA's")
  new("GsnapOutput", path = path)
}
