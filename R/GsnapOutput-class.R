### =========================================================================
### GsnapOutput class
### -------------------------------------------------------------------------
###
### Stores a path to a gSNAP output directory. This is based on the
### assumption that the alignments from a sample will land in a
### unique directory.
###

setOldClass("POSIXlt")
setClassUnion("POSIXltORNULL", c("POSIXlt", "NULL"))

setClassUnion("GsnapParamORNULL", c("GsnapParam", "NULL"))

setClass("GsnapOutput",
         representation(path = "character",
                        param = "GsnapParamORNULL",
                        version = "POSIXltORNULL"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

setMethod("path", "GsnapOutput", function(object) object@path)

## should contain any category ever known to be used by gsnap
.ALIGNMENT_CATEGORIES <- c("concordant_uniq", "concordant_mult",
                           "concordant_transloc", "paired_mult",
                           "paired_uniq_inv", "paired_uniq_long",
                           "paired_uniq_scr", "unpaired_uniq", "unpaired_mult",
                           "unpaired_transloc",  "halfmapping_uniq",
                           "halfmapping_mult", "halfmapping_transloc",
                           "nomapping")

### Q: How to handle the single BAM output case? We would still want a
### formal object, with a version, parameters, etc. If 'path' is a
### directory, we would assume multiple (split) BAMs in that
### directory. Otherwise, it is a single BAM. We provide the bamPath()
### convenience function that will always yield a single, sensible BAM
### for the output. In the multiple output case, this is the
### concordant_uniq file.

is_dir <- function(x) file.info(path(x))[,"isdir"]

setGeneric("bamPath", function(x) standardGeneric("bamPath"))

setMethod("bamPath", "GsnapOutput", function(x) {
  paths <- bamPaths(x)
  if (length(paths) > 1) # maybe this should be the analysis BAM if it exists?
    paths <- paths["concordant_uniq"]
  paths
})

setMethod("bamPaths", "GsnapOutput", function(x) {
  if (!is_dir(x))
    return(path(x))
  paths <- list_files_with_exts(path(x), "bam", full.names = TRUE)
  names(paths) <- sub(".*\\.([^.]*)\\.bam$", "\\1", paths)
  paths <- paths[names(paths) %in% .ALIGNMENT_CATEGORIES]
  paths
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

GsnapOutput <- function(path, param = NULL, version = NULL) {
  if (!is.character(path) || any(is.na(path)))
    stop("'path' must be a character vector without any NA's")
  if (!is(version, "POSIXltORNULL")) {
    if (is.character(version))
      version <- parseGsnapVersion(version)
    else stop("'version' must be a version string, POSIXlt object or NULL")
  }
  if (!is(param, "GsnapParamORNULL"))
    stop("'param' must be a GsnapParam object or NULL")
  new("GsnapOutput", path = path, version = version, param = param)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

setAs("GsnapOutput", "BamFileList", function(from) {
  BamFileList(bamPaths(from))
})

setAs("GsnapOutput", "BamFile", function(from) {
  BamFile(bamPath(from))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

setMethod("merge", c("GsnapOutput", "GsnapOutput"), function(x, y) {
### TODO: get Cory's lane merging code for BAMs in here
})

setGeneric("consolidate", function(x, ...) standardGeneric("consolidate"))

setMethod("consolidate", "GsnapOutput", function(x) {
### TODO: this should produce the pipeline's analyzed BAM
### TODO: the 'gsnap' function should probably consolidate by default.
})

##converts all gsnap SAM files to BAM files and creates the .bai index files
setMethod("asBam", "GsnapOutput",
          function(file) {
            gsp <- file ##the input is not a file. It is a GsnapOutput, but param name is enforced

            ##files other than those produced by gsnap maybe be in the
            ##output directory. Only take those produced by gsnap
            samFiles <- dir(path(gsp), full.names=TRUE)
            isFromGsnap <- sapply(samFiles,
                                  function(fileName) {
                                    any(sapply(.GSNAP_OUTPUT_FILE_STRINGS,
                                           function(pat) grepl(pat, fileName)))
                                  })
            samFiles <- samFiles[isFromGsnap]
            sapply(samFiles,
                   function(samFile) {
                     ##.bam will the added to the string supplied via
                     ##dest param
                     asBam(samFile, dest=samFile, indexDestination=TRUE)
                   })

            unlink(samFiles)            
})

.GSNAP_OUTPUT_FILE_STRINGS <- c("\\.nomapping$",
                                "\\.unpaired_mult$",
                                "\\.unpaired_transloc$",
                                "\\.unpaired_uniq$")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### List class
###

### Q: How to handle multiple output directories? We could vectorize
### GsnapOutput, but is it not simpler to have a GsnapOutputList
### object? If this object is vectorized, then every method would need
### to handle the multiple element case, and accessor methods like
### bamPaths would need to be CharacterList instead of simple
### character vectors. Would we vectorize the version or not? It makes
### sense to have an object that represents a single run of gsnap.

setClass("GsnapOutputList",
         prototype = prototype(elementType = "GsnapOutput"),
         contains = "List")

setClass("SimpleGsnapOutputList",
         contains = c("GsnapOutputList", "SimpleList"))
