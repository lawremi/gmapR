### =========================================================================
### GsnapOutput class
### -------------------------------------------------------------------------
###
### Stores a path to a gSNAP output directory. This is based on the
### assumption that the alignments from a sample will land in a
### unique directory.
###

setClassUnion("POSIXlt_OR_NULL", c("POSIXlt", "NULL"))

.valid_GsnapOutput <- function(object) {
  if (length(samPaths(object)) == 0L && length(bamPaths(object)) == 0L)
    paste0("No GSNAP output at '", object@path, "'")
}

setClassUnion("GmapAlignerParam_OR_NULL", c("GmapAlignerParam", "NULL"))

setClassUnion("GsnapParam_OR_NULL", c("GsnapParam", "NULL"))
setIs("GsnapParam_OR_NULL", "GmapAlignerParam_OR_NULL")

setClass("GmapAlignerOutput",
         representation(path = "character",
                        param = "GmapAlignerParam_OR_NULL",
                        version = "POSIXlt_OR_NULL"))
         
setClass("GsnapOutput",
         representation(param = "GsnapParam_OR_NULL"),
         contains="GmapAlignerOutput",
         validity = .valid_GsnapOutput)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

setMethod("path", "GmapAlignerOutput", function(object) object@path)

setGeneric("alignmentCategories",
           function(x) standardGeneric("alignmentCategories"))

setMethod("alignmentCategories", "GsnapOutput", function(x) {
              c("concordant_uniq", "concordant_mult",
                "concordant_transloc", "paired_mult",
                "paired_uniq_inv", "paired_uniq_long",
                "paired_uniq_scr", "unpaired_uniq", "unpaired_mult",
                "unpaired_transloc",  "halfmapping_uniq",
                "halfmapping_mult", "halfmapping_transloc",
                "nomapping", "concordant_circular",
                "halfmapping_circular", "paired_uniq_circular",
                "unpaired_circular")
          })

### Q: How to handle the single BAM output case? We would still want a
### formal object, with a version, parameters, etc. If 'path' is a
### directory, we would assume multiple (split) BAMs in that
### directory. Otherwise, it is a single BAM. We provide the bamPath()
### convenience function that will always yield a single, sensible BAM
### for the output. In the multiple output case, this is the
### concordant_uniq file.

is_dir <- function(x) {
  file.exists(path(x)) && file.info(path(x))[,"isdir"]
}

setGeneric("bamPath", function(x) standardGeneric("bamPath"))

setMethod("bamPath", "GmapAlignerOutput", function(x) {
              paths <- bamPaths(x)
              if (length(paths) > 1L)
                  paths <- paths[alignmentCategories(x)[1L]]
              paths
          })

## file_ext does not allow '_'
file_ext2 <- function(x) {
  gsub(".*\\.", "", x)
}

paths <- function(x) {
    if (!is_dir(x)) {
        return(path(x))
    }
    paths <- list_files_with_exts(path(x), alignmentCategories(x),
                                  full.names = TRUE)
    names(paths) <- file_ext2(paths)
    paths[alignmentCategories(x)]
}

samPaths <- function(x) {
  if (!is_dir(x)) {
    if (file_ext(path(x)) == "sam")
      return(path(x))
    return(character())
  }
  paths(x)
}

setMethod("bamPaths", "GmapAlignerOutput", function(x) {
  if (!is_dir(x)) {
    if (file_ext(path(x)) == "bam")
      return(path(x))
    return(character())
  }
  paths <- list_files_with_exts(path(x), "bam", full.names = TRUE)
  names(paths) <- file_ext2(file_path_sans_ext(paths))
  paths <- paths[names(paths) %in% alignmentCategories(x)]
  paths
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

newGmapAlignerOutput <- function(Class, path, param = NULL, version = NULL) {
    if (!is.character(path) || any(is.na(path)))
        stop("'path' must be a character vector without any NA's")
    if (!is(version, "POSIXlt_OR_NULL")) {
        if (is.character(version))
            version <- parseGsnapVersion(version)
        else stop("'version' must be a version string, POSIXlt object or NULL")
    }
    path <- file_path_as_absolute(path)
    new(Class, path = path, version = version, param = param)
}

GsnapOutput <- function(path, param = NULL, version = NULL) {
  newGmapAlignerOutput("GsnapOutput", path, param, version)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

setAs("GmapAlignerOutput", "BamFileList", function(from) {
  BamFileList(bamPaths(from))
})

setAs("GmapAlignerOutput", "BamFile", function(from) {
  BamFile(bamPath(from))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

setMethod("merge", c("GmapAlignerOutput", "GmapAlignerOutput"), function(x, y) {
### TODO: get Cory's lane merging code for BAMs in here
})

setGeneric("consolidate", function(x, ...) standardGeneric("consolidate"))

setMethod("consolidate", "GmapAlignerOutput", function(x) {
### TODO: this should produce the pipeline's analyzed BAM
  x
})

##converts all gsnap SAM files to BAM files and creates the .bai index files
setMethod("asBam", "GmapAlignerOutput",
          function(file) {
            ## the input is not a file. It is a GsnapOutput,
            ## but param name is enforced
            gsp <- file 

            ##files other than those produced by gsnap maybe be in the
            ##output directory. Only take those produced by gsnap.
            samFiles <- samPaths(gsp)
            bamFiles <- mapply(asBam, file = samFiles,
                               dest = samFiles,
                               MoreArgs = list(overwrite = TRUE))
            unlink(samFiles)

            if (!is_dir(file))
              file@path <- bamFiles
            
            file
})

setMethod("import", c("GsnapOutput", "missing", "missing"),
          function(con, format, text, ...) {
              import(as(con, "BamFile"), ...)
          })

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

setClass("GmapAlignerOutputList",
         prototype = prototype(elementType = "GmapAlignerOutput"),
         contains = "List")

setClass("GsnapOutputList",
         prototype = prototype(elementType = "GsnapOutput"),
         contains = "GmapAlignerOutputList")

setClass("SimpleGsnapOutputList",
         contains = c("GsnapOutputList", "SimpleList"))

GsnapOutputList <- function(...) {
  args <- list(...)
  if (length(args) == 1 && is.list(args[[1]])) 
    args <- args[[1]]
  S4Vectors:::new_SimpleList_from_list("SimpleGsnapOutputList", args)
}

setAs("GmapAlignerOutputList", "BamFileList", function(from) {
  BamFileList(lapply(from, as, "BamFile"))
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

setMethod("show", "GmapAlignerOutput", function(object) {
  cat("A ", class(object), " Object\n", "path: ", path(object), "\n", sep = "")
})
