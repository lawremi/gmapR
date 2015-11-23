### =========================================================================
### GmapOutput class
### -------------------------------------------------------------------------
###
### Stores a path to a GMAP output directory. This is based on the
### assumption that the alignments from a sample will land in a
### unique directory.
###

setClassUnion("GmapParamORNULL", c("GmapParam", "NULL"))
setIs("GmapParamORNULL", "GmapAlignerParamORNULL")

.valid_GmapOutput <- function(object) {
    if (length(paths(object)) == 0L)
        paste0("No GSNAP output at '", object@path, "'")
}

setClass("GmapOutput",
         representation(path = "character",
                        param = "GmapParamORNULL",
                        version = "POSIXltORNULL"),
         contains="GmapAlignerOutput",
         validity = .valid_GmapOutput)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

setMethod("alignmentCategories", "GmapOutput", function(x) {
              c("uniq", "mult", "chimera", "nomapping")
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

GmapOutput <- function(path, param = NULL, version = NULL) {
    newGmapAlignerOutput("GmapOutput", path, param, version)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

formatToExt <- function(x) {
    switch(x,
           gff3_gene=, gff3_match_cdna=, gff3_match_est="gff3",
           sampe=, samse="sam",
           psl="psl", splicesites="splicesites", introns="introns",
           map_exons="iit", map_ranges="iit", coords="coords")
}

setAs("GmapOutput", "RTLFileList", function(from) {
          as(lapply(paths(from), FileForFormat, formatToExt(from@param@format)),
             "List")
      })

setAs("GmapOutput", "RTLFile", function(from) {
          p <- paths(from)
          if (length(p) == 0L) {
              stop("no output files at: ", path(from))
          }
          FileForFormat(p[1L], formatToExt(from@param@format))
      })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

setMethod("import", c("GmapOutput", "missing", "missing"),
          function(con, format, text, ...) {
              if (con@param@format %in% c("samse", "sampe")) {
                  f <- as(con, "BamFile")
              } else {
                  f <- as(con, "RTLFile")
              }
              ans <- import(f, ...)
              seqinfo(ans) <- seqinfo(con@param@genome)
              ans
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### List class
###

setClass("GmapOutputList",
         prototype = prototype(elementType = "GmapOutput"),
         contains = "GmapAlignerOutputList")

setClass("SimpleGmapOutputList",
         contains = c("GmapOutputList", "SimpleList"))

GmapOutputList <- function(...) {
    args <- list(...)
    if (length(args) == 1 && is.list(args[[1]])) 
        args <- args[[1]]
    S4Vectors:::new_SimpleList_from_list("SimpleGmapOutputList", args)
}
