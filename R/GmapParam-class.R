### =========================================================================
### GmapParam class
### -------------------------------------------------------------------------
###
### High-level interface to gmap.
###

setClassUnion("GmapSnps_OR_NULL", c("GmapSnps", "NULL"))

setClass("GmapAlignerParam",
         representation(genome = "GmapGenome",
                        part = "character_OR_NULL",
                        batch = "character",
                        snps = "GmapSnps_OR_NULL",
                        mode = "character",
                        nthreads = "integer",
                        npaths = "integer",
                        quiet_if_excessive = "logical",
                        nofails = "logical", 
                        split_output = "logical",
                        extra = "list"))

setClass("GmapParam",
         representation(suboptimal_score = "integer_OR_NULL",
                        splicing = "logical",
                        format = "character"),
         contains="GmapAlignerParam")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

GmapParam <- function(genome, unique_only = FALSE,
                      suboptimal_score = NULL, mode = "standard",
                      snps = NULL,
                      npaths = if (unique_only) 1L else 100L,
                      quiet_if_excessive = unique_only, nofails = unique_only,
                      split_output = !unique_only,
                      splicing = TRUE, 
                      nthreads = 1L, part = NULL, batch = "2",
                      format = c("gff3_gene", "gff3_match_cdna",
                          "gff3_match_est", "sampe", "samse", "psl",
                          "splicesites", "introns",
                          "map_exons", "map_ranges", "coords"),
                      ...)
{
    format <- match.arg(format)
    newGmapAlignerParam("GmapParam", genome, snps)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

setAs("GmapParam", "list", function(from) {
          to <- GmapAlignerParam_asList(from)
          to$nosplicing <- !to$splicing
          to$splicing <- NULL
          to
      })
