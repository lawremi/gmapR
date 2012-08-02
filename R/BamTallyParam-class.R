### =========================================================================
### BamTallyParam class
### -------------------------------------------------------------------------
###
### Collects parameters for the bam_tally command
###

setClass("BamTallyParam",
         representation(which = "RangesList",
                        cycle_breaks = "integerORNULL",
                        high_quality_cutoff = "integer",
                        minimum_mapq = "integer",
                        concordant_only = "logical",
                        unique_only = "logical",
                        primary_only = "logical",
                        min_depth = "integer",
                        variant_strand = "integer",
                        ignore_query_Ns = "logical",
                        indels = "logical"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

BamTallyParam <- function(which = RangesList(), cycle_breaks = NULL,
                          high_quality_cutoff = 0L,
                          minimum_mapq = 0L,
                          concordant_only = FALSE, unique_only = FALSE,
                          primary_only = FALSE,
                          min_depth = 0L, variant_strand = 0L,
                          ignore_query_Ns = FALSE,
                          indels = FALSE)
{
  args <- names(formals(sys.function()))
  params <- mget(args, environment())
  params$which <- as(which, "RangesList")
  do.call(new, c("BamTallyParam", params))  
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

setAs("BamTallyParam", "list", function(from) {
  sapply(slotNames(from), slot, object = from, simplify = FALSE)
})

setMethod("as.list", "BamTallyParam", function(x) as(x, "list"))
