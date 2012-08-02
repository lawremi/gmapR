### =========================================================================
### iit_store command
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### High-level wrapper
###

setGeneric("iit_store", function(x, dest, ...) standardGeneric("iit_store"))

gmapRange <- function(x) {
  pos <- strand(x) == "+"
  range <- paste0(seqnames(x), ":", ifelse(pos, start(x), end(x)))
  if (!all(width(x) == 1L)) {
    range <- paste0(range, "..", ifelse(pos, end(x), start(x)))
  }
  range
}

setMethod("iit_store", c("GenomicRanges", "character"),
          function(x, dest, info = colnames(values(x))[1]) {
            lines <- paste0(">", names(x), " ", gmapRange(x), " ",
                            values(x)[[info]])
            p <- .iit_store(sort = "none", output = dest)
            writeLines(lines, p)
            close(p)
            dest
          })

setMethod("iit_store", c("character", "character"),
          function(x, dest, gff = file_ext(x) == "gff",
                   label = if (gff) "ID" else NULL)
          {
            .iit_store(gff = gff, label = label, sort = "none",
                       output = dest, inputfile = x)
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level wrapper
###

.iit_store <- function(fields = FALSE, gff = FALSE,
                       label = if (gff) "ID" else NULL,
                       sort = c("chrom", "none", "alpha", "numeric-alpha"),
                       output, .inputfile = NULL)
{
### TODO: assertions
  sort <- match.arg(sort)
  cl <- commandLine("iit_store")
  if (is.null(.inputfile))
    pipe(cl, open = "w")
  else .system(cl)
}
