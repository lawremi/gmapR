### =========================================================================
### gmap_build program
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### High-level wrapper
###

setGeneric("gmap_build", function(x, genome, ...) standardGeneric("gmap_build"))

setMethod("gmap_build", c("ANY", "GmapGenome"), function(x, genome, ...) {
  if (length(x) > 1 && is.null(names(x)))
    stop("if 'x' has multiple elements, it must have names")
  tmp_file <- tempfile("gmap_build_fasta")
  export(x, tmp_file, format = "fasta")
  gmap_build(tmp_file, genome = genome, ...)
})

setMethod("gmap_build", c("FastaFile", "GmapGenome"), function(x, genome, ...) {
  gmap_build(path(x), genome, ...)
})

setMethod("gmap_build", c("DNAStringSet", "GmapGenome"), function(x, genome, ...) {

  tmpfile <- tempfile()
  export(object=x, con=tmpfile, format="fasta")
  x <- FastaFile(tmpfile)
  
  gmap_build(x, genome, ...)
})


setMethod("gmap_build", c("character", "GmapGenome"),
          function(x, genome, kmer = 15L) {
            gz <- file_ext(x) == "gz"
            .gmap_build(d = genome(genome), D = path(directory(genome)), x,
                        k = kmer, s = "chrom", g = gz)
            genome
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level wrapper
###


##gmap_build args from command line
##Options:
##    -d STRING   genome name
##    -D STRING   destination directory for installation (defaults to gmapdb
##                directory specified at configure time)
##    -k INT      k-mer value for genomic index (allowed: 12..15, default 14)
##    -S          do not order chromosomes in numeric/alphabetic order, but
##                use order in FASTA file(s)
##    -g          files are gzipped, so need to gunzip each file first
##    -w INT      wait (sleep) this many seconds after each step

##NOTE: these options must be from a different version of
##gmap_build. -S is now -s and wants a string instead of being boolean
.gmap_build <- function(d, D = NULL, k = 14L,
                        s = c("none", "alpha", "numeric-alpha", "chrom"),
                        g = FALSE, .fasta, B = NULL) {
  s <- match.arg(s)
  if (!isSingleString(d)) {
    stop("'d' parameter (genome name) must be a single, non-NA string")
  }
  if (!is.null(D) && !isSingleString(D)) {
    stop("'D' parameter (destination directory), if not NULL, must be a ",
         "'single, non-NA string")
  }
  k <- as.integer(k)
  if(!isSingleInteger(k) || k > 15 || k < 12) {
    stop("'k' parameter (k-mer value for genomic index) must be an integer ",
         "from 12 to 15")
  }
  if (!isTRUEorFALSE(g)) {
    stop("'g' parameter (files are gzipped, so need to gunzip each file ",
         "first) must be TRUE or FALSE")
  }
  
  ##TODO: provide default
  if (is.null(D)) {
    stop("destination directory must be explicitly specified")
  }

  B <- system.file("usr/bin/", package="gmapR",
                            mustWork=TRUE)

  cl <- commandLine("gmap_build")
  res <- .system(cl)
  if (res != 0L) {
    stop("system call returned a non-0 status: ", cl)
  }
  invisible(TRUE)
}
