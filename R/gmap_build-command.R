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

setMethod("gmap_build", c("BSgenome", "GmapGenome"), function(x, genome, ...) {
  circular <- seqnames(x)[isCircular(x)]
  callNextMethod(x, genome, ..., circular=circular)
})

setMethod("gmap_build", c("FastaFile", "GmapGenome"), function(x, genome, ...) {
  gmap_build(path(x), genome, ...)
})

setMethod("gmap_build", c("RTLFile", "GmapGenome"), function(x, genome, ...) {
  gmap_build(import(x), genome, ...)
})

setMethod("gmap_build", c("FaFile", "GmapGenome"),
          function(x, genome, param, ...) {
              gmap_build(getSeq(x, param), genome, ...)
          })

setMethod("gmap_build", c("DNAStringSet", "GmapGenome"),
          function(x, genome, ...) {
            tmpfile <- tempfile()
            export(object=x, con=tmpfile, format="fasta")
            x <- FastaFile(tmpfile)
            gmap_build(x, genome, ...)
          })

setMethod("gmap_build", c("character", "GmapGenome"),
          function(x, genome, kmer = 15L, ...) {
            gmap_db_tmp_dir <- file.path(tempdir(), "gmap_db_tmp_dir")
            dir.create(gmap_db_tmp_dir, recursive=TRUE)
            cur_wd <- getwd()
            on.exit({unlink(gmap_db_tmp_dir, recursive=TRUE)
                setwd(cur_wd)})
            setwd(gmap_db_tmp_dir)

            gz <- file_ext(x) == "gz"
            x <- normalizePath(x, mustWork=TRUE)
            .gmap_build(db = genome(genome), dir = path(directory(genome)),
                        kmer = kmer, sort = "none", gunzip = gz, .fasta = x, ...)
            genome
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level wrapper
###

.gmap_build <- function(db, dir = NULL, kmer = 14L,
                        sort = c("chrom", "none", "alpha", "numeric-alpha"),
                        gunzip = FALSE, circular = NULL, .fasta, B = NULL) {
  sort <- match.arg(sort)
  if (!isSingleString(db)) {
    stop("'db' parameter (genome name) must be a single, non-NA string")
  }
  if (!is.null(dir) && !isSingleString(dir)) {
    stop("'dir' parameter (destination directory), if not NULL, must be a ",
         "'single, non-NA string")
  }
  kmer <- as.integer(kmer)
  if(!isSingleInteger(kmer) || kmer > 15 || kmer < 12) {
    stop("'kmer' parameter (k-mer value for genomic index) must be an integer ",
         "from 12 to 15")
  }
  if (!isTRUEorFALSE(gunzip)) {
    stop("'gunzip' parameter (files are gzipped, so need to gunzip each file ",
         "first) must be TRUE or FALSE")
  }
  if (!is.null(circular)) {
    circular <- paste(circular, collapse=",")
  }
  
  ##TODO: provide default
  if (is.null(dir)) {
    stop("destination directory must be explicitly specified")
  }

  B <- system.file("usr/bin/", package="gmapR", mustWork=TRUE)

  cl <- commandLine("gmap_build")
  res <- .system(cl)
  if (res != 0L) {
    stop("system call returned a non-0 status: ", cl)
  }
  invisible(TRUE)
}
