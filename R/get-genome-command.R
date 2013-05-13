### =========================================================================
### get-genome program
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### High-level wrapper
###

setGeneric("get_genome", function(x, ...) standardGeneric("get_genome"))

setMethod("get_genome", "GmapGenome",
          function(x, which = NULL, snps = NULL, ...)
          {
            if (!is.null(which))
              range <- gmapRange(which)
            else range <- ""
            usesnps <- snpsdir <- NULL
            if (!is.null(snps)) {
              usesnps <- name(snps)
              snpsdir <- directory(snps)
            }
            tmpfile <- file.path(tempdir(), "get_genome.fasta")
            .get_genome(path(directory(db)), genome(db), snpsdir = snpsdir,
                        usesnps = usesnps, snpformat = "2", ...,
                        .range = range, .redirect = tmpfile)
            readDNAStringSet(tmpfile)
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level wrapper
###

.get_genome <- function(dir, db, dibase = FALSE, strain = NULL, coords = FALSE,
                        uppercase = FALSE, wraplength = 60L, fullgenome = FALSE,
                        header = NULL, snpsdir = NULL, usesnps = NULL,
                        snpformat = c("3", "0", "1", "2"), mapdir = NULL,
                        map = NULL, relative = FALSE, ranks = FALSE,
                        raw = FALSE, flanking = 0L, dump = FALSE,
                        chromosomes = FALSE, contigs = FALSE, .range = "",
                        .redirect = NULL)
{
### TODO: assertions
  snpformat <- match.arg(snpformat)
  command <- commandLine("get-genome")
  if (!is.null(.redirect)) {
    paste(command, ">", .redirect)
    .system(command)
  } else pipe(command)
}
