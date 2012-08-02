### =========================================================================
### atoiindex program
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### High-level wrapper
###

atoiindex <- function(db, use_snps = NULL) {
  if (!is(db, "GmapGenome"))
    stop("'db' must be a GmapGenome object")
  .atoiindex(genome(db), directory(db), directory(db), use_snps = use_snps)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level wrapper
###

.atoiindex <- function(db = "hg19", sourcedir = NULL, destdir = NULL,
                       use_snps = NULL)
{
### TODO: assertions
  .system(commandLine("atoiindex"))
}
