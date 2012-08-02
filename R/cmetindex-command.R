### =========================================================================
### cmetindex command
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### High-level wrapper
###

cmetindex <- function(db, use_snps = NULL) {
  if (!is(db, "GmapGenome"))
    stop("'db' must be a GmapGenome object")
  .cmetindex(genome(db), directory(db), directory(db), use_snps = use_snps)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level wrapper
###

.cmetindex <- function(db = "hg19", sourcedir = NULL, destdir = NULL,
                       use_snps = NULL)
{
### TODO: assertions
  .system(commandLine("cmetindex"))
}
