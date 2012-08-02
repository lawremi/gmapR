### =========================================================================
### snpindex command
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### High-level wrapper
###

snpindex <- function(name, genome, destdir = NULL,
                     iitfile = paste(name, "iit", sep = "."))
{
  .snpindex(sourcedir = path(directory(genome)), db = genome(genome),
            destdir = destdir, snpsdb = name, .iitfile = iitfile)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level wrapper
###

.snpindex <- function(sourcedir = NULL, db, destdir = NULL, snpsdb,
                      .iitfile = NULL)
{
  .system(commandLine("snpindex"))
}
