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
  if (is.null(destdir) && !file.exists(destdir)) {
    stop("The destination directory for the SNPs has not been created.",
         " Perhaps run GmapSnpDirectory with the \"create\" arg set to TRUE.")
  }
  .system(commandLine("snpindex"))
}
