## gathers command line from non-NULL args to parent frame
commandLine <- function(binary = "gsnap") {
  args <- mget(names(formals(sys.function(-1L))), parent.frame())
  args <- Filter(Negate(is.null), args)
  paste(binary, paste(paste("--", gsub("_", "-", names(args)), sep = ""),
                      sapply(args, as.character), collapse = " "))
}

## low-level wrapper
.atoiindex <- function(db = "hg19", sourcedir = NULL, destdir = NULL,
                       use_snps = NULL)
{
  ## various assertions here
  system(commandLine("atoiindex"))
}

## high-level wrapper; could be a method on the Db object
atoiindex <- function(db, use_snps = NULL) {
  .atoiindex(name(db), directory(db), directory(db), use_snps = use_snps)
}
