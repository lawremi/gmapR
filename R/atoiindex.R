## gathers command line from non-NULL args to parent frame
commandLine <- function(binary = "gsnap") {
  args <- mget(names(formals(sys.function(sys.parent()))), parent.frame())
  args <- Filter(Negate(is.null), args)
  args <- Filter(function(x) !identical(x, FALSE), args)
  named <- !grepl("^\\.", names(args))
  args[!named] <- paste(sapply(args[!named], as.character), collapse = " ")
  args[named] <- paste(paste(ifelse(nchar(names(args[named])) > 1, "--", "-"),
                             gsub("_", "-", names(args[named])), sep = ""),
                       sapply(args[named], as.character))
  binary <- file.path(globals()['gmap_bin'], binary)
  paste(binary, paste(args, collapse = " "))
}

## at some point, may want to customize this
.system <- function(...) system(...)

## low-level wrapper
.atoiindex <- function(db = "hg19", sourcedir = NULL, destdir = NULL,
                       use_snps = NULL)
{
  ## various assertions here
  .system(commandLine("atoiindex"))
}

## high-level wrapper; could be a method on the Db object
atoiindex <- function(db, use_snps = NULL) {
  .atoiindex(name(db), directory(db), directory(db), use_snps = use_snps)
}
