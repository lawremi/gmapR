### =========================================================================
### System call utilities
### -------------------------------------------------------------------------
###

getDefaultGmapPath <- function() {
  system.file("usr", "bin", package = "gmapR")
}

## gathers command line from non-NULL args to parent frame
commandLine <- function(binary = "gsnap",
                        path = getOption("gmap.path", getDefaultGmapPath()))
{

  #Some args cannot be passed simultaneously to gsnap. Remove arg if
  #it is the default.
  defaultArgs <- formals(sys.function(sys.parent()))
  ##some args have a default vector with more than one element (for
  ##use with match.args). Taking first element:
  multis <- which(elementLengths(defaultArgs) > 1)
  multiVals <- lapply(defaultArgs[multis],
                      function(x) {
                        retval <- switch(class(x),
                                         list=unlist(x)[1],
                                         call=eval(x)[1],
                                    stop("Cannot get default value for an object of class",
                                         class(x)))                                                
                        return(retval)
                      })
  defaultArgs[multis] <- multiVals
  
  args <- mget(names(formals(sys.function(sys.parent()))), parent.frame())

  ##remove defaults
  isDefault <- mapply(identical, args, defaultArgs)
  if (sum(!isDefault > 0L)) {
    args <- args[!isDefault]
  }
  
  args <- Filter(Negate(is.null), args)
  args <- Filter(function(x) !identical(x, FALSE), args)
  named <- !grepl("^\\.", names(args))
  unnamed_args <- sapply(args[!named], as.character)
  named_args <- paste(ifelse(nchar(names(args[named])) > 1, "--", "-"),
                            gsub("_", "-", names(args[named])), sep = "")
  toggle_arg <- sapply(args[named], isTRUE)
  named_args[!toggle_arg] <-
    paste(named_args[!toggle_arg],
          sapply(args[named][!toggle_arg], as.character))
  if (!is.null(path))
    binary <- file.path(path, binary)
  paste(binary, paste(c(named_args, unnamed_args), collapse = " "))
}

## at some point, mxbay want to customize this
.system <- function(...) system(...)
