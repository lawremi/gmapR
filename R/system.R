### =========================================================================
### System call utilities
### -------------------------------------------------------------------------
###

getDefaultGmapPath <- function() {
  system.file("usr", "bin", package = "gmapR")
}

.commandLine.handleNameArgs <- function(defaultArgs, args) {
  ##arguments with no default need to be handled explicity, since
  ##objects of class "missing" cannot be eval'ed.
  classes <- sapply(defaultArgs, class)
  handledArgs <- defaultArgs
  handledArgs[classes == "name"] <- "NoT_pRoViDeD"

  ##verify function calling commandLine provided values for any
  ##argument with no default
  callerMustSupply <- names(classes[classes == "name"])
  if (any(class(args[callerMustSupply]) == "name")) {
    stop("An argument with a required value was not passed to commandLine()")
  }

  return(handledArgs)
}

## gathers command line from non-NULL args to parent frame
commandLine <- function(binary = "gsnap",
                        path = getOption("gmap.path", getDefaultGmapPath()))
{
  ##get values of arguments of function that called this function
  args <- mget(names(formals(sys.function(sys.parent()))), parent.frame())

  ##Some args cannot be passed simultaneously to gsnap. Remove arg if
  ##it is the default.
  defaultArgs <- formals(sys.function(sys.parent()))
  ##objects of class "name" must be handled specially since they
  ##cannot be eval'ed
  defaultArgs <- .commandLine.handleNameArgs(defaultArgs, args)
  parentEnv <- sys.parent()
  defaultArgs <- lapply(defaultArgs, eval, envir=parentEnv)  
  ##some args have a default vector with more than one element (for
  ##use with match.args). Taking first element:
  defaultArgs <- lapply(defaultArgs,
                        function(x) {if (length(x) > 1L) x[[1L]] else x}
                        )
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
