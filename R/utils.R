### =========================================================================
### Utilities
### -------------------------------------------------------------------------

checkPackageVersion <- function(pkg, version, required = FALSE) {
  installed.version <- try(packageVersion(pkg), silent = TRUE)
  success <- if (is(installed.version, "try-error"))
    FALSE
  else installed.version >= version
  if (required && !success)
    stop("Package '", pkg, "' (version >= '", version, "') required")
  success
}

##' Checks if a given package is installed.
##'
##' @title Checks if a given package is installed.
##' @param pkg A character string containing a package name.
##' @param required A boolean. The functions stops if set to \code{TRUE} and if the required package is not present. Default is \code{FALSE}.
##' @return A boolean.
##' @author Cory Barr
##' @export
checkPackageInstalled <- function(pkg, required = FALSE) {
  checkPackageVersion(pkg, "0.0.0", required)
}
