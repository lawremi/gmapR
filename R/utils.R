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

checkPackageInstalled <- function(pkg, required = FALSE) {
  checkPackageVersion(pkg, "0.0.0", required)
}
