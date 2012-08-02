### =========================================================================
### Routines for retrieving GMAP installation information
### -------------------------------------------------------------------------
###

gsnapInfo <- function() {
  info <- .gsnap(version = TRUE)
  getField <- function(name) {
    sub(".*?: ", "", info[grep(paste("^", name, ":", sep = ""), info)])
  }
  list(defaultDirectory = getField("Default gmap directory"))
}
