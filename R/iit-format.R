### =========================================================================
### IIT support
### -------------------------------------------------------------------------
###
### This is the primary format for intervals for GMAP. Simple wrappers
### around the GMAP commands.
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### IITFile class
###

setClass("IITFile", contains = "RTLFile")

IITFile <- function(resource) {
  new("IITFile", resource = resource)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Export
###

setGeneric("export.iit",
           function(object, con, ...) standardGeneric("export.iit"))

setMethod("export.iit", "ANY",
          function(object, con, ...)
          {
            export(object, con, "iit", ...)
          })

setMethod("export", c("ANY", "IITFile"), function(object, con, ...) {
  iit_store(object, path(con), ...)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Import
###
