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

setClass("IIT", representation(ptr = "externalptr"))
setClassUnion("IITORNULL", c("IIT", "NULL"))

setRefClass("IITFile",
            fields=c(iit="IITORNULL"),
            contains="RTLFile")

IITFile <- function(resource) {
  new("IITFile", resource = resource)
}

IIT <- function(ptr) {
    new("IIT", ptr=ptr)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

iit <- function(x) x$iit
`iit<-` <- function(x, value) {
    x$iit <- value
    x
}

fieldNames <- function(x) {
    .Call(R_iit_fieldNames(x), iit(x))
}

typeNames <- function(x) {
    .Call(R_iit_typeNames(x), iit(x))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Open/close
###

open.IITFile <- function(con, ranges = TRUE, labels = TRUE) {
    stopifnot(isTRUEorFALSE(ranges),
              isTRUEorFALSE(labels),
              !isOpen(con))
    iit(con) <- IIT(.Call(R_iit_read, resource(con), ranges, labels))
    invisible(con)
}

close.IITFile <- function(con) {
    iit(con) <- NULL
    invisible(con)
}

setMethod("isOpen", "IITFile", function(con, rw=c("", "read", "write")) {
    rw <- match.arg(rw)
    if (rw == "write")
        FALSE
    else !is.null(iit(con))
})

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

normArgWhich <- function(which) {
    if (!is.null(which) && !is.character(which)) {
        which <- as(which, "GRanges")
        sign <- as.integer(strand(which))
        sign[sign == 1L] <- 4L
        sign <- sign - 3L
        which <- list(as.character(seqnames(which)), start(which), end(which),
                      sign)
    }
    which
}

openForWhich <- function(con, which) {
    labels <- is.character(which)
    if (labels && anyNA(which))
        stop("query labels cannot be NA")
    ranges <- !labels && !is.null(which)
    if (ranges)
        which <- as(which, "GRanges")
    open(con, ranges=ranges, labels=labels)
    which
}

normArgFields <- function(con, fields) {
    iit_fields <- fieldNames(con)
    bad_fields <- setdiff(fields, iit_fields)
    if (length(bad_fields) > 0L)
        stop("invalid 'fields' requested: ", paste(bad_fields, collapse=", "))
    fields
}

normArgType <- function(con, type) {
    if (is.null(type))
        return(NULL)
    if (!isSingleString(type))
        stop("'type' must be NULL or a single, non-NA string")
    iit_types <- typeNames(con)
    if (!(type %in% iit_types))
        stop("'type' must be one of ", paste(iit_types, collapse=", "))
    type
}

### NOTE: if 'exact=TRUE' and is.null(type), only gets intervals with NO type
setMethod("import", "IITFile",
          function(con, format, text, which = NULL, type = NULL,
                   fields = fieldNames(con),
                   ignore.strand = FALSE, exact = FALSE,
                   as = c("GRanges", "GRangesList", "DataFrame"))
{
    if (!missing(format))
        checkArgFormat(con, format)
    if (!missing(text))
        stop("'text' not supported")
    if (!isOpen(con)) {
        openForWhich(con, which)    
    }
    type <- normArgType(con, type)
    if (!is.null(type) && is.character(which))
        stop("cannot restrict by 'type' when 'which' is character")
    fields <- normArgFields(con, fields)
    stopifnot(isTRUEorFALSE(ignore.strand),
              isTRUEorFALSE(exact))
    as <- match.arg(as)
    exact.strand <- exact && !ignore.strand && !is.null(which)
    include.ranges <- as %in% c("GRanges", "GRangesList") || exact.strand
    
    result <- .Call(R_iit_read, iit(con), which, type, fields, ignore.strand,
                    exact, include.ranges)
    names(ans) <- c(if (include.ranges)
                        c("seqnames", "start", "width", "strand"),
                    "which", "mcols")
    mcols <- DataFrame(result["which"], setNames(result[["mcols"]], fields))
    result[["mcols"]] <- NULL
    result[["which"]] <- NULL
    
    if (as %in% c("GRanges", "GRangesList")) {
        ans <- makeGRangesFromDataFrame(DataFrame(result, mcols))
        if (as == "GRangesList")
            ans <- split(ans, mcols[["which"]])
    } else {
        ans <- mcols
    }
    
    if (exact.strand) {
        which <- which[mcols[["which"]]]
        ans <- extractROWS(ans, strand(which) == "*" |
                                result[["strand"]] == "*" |
                                ans[["strand"]] == strand(which))
    }

    ans
})
