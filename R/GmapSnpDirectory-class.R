### =========================================================================
### SnpDirectory class
### -------------------------------------------------------------------------
###
### Database of SNPs used by the GMAP suite.
###

setClass("GmapSnpDirectory", representation(path = "character"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

setMethod("path", "GmapSnpDirectory", function(object) object@path)

setMethod("names", "GmapSnpDirectory", function(x) {
  snp_files <- dir(path(x), ".*genomecomp\\.")
  sub("\\.genomecomp\\.", ":", snp_files)
})

setMethod("length", "GmapSnpDirectory", function(x) {
  length(names(x))
})

setReplaceMethod("snps", c("GmapSnpDirectory", "VCF"),
                 function(x, name, genome = GmapGenome(genome(x)),
                          iitPath = tempdir(), ..., value)
                 {
                   gr <- rowData(value)
                   values(gr) <- values(fixed(value))[c("REF", "ALT")]
                   iitFile <- file.path(iitPath, paste(name, "iit", sep = "."))
                   alt <- values(gr)$ALT
                   if (is(alt, "List")) {
                     gr <- rep(gr, elementLengths(alt))
                     alt <- unlist(alt)
                   }
                   ref <- values(gr)$REF
                   single <- nchar(alt) == 1L & nchar(ref) == 1L
                   change <- paste(ref[single], alt[single], sep = "")
                   gr <- gr[single]
                   values(gr) <- DataFrame(change)
                   export.iit(gr, iitFile)
                   snpindex(name, genome, path(x), iitFile)
                   x
                 })

setReplaceMethod("snps", c("GmapSnpDirectory", "character"),
                 function(x, name, genome, which, ..., value)
                 {
                   if (missing(genome)) {
                     stop("Please supply the \"genome\" argument")
                   }

                   param <- ScanVcfParam(fixed = "ALT", info = NA, geno = NA)
                   if (!missing(which)) # FIXME: waiting for vcfWhich<-
                     param@which <- as(which, "RangesList")
                   snps(x, name = name, genome = genome, ...) <-
                     readVcf(value, genome(genome), param)
                   x
                 })

setMethod("[[<-", c("GmapSnpDirectory", value="ANY"),
          function(x, i, j, ..., value) {
            if (!missing(j))
              warning("argument 'j' ignored")
            snps(x, name = i, ...) <- value
            x
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

GmapSnpDirectory <- function(path, create = FALSE) {
  if (is(path, "GmapGenome"))
    path <- path(path)
  if (!isSingleString(path))
    stop("'path' must be a single, non-NA string")
  if (!isTRUEorFALSE(create))
    stop("'create' must be TRUE or FALSE")
  if (create) {
    if (file.exists(path))
      message("NOTE: snp directory '", path, "' already exists, not recreating")
    else dir.create(path, recursive = TRUE)
  }
  new("GmapSnpDirectory", path = path)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

setMethod("show", "GmapSnpDirectory", function(object) {
  cat("GmapSnpDirectory object\n", "path: ", path(object),
      "\nnames: ", names(object), "\n", sep = "")
})
