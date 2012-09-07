### =========================================================================
### GmapGenome class
### -------------------------------------------------------------------------
###
### Database of reference sequence used by the GMAP suite.
###

setClass("GmapGenome", representation(name = "character",
                                      directory = "GmapGenomeDirectory"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

directory <- function(x) {
  retVal <- NULL
  if (!is.null(x)) {
    retVal <- x@directory ## 'dir' is already taken
  }
  return(retVal)
}

setMethod("genome", "GmapGenome", function(x) x@name)

setMethod("path", "GmapGenome",
          function(object) file.path(path(directory(object)), genome(object)))

mapsDirectory <- function(x) {
  file.path(path(x), paste(genome(x), "maps", sep = "."))
}

setMethod("seqinfo", "GmapGenome", function(x) {
  tab <- read.table(.get_genome(path(directory(x)), genome(x),
                                chromosomes = TRUE),
                    colClasses = c("character", "NULL", "integer"))
  Seqinfo(tab[,1], tab[,2], genome = genome(x))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

setGeneric("genomeName", function(x) standardGeneric("genomeName"))

setMethod("genomeName", "character", function(x) x)
setMethod("genomeName", "BSgenome", function(x) providerVersion(x))
setMethod("genomeName", "FastaFile",
          function(x) file_path_sans_ext(basename(path(x)), TRUE))

GmapGenome <- function(genome,
                       directory = GmapGenomeDirectory(create = create), 
                       name = genomeName(genome), create = FALSE, ...)
{
  if (!isTRUEorFALSE(create))
    stop("'create' must be TRUE or FALSE")
  if (isSingleString(directory))
    directory <- GmapGenomeDirectory(directory, create = create)
  if (is(genome, "DNAStringSet")) {
    if (missing(name))
      stop("If the genome argument is a DNAStringSet object",
           "the name argument must be provided") 
  }
  if (!isSingleString(name))
    stop("'name' must be a single, non-NA string")
  if (!is(directory, "GmapGenomeDirectory"))
    stop("'directory' must be a GmapGenomeDirectory object or path to one")
  db <- new("GmapGenome", name = name, directory = directory)
  if (create) {
    if (name %in% genome(directory))
      message("NOTE: genome '", name, "' already exists, not overwriting")
    else referenceSequence(db, ...) <- genome
  }
  db
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Building
###

setReplaceMethod("referenceSequence",
                 signature(x = "GmapGenome", value = "ANY"),
                 function(x, name, ..., value)
                 {
                   gmap_db_tmp_dir <- file.path(tempdir(), "gmap_db_tmp_dir")
                   dir.create(gmap_db_tmp_dir, recursive=TRUE)
                   cur_wd <- getwd()
                   on.exit({unlink(gmap_db_tmp_dir, recursive=TRUE)
                            setwd(cur_wd)})
                   setwd(gmap_db_tmp_dir)

                   db <- gmap_build(value, x, ...)
                   
                   db
                 })

setGeneric("snps<-", function(x, name, ..., value) standardGeneric("snps<-"))

setReplaceMethod("snps", c("GmapGenome", "ANY"),
                 function(x, name, ..., value) {
                   snpDir <- GmapSnpDirectory(x)
                   snps(snpDir, name = name, genome = x,
                        iitPath = mapsDirectory(x), ...) <- value
                   x
                 })

setGeneric("spliceSites<-",
           function(x, ..., value) standardGeneric("spliceSites<-"))

setReplaceMethod("spliceSites", c("GmapGenome", "GRangesList"),
                 function(x, name, value) {
                   exonsFlat <- unlist(value, use.names=FALSE)
                   exonsPart <- PartitioningByWidth(exons)
                   exonsHead <- exonsFlat[-end(exonsPart)]
                   donors <- flank(exonsHead, 1L, start = FALSE)
                   exonsTail <- exonsFlat[-start(exonsPart)]
                   acceptors <- flank(exonsTail, 1L, start = TRUE)
                   sites <- c(resize(donors, 2L, fix = "end"),
                              resize(acceptors, 2L, fix = "start"))
                   names(sites) <- values(sites)$exon_id
                   info <- rep(c("donor", "acceptor"), each = length(donors))
                   intronWidths <- abs(acceptors - donors) + 1L
                   info <- paste(info, intronWidths)
                   values(sites) <- DataFrame(info)
                   iit_store(sites, file.path(mapsDirectory(x), name))
                   x
                 })

setReplaceMethod("spliceSites", c("GmapGenome", "TranscriptDb"),
                 function(x, name, value) {
                   spliceSites(x, name) <- exonsBy(value)
                   x
                 })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coerce
###

setAs("ANY", "GmapGenome", function(from) GmapGenome(from))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

setMethod("show", "GmapGenome", function(object) {
  cat("GmapGenome object\ngenome:", genome(object), "\ndirectory:",
      path(directory(object)), "\n")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

.normArgDb <- function(db, dir) {
  if (!is(db, "GmapGenome"))
    db <- GmapGenome(db, dir)
  db
}

