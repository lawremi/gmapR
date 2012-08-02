### =========================================================================
### GmapGenomeDirectory class
### -------------------------------------------------------------------------
###
### A directory containing one or more GMAP databases (GmapGenome objects)
###

setClass("GmapGenomeDirectory", representation(path = "character"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

setMethod("path", "GmapGenomeDirectory", function(object) object@path)
setMethod("path", "NULL", function(object) NULL)

setMethod("genome", "GmapGenomeDirectory", function(x) {
  paths <- dir(path(x), full.names = TRUE)
  is_dir <- file.info(paths)[,"isdir"]
  genomes <- basename(paths[is_dir])
  genomes[file.exists(file.path(paths[is_dir], paste0(genomes, ".version")))]
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

getDefaultGmapGenomePath <- function() {
  path.expand(file.path(Sys.getenv("XDG_DATA_HOME", "~/.local/share"), "gmap"))
}

GmapGenomeDirectory <- function(path = getDefaultGmapGenomePath(),
                                create = FALSE)
{
  if (!isSingleString(path))
    stop("'path' must be a single, non-NA string")
  if (create && !file.exists(path)) {
    message("Creating directory ", path)
    dir.create(path, recursive=TRUE)
  }
  new("GmapGenomeDirectory", path = file_path_as_absolute(path))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

setMethod("show", "GmapGenomeDirectory", function(object) {
  cat(class(object), "object\npath:", path(object), "\n")
})
