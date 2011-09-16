
gmapMake <- function(suffix, rule) {
  sys_command <- paste("make -f", paste('Makefile.', suffix, sep=''), rule)
  .system(sys_command)
}

## there could be a 'build' method on GmapDb
## this would then become a convenience wrapper
buildGmapDb <- function(x, db = deparse(substitute(x)), dir = NULL) {
  gmap_db_tmp_dir <- file.path(tempdir(), "gmap_db_tmp_dir")
  dir.create(gmap_db_tmp_dir, recursive=TRUE)
  cur_wd <- getwd()
  on.exit({unlink(gmap_db_tmp_dir, recursive=TRUE)
           setwd(cur_wd)})
  setwd(gmap_db_tmp_dir)

  db <- gmap_setup(x, db = db, dir = dir)
    
  gmapMake(name(db), "coords")
  gmapMake(name(db), "gmapdb")
  gmapMake(name(db), "install")

  db
}

setClass("GmapDb", representation(name = "character", directory = "character"))

directory <- function(x) x@directory ## 'dir' is already taken
name <- function(x) x@name

GmapDb <- function(name, directory) {
  if (!file.exists(directory)) {
    message("Creating directory ", directory)
    dir.create(directory, recursive=TRUE)
  }
  new("GmapDb", name = name, directory = file_path_as_absolute(directory))
}

.normArgDb <- function(db, dir) {
  if (!is(db, "GmapDb"))
    db <- GmapDb(db, dir)
  db
}

### FIXME: in latest gsnap this has already been renamed to gmap_build
setGeneric("gmap_setup", function(x, ...) standardGeneric("gmap_setup"))

setMethod("gmap_setup", "DNAStringSet",
          function(x, db = deparse(substitute(x)), dir = NULL)
          {
            if (length(x) > 1 && is.null(names(x)))
              stop("if 'x' has multiple elements, it must have names")
            tmp_file <- tempfile("gmap_setup_fasta")
            ## on.exit({
            ##   message("DELETING THE FILE")
            ##   unlink(tmp_file, recursive=TRUE)
            ## })
            write.XStringSet(x, tmp_file)
            gmap_setup(tmp_file, db = db, dir = dir)
          })

setMethod("gmap_setup", "character",
          function(x, db = deparse(substitute(x)), dir = NULL)
          {
            db <- .normArgDb(db, dir)
            .gmap_setup(d = name(db), D = directory(db), x)
            db
          })

.gmap_setup <- function(d, D = NULL, .fasta) {
  .system(commandLine("gmap_setup"))
}
