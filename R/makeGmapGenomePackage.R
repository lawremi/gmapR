##This code borrows very heavily from GenomicFeatures::makeTxDbPackage
makeGmapGenomePackage <- function(gmapGenome,
                                  version,
                                  maintainer,
                                  author,
                                  destDir=".",
                                  license="Artistic-2.0",
                                  pkgName){

  if (missing(pkgName)) {
    stop("Please supply the package name")
  }

  template_path <- system.file("GmapGenomePackage-template",
                               package="gmapR",
                               mustWork=TRUE)
  ## We need to define some symbols in order to have the 
  ## template filled out correctly.
  symvals <- list(PKGTITLE=paste("GmapGenome package for the genome", genome(gmapGenome)),
                  PKGNAME=pkgName,
                  PKGDESCRIPTION=paste("Contains the GmapGenome object and IIT files for the genome",
                    genome(gmapGenome)),
                  PKGVERSION=version,
                  AUTHOR=author,
                  MAINTAINER=maintainer,
                  GMAPRVERSION=as.character(packageVersion("gmapR")),
                  LIC=license,
                  GMAPOBJNAME=gsub("\\.", "_", pkgName))
  
  ## Should never happen
  if (any(duplicated(names(symvals)))) {
    str(symvals)
    stop("'symvals' contains duplicated symbols")
  }

  ## All symvals should by single strings (non-NA)
  is_OK <- sapply(symvals, isSingleString)
  if (!all(is_OK)) {
    bad_syms <- paste(names(is_OK)[!is_OK], collapse=", ")
    stop("values for symbols ", bad_syms, " are not single strings")
  }
  
  if (!file.exists(destDir)) dir.create(destDir, recursive=TRUE)
  createPackage(pkgname=pkgName,
                destinationDir=destDir,
                originDir=template_path,
                symbolValues=symvals)
  
  ## then copy the contents of the GmapGenome into the extdata dir
  iit_from_path <- path(gmapGenome)
  iit_dest_path <- file.path(destDir, pkgName, "inst", "extdata")
  if (!file.exists(iit_dest_path)) dir.create(iit_dest_path, recursive=TRUE)
  success <- file.copy(from=iit_from_path,
                       to=iit_dest_path, recursive=TRUE)
  if (!success) stop("Could not create package.",
                     "IIT files would not copy from source to destination")
  
  invisible(TRUE)
}
