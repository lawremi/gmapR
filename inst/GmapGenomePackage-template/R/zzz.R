###
### Load the GmapGenome object whenever the package is loaded.
###

.onLoad <- function(libname, pkgname)
{
  ns <- asNamespace(pkgname)

  iit_path <- system.file("extdata", package=pkgname, lib.loc=libname)
  genome <- GmapGenome(genome = pkgname,
                       directory = GmapGenomeDirectory(iit_path),
                       name = pkgname,
                       create = FALSE)

  objname <- pkgname
  assign(objname, genome, envir=ns)
  namespaceExport(ns, objname)
}
