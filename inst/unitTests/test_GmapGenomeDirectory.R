test_GmapGenomeDirectory <- function() {

  genomeDir <- tempdir()
  if (!file.exists(genomeDir)) dir.create(genomeDir)
  on.exit(unlink(genomeDir, recursive=TRUE))
  
  ##test constructor
  genomeDir <- tools::file_path_as_absolute(genomeDir)
  ggd <- GmapGenomeDirectory(path=genomeDir, create=TRUE)
  checkTrue(is(ggd, "GmapGenomeDirectory"))
  
  ##test methods
  checkIdentical(genome(ggd), character(0))

  ##for compatibility with Macs
  .unSymLink <- function(d) {
    pieces <- unlist(strsplit(d, "/"))
    unSymed <- Sys.readlink(file.path(pieces[1], pieces[2]))
    if (unSymed != "") {
      pieces[2] <- sub("^/", "", unSymed)
      d <- paste(pieces, collapse="/")
    }
    return(d)
  }
  checkIdentical(.unSymLink(path(ggd)), .unSymLink(genomeDir))
}
