test_GmapGenomeDirectory <- function() {

  genomeDir <- tempdir()
  if (!file.exists(genomeDir)) dir.create(genomeDir)
  on.exit(unlink(genomeDir, recursive=TRUE))
  
  ##test constructor
  ggd <- GmapGenomeDirectory(path=genomeDir,
  create=TRUE)
  checkTrue(is(ggd, "GmapGenomeDirectory"))
  
  ##test methods
  checkIdentical(genome(ggd), character(0))
  checkIdentical(path(ggd), genomeDir)
}
