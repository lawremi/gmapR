test_GmapGenome_constructor_DNAStringSet_noCreate <- function() {
  dna <- Biostrings::DNAStringSet("ACTGTGTCAG")
  gmapGenome <- GmapGenome(genome=dna, name="thing")
  checkTrue(is(gmapGenome, "GmapGenome"))
}

test_GmapGenome_constructor_DNAStringSet_create <- function() {
  dna <- Biostrings::DNAStringSet("ACTGTGTCAG")
  genomeDir <- file.path(tempdir(), as.integer(runif(1) * 1000000000))
  if (file.exists(genomeDir)) unlink(genomeDir, recursive=TRUE)
  dir.create(genomeDir, recursive=TRUE)
  
  on.exit(unlink(genomeDir, recursive=TRUE))
  gmapGenome <- GmapGenome(genome=dna, directory=genomeDir,
                           name="thing", create=TRUE)                           
  checkTrue(is(gmapGenome, "GmapGenome"))
}

test_GmapGenome_constructor_BSgenome_create <- function() {
  library("BSgenome.Scerevisiae.UCSC.sacCer3")
  genomeName <- "yeast"  
  genomeDir <- file.path(tempdir(), as.integer(runif(1) * 1000000000))
  if (file.exists(genomeDir)) unlink(genomeDir, recursive=TRUE)
  dir.create(genomeDir, recursive=TRUE)  
  on.exit(unlink(genomeDir, recursive=TRUE))
  gmapGenome <- GmapGenome(genome=Scerevisiae, directory=genomeDir,
                           name=genomeName, create=TRUE)                           
  checkTrue(is(gmapGenome, "GmapGenome"))
}

test_GmapGenome_accessors <- function() {
  genomeName <- "testGenome"
  dna <- Biostrings::DNAStringSet("ACTGTGTCAG")
  genomeDir <- file.path(tempdir(), as.integer(runif(1) * 1000000000))
  if (file.exists(genomeDir)) unlink(genomeDir, recursive=TRUE)
  dir.create(genomeDir, recursive=TRUE)
  
  on.exit(unlink(genomeDir, recursive=TRUE))
  gmapGenome <- GmapGenome(genome=dna, directory=genomeDir,
                           name=genomeName, create=TRUE)

  checkIdentical(path(gmapGenome), file.path(genomeDir, genomeName))
  checkTrue(is(directory(gmapGenome), "GmapGenomeDirectory"))
  checkIdentical(genome(gmapGenome), genomeName)
}