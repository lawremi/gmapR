test_GmapGenome_constructor_DNAStringSet_noCreate <- function() {
  dna <- Biostrings::DNAStringSet("ACTGTGTCAG")
  names(dna) <- "test"
  gmapGenome <- GmapGenome(genome=dna, name="thing", create=FALSE)
  checkTrue(is(gmapGenome, "GmapGenome"))
}

test_GmapGenome_constructor_DNAStringSet_create <- function() {
  set.seed(1)
  seq <- paste0(sample(c("A", "C", "G", "T"), 2000, replace=TRUE),
                collapse="")
  dna <- Biostrings::DNAStringSet(seq)
  genomeDir <- file.path(tempdir(), as.integer(runif(1) * 1000000000))
  if (file.exists(genomeDir)) unlink(genomeDir, recursive=TRUE)
  dir.create(genomeDir, recursive=TRUE)
  on.exit(unlink(genomeDir, recursive=TRUE))
  checkException(GmapGenome(genome=dna, directory=genomeDir,
                            name="thing", create=TRUE))
  names(dna) <- "sampleDNAStringSet"
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

testGmapGenome_constructor_FastaFile_create <- function() {
  fa <- system.file("extdata/hg19.p53.fasta", package="gmapR", mustWork=TRUE)
  fastaFile <- rtracklayer::FastaFile(fa)
  gmapGenome <- GmapGenome(fastaFile, create=TRUE)
  checkTrue(is(gmapGenome, "GmapGenome"))
}

test_GmapGenome_accessors <- function() {
  genomeName <- "testGenome"
  dna <- Biostrings::DNAStringSet("ACTGTGTCAG")
  names(dna) <- "testDNAString"
  genomeDir <- file.path(tempdir(), as.integer(runif(1) * 1000000000))
  if (file.exists(genomeDir)) unlink(genomeDir, recursive=TRUE)
  dir.create(genomeDir, recursive=TRUE)
  on.exit(unlink(genomeDir, recursive=TRUE))
  gmapGenome <- GmapGenome(genome=dna, directory=genomeDir,
                           name=genomeName, create=FALSE)
  checkIdentical(path(gmapGenome), file.path(genomeDir, genomeName))
  checkTrue(is(directory(gmapGenome), "GmapGenomeDirectory"))
  checkIdentical(genome(gmapGenome), genomeName)
}
