test_GmapGenome_constructor_DNAStringSet_noCreate <- function() {
  dna <- Biostrings::DNAStringSet("ACTGTGTCAG")
  names(dna) <- "test"
  gmapGenome <- GmapGenome(genome=dna, name="thing", create=FALSE, k = 12)
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
                            name="thing", create=TRUE, k=12))
  names(dna) <- "sampleDNAStringSet"
  gmapGenome <- GmapGenome(genome=dna, directory=genomeDir,
                           name="thing", create=TRUE, k=12)
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
                           name=genomeName, create=TRUE, k=12) 
  checkTrue(is(gmapGenome, "GmapGenome"))
}

testGmapGenome_constructor_FastaFile_create <- function() {
  fa <- system.file("extdata/hg19.p53.fasta", package="gmapR", mustWork=TRUE)
  fastaFile <- rtracklayer::FastaFile(fa)
  gmapGenome <- GmapGenome(fastaFile, create=TRUE, k=12)
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
  genomeDir <- normalizePath(genomeDir)
  gmapGenome <- GmapGenome(genome=dna, directory=genomeDir,
                           name=genomeName, create=FALSE, k=12)
  checkIdentical(path(gmapGenome), file.path(genomeDir, genomeName))
  checkTrue(is(directory(gmapGenome), "GmapGenomeDirectory"))
  checkIdentical(genome(gmapGenome), genomeName)
}

test_GmapGenome_spliceSites_replacement <- function() {
  library("TxDb.Hsapiens.UCSC.hg19.knownGene")
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  getTP53Range <- function() {
    library(org.Hs.eg.db)
    eg <- org.Hs.egSYMBOL2EG[["TP53"]]
    txTP53 <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene,
                          filter = list(gene_id = eg))
    rngs <- GRanges(ranges=IRanges(start(range(txTP53)), end(range(txTP53))),
                    seqnames="chr17")
    rngs + 1e6
  }
  rngTP53 <- getTP53Range()
  
  exonsByTx <- exonsBy(txdb, by="tx")
  exonsInRegion <- exonsByTx[exonsByTx %over% rngTP53]
  
  ##shift coords of retrieved exons so the ranges match the 
  ##region of the genome used for this example
  shiftCoords <- function(x) {
    x <- exonsInRegion
    w <- width(x)
    r <- ranges(x)
    r <- r + start(rngTP53)
    width(r) <- w
    ranges(x) <- r
    return(x)
  }
  shiftedExons <- shiftCoords(exonsInRegion)
  genome <- TP53Genome()
  x <- spliceSites(genome, name="dbSnp") <- shiftedExons
  checkIdentical(class(x), class(GRangesList()))
}

test_if_GmapGenome_dir_does_not_exist <- function() {
  checkException(GmapGenome(genome="NoGenome", directory = file.path(tempdir(), "DoesNotExist")))
}

test_GmapGenome_getSeq <- function() {
  genomeName <- "testGenome"
  dna <- Biostrings::DNAStringSet(c(testA = "ACTGTGTCAGTTCATGGGACCGTTGC",
                                    testB = "CAACAAATCCGGG"))
  genomeDir <- tempfile()
  if (file.exists(genomeDir))
    unlink(genomeDir, recursive=TRUE)
  dir.create(genomeDir, recursive=TRUE)
  on.exit(unlink(genomeDir, recursive=TRUE))
  genome <- GmapGenome(genome=dna, directory=genomeDir,
                       name=genomeName, create=TRUE, k=12)
  
  gr <- GRanges(rep(c("testA", "testB"), 2),
                IRanges(c(1, 5, 5, 11), c(1, 10, 7, 10)),
                c("+", "-", "*", "+"))
  seqs <- getSeq(genome, gr)
  checkIdentical(seqs, c("A", "GGATTT", "TGT", ""))

  seqs <- getSeq(genome, GRanges())
  checkIdentical(seqs, character())
}
