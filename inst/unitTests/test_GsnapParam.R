test_GsnapParam_constructor <- function() {
  
  fa <- system.file("extdata/hg19.p53.fasta", package="gmapR", mustWork=TRUE)
  fastaFile <- rtracklayer::FastaFile(fa)
  gmapGenome <- GmapGenome(fastaFile, create=TRUE, k = 12)

  gsnapParam <- GsnapParam(genome = gmapGenome,
                           unique_only = FALSE,
                           max_mismatches = NULL,
                           suboptimal_levels = 0L, mode = "standard",
                           npaths = 10L,
                           novelsplicing = FALSE, splicing = NULL, 
                           nthreads = 1L,
                           batch = "2")

  checkTrue(is(gsnapParam, "GsnapParam"))
}
