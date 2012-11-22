test_IITs_not_created <- function() {
  gmapGenome <- GmapGenome(genome="NoGenome", directory = getwd())
  btp <- BamTallyParam(genome=gmapGenome,
                       which = RangesList())
  bam_file <- system.file("extdata/test_data_aln/test_data_aln.concordant_uniq.bam", package="gmapR", mustWork=TRUE)
  checkException(bam_tally(x=bam_file, param=btp))
}
