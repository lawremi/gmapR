library(gmapR)

test_bam_tally <- function() {
  param <- BamTallyParam(TP53Genome(), TP53Which(), indels = TRUE)
  bam <- system.file("extdata/H1993.analyzed.bam", 
                     package="LungCancerLines", mustWork=TRUE)
  tallies <- bam_tally(bam, param)
  variants <- variantSummary(tallies, NULL, 0L)
}

test_IITs_not_created <- function() {
  gmapGenome <- GmapGenome(genome="NoGenome", directory = getwd())
  btp <- BamTallyParam(genome=gmapGenome,
                       which = RangesList())
  bam_file <-
    system.file("extdata/test_data_aln/test_data_aln.concordant_uniq.bam",
                package="gmapR", mustWork=TRUE)
  checkException(bam_tally(x=bam_file, param=btp))
}
