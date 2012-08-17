test_bam_tally <- function() {
  f <- system.file("extdata", "test_data_aln",
                   "test_data_aln.concordant_uniq.bam", package = "gmapR")
  bf <- BamFile(f)
  
  bam_tally(bf, BamTallyParam(genome))
}

test_bam_tally_C <- function() {
  library(BSgenome.Dmelanogaster.UCSC.dm3)

  genome <- GmapGenome(Dmelanogaster, create = TRUE)

  library(pasillaBamSubset)
  bam <- untreated3_chr4()
  
  bf <- Rsamtools::BamFile(bam)
  which <- RangesList(chr4 = IRanges(1e6, 2e6))
  gr <- bam_tally(bf, BamTallyParam(genome))
  gr <- bam_tally(bf, BamTallyParam(genome, variant_strand = 1L))
  gr <- bam_tally(bf, BamTallyParam(genome, which = which, variant_strand = 1L))
  empty <- RangesList(chr2L = IRanges(1e6, 2e6))
  gr <- bam_tally(bf, BamTallyParam(genome, which = empty))
  
  gr <- bam_tally(bf, BamTallyParam(genome, cycle_breaks = c(1, 15, 30, 40)))
}
