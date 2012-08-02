test_snps <- function() {
  ## FIXME: rewrite this test to use the output of bam_tally, once we
  ## can convert those results to a VCF.
  library(gmapR)
  genome <- GmapGenome("hg19_IGIS21", create = TRUE)
  dir <- GmapSnpDirectory(genome)
  gr <- as(seqinfo(genome), "GenomicRanges")
  gr <- GRanges("3", IRanges(3e6,4e6))
  options(error=recover)
  snps(dir, genome = genome, which = gr) <- "~/share/tracks/00-All.vcf.gz"
  snps(genome, which = gr, name = "snps2") <-
    "~/share/tracks/00-All.vcf.gz"
}
