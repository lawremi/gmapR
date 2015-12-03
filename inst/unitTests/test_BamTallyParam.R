test_BamTallyParam <- function() {
  options(warn = 2)
  on.exit(options(warn = 0))
  genome <- TP53Genome()
  genome.name <- genome(genome)
  param <- new("BamTallyParam", genome = genome,
               which = gmapR:::normArgWhich(GRanges(), genome),
               minimum_mapq = 0L,
               concordant_only = FALSE, unique_only = FALSE,
               primary_only = FALSE, ignore_duplicates = FALSE,
               min_depth = 0L, variant_strand = 0L,
               ignore_query_Ns = FALSE,
               indels = FALSE, include_soft_clips = 0L,
               xs = FALSE, read_pos = FALSE, min_base_quality = 0L,
               noncovered = FALSE, nm = FALSE)
  which <- TP53Which()
  wicked.which <- renameSeqlevels(which, c(TP53 = "chr1"))
  
  checkException(BamTallyParam(), silent = TRUE)
  checkException(BamTallyParam(5), silent = TRUE)
  
  checkIdentical(BamTallyParam(genome.name), param)
  checkException(BamTallyParam(genome, IRanges(1, 10)), silent = TRUE)
  checkException(BamTallyParam(genome, wicked.which), silent = TRUE)
  param@which <- gmapR:::normArgWhich(which, genome)
  checkIdentical(BamTallyParam(genome, which), param)
  checkException(BamTallyParam(genome, concordant_only = NA), silent = TRUE)
  checkException(BamTallyParam(genome, unique_only = rep(FALSE, 2)),
                 silent = TRUE)
  checkException(BamTallyParam(genome, min_depth = -1), silent = TRUE)
  checkException(BamTallyParam(genome, variant_strand = 5), silent = TRUE)
  param@variant_strand <- 2L
  checkIdentical(BamTallyParam(genome, which, variant_strand = 2), param)
}
