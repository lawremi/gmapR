TP53Genome <- function() {
  checkPackageInstalled("org.Hs.eg.db", required = TRUE)
  checkPackageInstalled("TxDb.Hsapiens.UCSC.hg19.knownGene", required = TRUE)
  checkPackageInstalled("BSgenome.Hsapiens.UCSC.hg19", required = TRUE)
  gene <- "TP53"
  eg <- org.Hs.eg.db::org.Hs.egSYMBOL2EG[[gene]]
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  tx <- transcripts(txdb, vals = list(gene_id = eg))
  p53Seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, range(tx) + 1e6,
                   as.character = FALSE)
  names(p53Seq) <- gene
  genome <- GmapGenome(genome = p53Seq,
                       name = paste0(gene, "_demo"),
                       create = TRUE)
  spliceSites(genome, "") <- txdb
  genome
}
