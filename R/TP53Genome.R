TP53Genome <- function() {
  gene <- "TP53"
  genomeName <- paste0(gene, "_demo")
  
  if (genomeName %in% genome(GmapGenomeDirectory())) {
    GmapGenome(genomeName)
  } else {
    checkPackageInstalled("org.Hs.eg.db", required = TRUE)
    checkPackageInstalled("TxDb.Hsapiens.UCSC.hg19.knownGene", required = TRUE)
    checkPackageInstalled("BSgenome.Hsapiens.UCSC.hg19", required = TRUE)
    eg <- org.Hs.eg.db::org.Hs.egSYMBOL2EG[[gene]]
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    tx <- transcripts(txdb, vals = list(gene_id = eg))
    roi <- range(tx) + 1e6
    
    p53Seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, roi,
                     as.character = FALSE)
    names(p53Seq) <- gene
    genome <- GmapGenome(genome = p53Seq,
                         name = genomeName,
                         create = TRUE)
    
    exons <- shift(subsetByOverlaps(exonsBy(txdb), roi), 1L - start(roi))
    exons <- renameSeqlevels(exons, setNames(gene, seqnames(roi)))
    spliceSites(genome, "knownGene") <- exons
    
    genome
  }
}
