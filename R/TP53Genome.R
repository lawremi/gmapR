TP53Genome <- function() {
  gene <- "TP53"
  genomeName <- paste0(gene, "_demo")
  
  if (genomeName %in% genome(GmapGenomeDirectory(create=TRUE))) {
    GmapGenome(genomeName)
  } else{ 
    checkPackageInstalled("org.Hs.eg.db", required = TRUE)
    checkPackageInstalled("TxDb.Hsapiens.UCSC.hg19.knownGene", required = TRUE)
    checkPackageInstalled("BSgenome.Hsapiens.UCSC.hg19", required = TRUE)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene

    ## get region of interest
    roi <- getGeneRoi(gene)
    
    p53Seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, roi,
                     as.character = FALSE)
    names(p53Seq) <- gene
    genome <- GmapGenome(genome = p53Seq,
                         name = genomeName,
                         create = TRUE)
    
    exons <- subsetRegion(exonsBy(txdb), roi, gene)
    spliceSites(genome, "knownGene") <- exons
    
    genome
  }
}

getGeneRoi <- function(gene, extend=1e6) {
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  eg <- org.Hs.eg.db::org.Hs.egSYMBOL2EG[[gene]]
  tx <- transcripts(txdb, vals = list(gene_id=eg))
  range(tx) + extend
}

subsetRegion <- function(x, roi, newseqname) {
  x <- shift(subsetByOverlaps(x, roi), 1L - start(roi))
  x <- renameSeqlevels(x, setNames(newseqname, seqnames(roi)))
  seqlengths(x)[newseqname] <- width(roi)
  x
}
