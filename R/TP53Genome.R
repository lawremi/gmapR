geneGenomeName <- function(gene) {
  paste0(gene, "_demo_", packageVersion("TxDb.Hsapiens.UCSC.hg19.knownGene"))
}

TP53Genome <- function() {
  gene <- "TP53"
  genomeName <- geneGenomeName(gene)
  
  if (genomeName %in% genome(GmapGenomeDirectory(create=TRUE))) {
    GmapGenome(genomeName)
  } else{ 
    checkPackageInstalled("TxDb.Hsapiens.UCSC.hg19.knownGene")
    checkPackageInstalled("BSgenome.Hsapiens.UCSC.hg19")
    checkPackageInstalled("org.Hs.eg.db")    
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene

    ## get region of interest
    roi <- getGeneRoi(txdb, org.Hs.eg.db::org.Hs.eg.db, gene)

    strand(roi) <- "+"
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

exonsToGene <- range

getExons <- function(txdb, orgdb, gene) {
  eg <- select(orgdb, gene, "ENTREZID", "SYMBOL")$ENTREZID
  exons(txdb, vals = list(gene_id=eg))
}

getGeneRoi <- function(txdb, orgdb, gene, extend=1e6) {
  exons <- getExons(txdb, orgdb, gene)
  exonsToGene(exons) + extend
}

subsetRegion <- function(x, roi, newseqname) {
  if (all(!grepl("^chr", seqlevels(x))))
    seqlevels(roi) <- sub("chr", "", seqlevels(roi))
  x <- shift(subsetByOverlaps(x, roi, ignore.strand=TRUE), 1L - start(roi))
  x <- renameSeqlevels(x, setNames(newseqname, seqnames(roi)))
  x <- keepSeqlevels(x, newseqname)
  seqlengths(x) <- width(roi)
  x
}

translateToP53Genome <- function(x) {
  checkPackageInstalled("TxDb.Hsapiens.UCSC.hg19.knownGene")
  checkPackageInstalled("org.Hs.eg.db")    
  gene <- "TP53"
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  orgdb <- org.Hs.eg.db::org.Hs.eg.db
  roi <- getGeneRoi(txdb, orgdb, "TP53")
  subregion <- subsetRegion(x, roi, gene)
  genome(subregion) <- geneGenomeName("TP53")
  subregion
}

exonsOnTP53Genome <- function(gene) {
  checkPackageInstalled("TxDb.Hsapiens.UCSC.hg19.knownGene")
  checkPackageInstalled("org.Hs.eg.db")    
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  orgdb <- org.Hs.eg.db::org.Hs.eg.db
  translateToP53Genome(getExons(txdb, orgdb, gene))
}

geneOnTP53Genome <- function(gene) {
  exonsToGene(exonsOnTP53Genome(gene))
}

TP53Which <- function() {
  ##geneOnTP53Genome("TP53")
  ## for performance:
  gr <- GRanges("TP53", IRanges(1000001, 1025767))
  seqinfo(gr) <- seqinfo(TP53Genome())
  gr
}


