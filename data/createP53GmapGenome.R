p53Genome <- local({
  library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(gmapR)
  
  gene <- "TP53"
  eg <- org.Hs.egSYMBOL2EG[[gene]]
  tx <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene,
                    vals = list(gene_id = eg))
  p53Seq <- getSeq(Hsapiens, range(tx) + 1e6, as.character = FALSE)
  names(p53Seq) <- gene
  GmapGenome(genome = p53Seq,
             name = paste0(gene, "_demo"),
             create = TRUE)
})
