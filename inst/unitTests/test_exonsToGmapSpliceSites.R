test_exonsToGmapSpliceSites <- function() {

  exon_df <- data.frame(name=c("gene1", "gene2"),
                        strand=c("+", "-"),
                        chrom=c("chrA", "chrB"),
                        exonStarts=c(
                          paste(c(1, 10, 20), collapse=","),
                          paste(c(5, 150, 5000), collapse=",")),
                        exonEnds=c(
                          paste(c(5, 15, 150), collapse=","),
                          paste(c(51, 170, 5005), collapse=",")))

  expected_res <- c(">gene1.pos.exon1 chrA:5..6 donor",
                    ">gene1.pos.exon2 chrA:10..11 acceptor",    
                    ">gene1.pos.exon2 chrA:15..16 donor",       
                    ">gene1.pos.exon3 chrA:20..21 acceptor",    
                    ">gene2.neg.exon1 chrB:5001..5000 donor",   
                    ">gene2.neg.exon2 chrB:171..170 acceptor",
                    ">gene2.neg.exon2 chrB:151..150 donor",   
                    ">gene2.neg.exon3 chrB:52..51 acceptor")

  actual_res <- exonsToGmapSpliceSites(exon_df)

  checkIdentical(expected_res, actual_res)
}
