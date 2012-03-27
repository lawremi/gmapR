##############################################################
### Creates a GRanges of chromosome length from the gmap genome
##############################################################

getChrGRangesFromGmap <- function(genome, genome_dir) {
  tab <-read.table(pipe(paste("get-genome -L -d", genome, "-D", genome_dir)))
  chr_gr <- GRanges(seqnames = tab[,1],
                    IRanges(start = 1,
                            end = as.numeric(tab[,3])))
  seqlengths(chr_gr)[tab[,1]] <- as.numeric(tab[,3])
  return(chr_gr)
}
