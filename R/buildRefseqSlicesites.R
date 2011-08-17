library(rtracklayer)
library(HTSeqGenieBase)
##' add refSeq splice sites to gsnap
##'
##' downloads refSeq table from UCSC, parses out splice site donors
##' and acceptors, creates IIT for gsnap, then copies IIT to genome's
##' shared gsnap dir
##' @title Build RefSeq Splice Site IIT for gsnap
##' @param genome Genome to construct split sites for (only hg19 supported at the moment)
##' @return 1
##' @author Cory Barr
gsnapBuildRefseqSplicesites <- function(genome) {

  if(!(genome == 'hg19' || genome == 'mm9'))
    stop("TODO: only hg19 and mm9 supported. Have to restructure gsnap's shared directory to directly match assembly name")

  session <- browserSession()
  genome(session) <- genome
  query <- ucscTableQuery(session, "refGene")
  refGeneTable <- getTable(query)

  gene_names <- as.character(refGeneTable$name)
  chroms <- as.character(refGeneTable$chrom)
  strands <- as.character(refGeneTable$strand)
  exon_starts <- as.character(refGeneTable$exonStarts)
  exon_ends <- as.character(refGeneTable$exonEnds)
  
  ##exon_starts and exon_ends should be a single string comma-sep set
  ##of coordinates
  .exonPositionsToSpliceSites <- function(index) {

    .positions_to_int <- function(x) {
      as.integer(unlist(strsplit(x, ",")))
    }

    .int_to_half_open <- function(x, strand) {
      if (strand == '+') {
        paste(x, x + 1L, sep='..')
      } else if (strand == '-') {
        paste(x + 1L, x, sep='..')
      }
    }

    .positions_to_gsnap_input <- function(positions, gene_name, chr, type, strand) {
      strand_str <- 'pos'
      if (strand == '-')
        strand_str <- 'neg'
      paste(">",
            gene_names[index],
            ".",
            strand_str,
            ".",
            "exon",
            seq_len(length(positions)),
            " ",
            chroms[index],
            ":",
            .int_to_half_open(positions, strand),
            " ",
            type,
            sep="")
    }
    
    starts <- .positions_to_int(exon_starts[index])
    ends <- .positions_to_int(exon_ends[index])

    if(length(starts) != length(ends))
       stop("There are a different number of exon starts and exon ends")

    if(length(starts) < 2 )
      return(character(0))

    num_exons <- length(starts)
    if(strands[index] == '+') {
      donor_strings <- .positions_to_gsnap_input(ends[1:(num_exons-1)],
                                                 gene_names[index],
                                                 chroms[index],
                                                 "donor",
                                                 strands[index])
      acceptor_strings <- .positions_to_gsnap_input(starts[2:num_exons],
                                                    gene_names[index],
                                                    chroms[index],
                                                    "acceptor",
                                                    strands[index])                                
    } else {
      donor_strings <- .positions_to_gsnap_input(starts[2:num_exons],
                                                 gene_names[index],
                                                 chroms[index],
                                                 "donor",
                                                 strands[index])
      acceptor_strings <- .positions_to_gsnap_input(ends[1:(num_exons-1)],
                                                    gene_names[index],
                                                    chroms[index],
                                                    "acceptor",
                                                    strands[index])                                
    }

    pairUpResults <- function(strVec1, strVec2) {      
      
      paired <- rep(as.character(NULL), length(strVec1) * 2)
      evens <- 2 * seq_len(length(strVec1))
      odds <- evens - 1
      paired[odds] <- strVec1
      paired[evens] <- strVec2
      return(paired)      
    }
    
    paired_up <- pairUpResults(donor_strings, acceptor_strings)
  }

  res <- unlist(mclapply(seq_len(length(gene_names)), .exonPositionsToSpliceSites))
  scratch_dir <- file.path(getwd(), paste('tmp', ceiling(runif(1) * 1000000), sep='.'))
  dir.create(scratch_dir)
  splice_file <- file.path(scratch_dir, "splice_file.txt")

  writeLines(res, splice_file)

  iit_store <- file.path(globals()$gsnap_bin_dir, "iit_store")
  splice_site_file <- 'refSeq.splices'
  setwd(scratch_dir)
  sys_command <- paste("cat", splice_file, "|", iit_store, "-o", splice_site_file)
  system(sys_command)

  if (genome=='hg19') {
    iit_dest_dir <- file.path(globals()$gsnap_save_dir, "hg19_ucsc", "hg19_ucsc.maps") } else if (genome=='mm9') {
      iit_dest_dir <- file.path(globals()$gsnap_save_dir, "mm9", "mm9.maps") } else {
        stop("genome not supported")
    }

  iit_file <- 'refSeq.splices.iit'
  
  sys_command <- paste("cp", iit_file, iit_dest_dir)
  system(sys_command)

  setwd('..')
    
  unlink(scratch_dir, recursive=T)

  return(1)
}
