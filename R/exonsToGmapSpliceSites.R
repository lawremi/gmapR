##' Convert a data frame of exon information into input to build a splice site index for Gmap/Gsnap
##'
##' the input data frame needs to have columns named "chrom,"
##' "strand," "exonStarts", "exonEnds", and "name."  A common way to
##' obtain this is from UCSC via rtracklayer:
##'
##' library(rtracklayer)
##' session <- browserSession()
##' genome(session) <- genome
##' query <- ucscTableQuery(session, "refGene")
##' refGeneTable <- getTable(query)
##'
##' To create a splice-site index for Gmap/Gsnap, take the following
##' output (here called spliceSites) and pass it to buildGmapSpliceSites:
##'
##' spliceSites <- exonsToSpliceSites(refGeneTable)
##' 
##' @title Convert a data frame of exon info into input for incorporating known splice sites into Gmap/Gsnap
##' @param exons A data frame containing the fields described in the details section.
##' @return a character vector
##' @author Cory Barr
##' @export
exonsToGmapSpliceSites <- function(exons) {

  if(class(exons) != "data.frame")
    stop("input is not a data.frame")
  
  gene_names <- as.character(exons$name)
  chroms <- as.character(exons$chrom)
  strands <- as.character(exons$strand)
  exon_lcoords <- as.character(exons$exonStarts)
  exon_rcoords <- as.character(exons$exonEnds)
  
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
  
  exon_lcoords <- lapply(exon_lcoords, .positions_to_int)
  exon_rcoords <- lapply(exon_rcoords, .positions_to_int)
  if(any(elementLengths(exon_lcoords) != elementLengths(exon_rcoords)))
    stop("not the same number of exon starts and stops")

  single_exon_genes <- which(elementLengths(exon_lcoords) < 2)
  if(length(single_exon_genes) > 0) {
    exon_lcoords <- exon_lcoords[-single_exon_genes]
    exon_rcoords <- exon_rcoords[-single_exon_genes]
    gene_names <- gene_names[-single_exon_genes]
    chroms <- chroms[-single_exon_genes]
    strands <- strands[-single_exon_genes]
  }
  rm(single_exon_genes)

  .positions_to_gsnap_input <- function(positions, gene_name, chr, type, strand) {
    strand_str <- 'pos'
    if (strand == '-')
      strand_str <- 'neg'

    exon_nums <- seq_len(length(positions))    
    if (type == 'acceptor')
      exon_nums <- exon_nums + 1

    paste(">", gene_name, ".", strand_str, ".", "exon", exon_nums, " ", chr,
          ":", .int_to_half_open(positions, strand), " ", type, sep="")
  }
  
  .getGmapInputStrings <- function(i) {
    lcoords <- unlist(exon_lcoords[i])
    rcoords <- unlist(exon_rcoords[i])
    strand <- strands[i]
    gene_name <- gene_names[i]
    chrom <- chroms[i]
    num_exons <- length(lcoords)
    
    if (strand == '+') {
      donor_pos <- rcoords[1:(num_exons - 1)]
      acceptor_pos <- lcoords[2:num_exons]
    } else if (strand == '-') {
      donor_pos <- lcoords[num_exons:2]
      acceptor_pos <- rcoords[(num_exons - 1):1]
    } else {
      stop("Strand is not + or -")
    }

    ret_val <- list()
    ret_val$donor_str <- .positions_to_gsnap_input(donor_pos, gene_name, chrom, 'donor', strand)
    ret_val$acceptor_str <- .positions_to_gsnap_input(acceptor_pos, gene_name, chrom, 'acceptor', strand)
    return(ret_val)
  }

  if(is.loaded("mc_fork", PACKAGE = "multicore")) {
    apply_func <- mclapply
  } else {
    apply_func <- lapply
  } 
  gmapInputStrings <- apply_func(seq_len(length(exon_lcoords)), .getGmapInputStrings)
  
  .pairUpResults <- function(l) {      
    donors <- l$donor_str
    acceptors <- l$acceptor_str
    
    paired <- rep(as.character(NULL), length(donors) * 2)
    evens <- 2 * seq_len(length(donors))
    odds <- evens - 1
    paired[odds] <- donors
    paired[evens] <- acceptors
    return(paired)      
  }
  
  paired_up <- apply_func(gmapInputStrings, .pairUpResults)
  return(unlist(paired_up))
}
