##' add refSeq splice sites to gsnap
##'
##' downloads refSeq table from UCSC, parses out splice site donors
##' and acceptors, creates IIT for gsnap, then copies IIT to genome's
##' shared gsnap dir
##' @title Build RefSeq Splice Site IIT for gsnap
##' @param genome Genome to construct split sites for (only hg19 supported at the moment)
##' @return 0
##' @author Cory Barr
gsnapBuildRefseqSplicesites <- function(genome) {

  if(!(genome == 'hg19' || genome == 'mm9'))
    stop("TODO: only hg19 and mm9 supported. Have to restructure gsnap's shared directory to directly match assembly name")

  session <- browserSession()
  genome(session) <- genome
  query <- ucscTableQuery(session, "refGene")
  refGeneTable <- getTable(query)

  
  ##exon_starts and exon_ends should be a single string comma-sep set
  ##of coordinates
  splice_site_lines <- exonsToGmapSpliceSites(refGeneTable)

  scratch_dir <- file.path(getwd(), paste('tmp', ceiling(runif(1) * 1000000), sep='.'))
  dir.create(scratch_dir)
  splice_file <- file.path(scratch_dir, "splice_file.txt")

  writeLines(splice_site_lines, splice_file)

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

  return(0)
}
