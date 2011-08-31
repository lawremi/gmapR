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

  if(buildGmapSpliceSites(splice_sites) != 0)
    stop("Could not build Gmap splice site index")
  
  return(0)
}
