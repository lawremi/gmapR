##' add refSeq splice sites to gsnap
##'
##' downloads refSeq table from UCSC, parses out splice site donors
##' and acceptors, creates IIT for gsnap, then copies IIT to genome's
##' shared gsnap dir
##' @title Build RefSeq Splice Site IIT for gsnap
##' @param genome Genome to construct split sites for (only hg19 supported at the moment)
##' @return 0
##' @author Cory Barr
##' @export
buildRefseqSpliceSites <- function(genome) {

  if(!(genome == 'hg19' || genome == 'mm9'))
    stop("TODO: only hg19 and mm9 supported. Have to restructure gsnap's shared directory to directly match assembly name")

  print("Downloading gene models from UCSC...")
  session <- browserSession()
  genome(session) <- genome
  query <- ucscTableQuery(session, "refGene")
  refGeneTable <- getTable(query)
  print("...done")
  
  ##exon_starts and exon_ends should be a single string comma-sep set
  ##of coordinates
  splice_site_lines <- exonsToGmapSpliceSites(refGeneTable)

  res <- buildGmapSpliceSites(splice_sites=splice_site_lines,
                              gsnap_data_dir=globals()[['gsnap_save_dir']],
                              genome='hg19',
                              splices_name='refSeq.splices')
  if(res != 0)
    stop("Error building gmap splice sites for hg19")
}
