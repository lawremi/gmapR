##' This function helps create the files needed for SNP-tolerant gsnap
##' alignment.
##'
##' @title Build SNP-tolerant gsnap alignment
##' @param snp_file The fasta-like file containing the SNPs
##' @param gsnap_bin location of the gsnap executable. Defaults to
##' location at gsnap:::globals()$gsnap
##' @param gsnap_data_dir location of the gmap data directory where indexes are stored. Defaults
##' to location at gsnap:::globals()$gsnap_save_dir
##' @param genome name of the genome, which corresponds to
##' the subdirectory of save_dir. The IIT files go here. Typical
##' examples are "hg19" and "mm9"
##' @param snps_name This argument specifies the name of the
##' argument to call via gsnap's -v option for snp-tolerant alignment. For example, if snp_name
##' is db131, then gsnap will align via gsnap will incorporate your injected snps through the command 'gsnap -d genome -v db131'
##' @return 0
##` @author Cory Barr
##` rdname buildGsnapSNPIIT

buildGsnapSNPIIT <- function(snp_file,
                             gsnap_bin,
                             gsnap_data_dir,
                             genome,
                             snps_name) {

  if (missing(gsnap_bin))
    gsnap_bin <- globals()$gsnap
  if(missing(gsnap_data_dir))
    gsnap_data_dir <- globals()$gsnap_save_dir

  full_genome_dir <- file.path(gsnap_data_dir, genome)
  
  if(!file.exists(snp_file))
    stop(paste("Could not find the file", snp_file))
  if(!file.exists(gsnap_data_dir))
    stop(paste("The directory", gsnap_data_dir, "does not exist."))
  if(!file.exists(full_genome_dir))
    stop(paste("The genome", genome, "is not installed in", gsnap_data_dir))
    
  tmp_dir <- tempdir()
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive=TRUE))

  iit_store <- file.path(globals()$gsnap_bin_dir, "iit_store")

  iit_file <- file.path(tmp_dir, snps_name)
  iit_file <- paste(iit_file, ".iit", sep="")
  
  ##ex: cat snp131.converted | iit_store -o snp131
  sys_call <- paste("cat", snp_file, "|", iit_store, "-o", iit_file)
  res <- system(sys_call)
  if (res != 0) stop("Error performing system call to iit_store")

  ##ex: cp -f snp131.iit /gne/home/coryba/tools/gsnap/share/hg19_ucsc/hg19_ucsc.maps
  maps_dir <- file.path(full_genome_dir,
                        paste(genome, ".maps", sep=''))
  sys_call <- paste("cp -f", iit_file, maps_dir)
  res <- system(sys_call)
  if (res != 0) stop("Error performing system call to copy snp IIT")
  
  ##ex: snpindex -d hg19_ucsc -V /gne/home/coryba/tools/gsnap/share/hg19_ucsc -v snp131 snp131.iit
  snpindex <- file.path(globals()$gsnap_bin, "snpindex")
  sys_call <- paste(snpindex,
                    "-d",
                    genome,
                    "-D",
                    dirname(full_genome_dir),
                    "-V",
                    full_genome_dir,
                    "-v",
                    sub("\\.iit$", "", basename(iit_file)),
                    iit_file)
  res <- system(sys_call)
  if (res != 0) stop("Error performing system call to gsnap")
  
  return(0)
}
