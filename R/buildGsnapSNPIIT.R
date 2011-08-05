##TODO: put in gsnap package

##' This function helps create the files needed for SNP-tolerant gsnap
##' alignment.
##'
##' @title Build SNP-tolerant gsnap alignment
##' @param index_file The fasta-like file containing the SNPs
##' @param gsnap_bin location of the gsnap executable. Defaults to
##' location at gsnap:::globals()$gsnap
##' @param save_dir location of the default gmap directory. Defaults
##' to location at gsnap:::globals()$gsnap_save_dir
##' @param gsnap_genome_dir name of the genome, which corresponds to
##' the subdirectory of save_dir. The IIT files go here. Typical
##' examples are "hg19" and "mm9"
##' @param snp_file This argument will specify the name of the
##' resulting final IIT of SNPs and consequently the name of the
##' database specified via gsnap's -v option. For example, if snp_file
##' is db131, then gsnap will align via gsnap -d genome -v db131
##' input.fastq
##' @return 0
##` @author Cory Barr
##` rdname buildGsnapSNPIIT

buildGsnapSNPIIT <- function(index_file,
                             gsnap_bin,
                             save_dir,
                             gsnap_genome_dir,
                             snp_file) {

  if(!file.exists(index_file))
    stop(paste("Could not find the file", index_file))

  if (missing(gsnap_bin))
    gsnap_bin <- globals()$gsnap
  if (missing(save_dir))
    save_dir <- globals()$gsnap_save_dir

  full_save_dir <- file.path(save_dir, gsnap_genome_dir)

  if(!file.exists(full_save_dir))
    stop(paste("Error. The directory", full_save_dir, "does not exist."))

  tmp_dir <- paste("tmp.gsnap_index_creation.",
                   ceiling(runif(1) * 10000000),
                   sep="")

  dir.create(tmp_dir)

  iit_store <- file.path(globals()$gsnap_bin_dir, "iit_store")

  if(missing(snp_file)) {
    snp_file <- gsub("^.*/", "", index_file, perl=T)
  }
  snp_file <- file.path(tmp_dir, snp_file)
  iit_file <- paste(snp_file, ".iit", sep="")
  
  ##ex: cat snp131.converted | iit_store -o snp131
  sys_call <- paste("cat", index_file, "|", iit_store, "-o", iit_file)
  system(sys_call)

  maps_dir <- paste(globals()$gsnap_save_dir,
                    "/",
                    gsnap_genome_dir,
                    "/",
                    gsnap_genome_dir,
                    ".maps",
                    sep="")
  ##ex:cp -f snp131.iit /gne/home/coryba/tools/gsnap/share/hg19_ucsc/hg19_ucsc.maps
  sys_call <- paste("cp -f", iit_file, maps_dir)
  system(sys_call)
  
  snpindex <- file.path(globals()$gsnap_bin, "snpindex")
  gsnap_genome_path <- file.path(globals()$gsnap_save_dir, gsnap_genome_dir)
  ##ex:snpindex -d hg19_ucsc -V /gne/home/coryba/tools/gsnap/share/hg19_ucsc -v snp131 snp131.iit
  sys_call <- paste(snpindex,
                    "-d",
                    gsnap_genome_dir,
                    "-V",
                    gsnap_genome_path,
                    "-v",
                    strsplit(snp_file, "/")[[1]][2],
                    paste(snp_file, ".iit", sep=""))
  system(sys_call)
  unlink(tmp_dir, recursive=TRUE)
  return(0)
}
