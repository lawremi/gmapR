##' Downloads a version of dbSNP and integrates it into an installation of GMAP/GSNAP for SNP-tolerant alignment
##'
##' .. content for \details{} ..
##' @title Download and Integrate dbSNP into GMAP/GSNAP
##' @param version version of dbSNP to obtain
##' @return 0
##' @author Cory Barr
##' @export
buildGmapDbSNPIndex <- function(version) {

  if(version != 'snp131')
    stop("only supports dbSNP version snp131 at present")

  tmp_dir <- 'snp.tmp'
  if(file.exists(tmp_dir))
    unlink(tmp_dir, recursive=TRUE)
  dir.create(tmp_dir)
  
  start_dir <- file_path_as_absolute(getwd())
  setwd(tmp_dir)
  system('wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/snp131.txt.gz')
  system('gunzip snp131.txt.gz')
  
  ##all entries must be in this format: >rs62211261 21:14379270..14379270 CG
  converted_script <- system.file("scripts/parse_dbsnp.pl", package="gmapR")
  sys_call <- paste("cat",
                    "snp131.txt",
                    "| perl",
                    converted_script,
                    ">",
                    "snp131.converted")
  system(sys_call)

  ##cat snp131.converted | iit_store -o snp131
  iit_store <- file.path(start_dir,
                         "gmap/bin/iit_store")
  sys_call <- paste("cat",
                    "snp131.converted",
                    "|",
                    iit_store,
                    "-o",
                    "snp131")
  system(sys_call)

  ##cp -f snp131.iit /gne/home/coryba/tools/gsnap/share/hg19_ucsc/hg19_ucsc.maps
  sys_call <- paste("cp -f",
                    "snp131.iit",
                    file.path(start_dir,
                              "gmap",
                              "genomes",
                              "hg19",
                              "hg19.maps"))
  res <- system(sys_call, intern=TRUE)

  ##snpindex -d hg19_ucsc -V /gne/home/coryba/tools/gsnap/share/hg19_ucsc -v snp131 snp131.iit
  snpindex <- file.path(start_dir, "gmap/bin/snpindex")
  sys_call <- paste(snpindex,
                    "-d",
                    "hg19",
                    "-V",
                    file.path(start_dir,
                              "gmap/genomes",
                              "hg19"),
                    "-v snp131 snp131.iit")
  res <- system(sys_call, intern=TRUE)

  setwd(start_dir)
  unlink(tmp_dir, recursive=TRUE)

  return(0)
}
