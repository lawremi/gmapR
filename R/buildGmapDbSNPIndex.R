##' Downloads a version of dbSNP and integrates it into an installation of GMAP/GSNAP for SNP-tolerant alignment
##'
##' .. content for \details{} ..
##' @title Download and Integrate dbSNP into GMAP/GSNAP
##' @param version version of dbSNP to obtain
##' @return 0
##' @author Cory Barr
##' @export
buildGmapDbSNPIndex <- function(version) {

  ##TODO: use the dbsnp Bioc package instead of much of the code below
  
  if(version != 'snp131')
    stop("only supports dbSNP version snp131 at present")

  tmp_dir <- 'snp.tmp'
  if(file.exists(tmp_dir))
    unlink(tmp_dir, recursive=TRUE)
  dir.create(tmp_dir)
  
  start_dir <- file_path_as_absolute(getwd())
  setwd(tmp_dir)
  system('wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/snp131.txt.gz')
  print("unzipping downloaded snp dbSNP data...")
  system('gunzip snp131.txt.gz')
  print("...done")

  ##all entries must be in this format: >rs62211261 21:14379270..14379270 CG
  converted_script <- system.file("scripts/parse_dbsnp.pl", package="gmapR")
  sys_call <- paste("cat",
                    "snp131.txt",
                    "| perl",
                    converted_script,
                    ">",
                    "snp131.converted")
  res <- system(sys_call)
  if (res != 0)
    stop(paste("Error performing system call:", sys_call))

  ##cat snp131.converted | iit_store -o snp131
  iit_store <- file.path(globals()[['gsnap_bin_dir']], "iit_store")
  sys_call <- paste("cat",
                    "snp131.converted",
                    "|",
                    iit_store,
                    "-o",
                    "snp131")
  res <- system(sys_call)
  if (res != 0)
    stop(paste("Error performing system call:", sys_call))

  ##cp -f snp131.iit /gne/home/coryba/tools/gsnap/share/hg19_ucsc/hg19_ucsc.maps
  sys_call <- paste("cp -f",
                    "snp131.iit",
                    file.path(globals()[['gsnap_save_dir']],
                              "hg19",
                              "hg19.maps"))
  res <- system(sys_call)
  if (res != 0)
    stop(paste("Error performing system call:", sys_call))
  
  ##snpindex -d hg19_ucsc -V /gne/home/coryba/tools/gsnap/share/hg19_ucsc -v snp131 snp131.iit
  snpindex <- file.path(globals()[['gsnap_bin_dir']], "snpindex")  
  sys_call <- paste(snpindex,
                    "-d",
                    "hg19",
                    "-V",
                    file.path(globals()[['gsnap_save_dir']],
                              "hg19"),
                    "-v snp131 snp131.iit")
  res <- try(system(sys_call, intern=TRUE))
  if(class(res) == "try-error")
    stop(paste("Error performing system call:", sys_call))
    
  setwd(start_dir)

  return(0)
}
