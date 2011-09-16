##' Downloads the GMAP source code, builds it, then installs it to a specified directory
##'
##' 
##' @title Download and Install GMAP
##' @param install_dir directory where GMAP will be installed
##' @param samtools_dir directory where samtools is installed
##' @return 0
##' @author Cory Barr
##' @export
installGmap <- function(install_dir, samtools_dir=NULL) {

  if (is.null(samtools_dir))
    warning("Samtools directory not specified. Utilities to parse BAM files will be unavailable")
  
  if (!file.exists(install_dir))
    dir.create(install_dir, recursive=TRUE)
  install_dir <- file_path_as_absolute(install_dir)
  
  gmap_temp_dir <- file.path(tempdir(), "gmap_temp_dir")
  if(file.exists(gmap_temp_dir))
    unlink(gmap_temp_dir, recursive=TRUE)
  dir.create(gmap_temp_dir)
  on.exit(unlink(gmap_temp_dir, recursive=TRUE))
  setwd(gmap_temp_dir)

  gmap_zip <- 'gmap-gsnap-2011-09-09.tar.gz'
  system(paste('wget',
               file.path('http://research-pub.gene.com/gmap/src', gmap_zip)))
  dir.create("unzip_dir")
  system(paste('tar xzvf', gmap_zip, "-C unzip_dir"))  
  orig_dir <- getwd()
  gmap_src_dir <- dir("unzip_dir", full.names=TRUE)
  setwd(gmap_src_dir)
  
  gsnap_dir <- file.path(install_dir, "gmap")
  gsnap_genomes_dir <- file.path(install_dir, "genomes")
  dir.create(gsnap_genomes_dir, recursive=TRUE)

  sys_call <- paste("./configure",
                    paste("prefix=", install_dir, sep=''))#,
                    ##paste("with_gmapdb=", gsnap_genomes_dir, sep=''))
  if (!is.null(samtools_dir))
    sys_call <- paste(sys_call,
                      paste("with_samtools=", samtools_dir, sep=""))
  system(sys_call)                  
  system('make')
  system('make install')

  unlink(gmap_temp_dir, recursive=TRUE)
  return(0)
}
