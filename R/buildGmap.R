##' Downloads the GMAP source code, builds it, then installs it to a specified directory
##'
##' 
##' @title Download and Install GMAP
##' @param install_dir directory where GMAP will be installed
##' @return 0
##' @author Cory Barr
##' @export
installGmap <- function(install_dir) {

  install_dir <- file_path_as_absolute(install_dir)
  
  temp_dir <- 'temp'
  if(file.exists(temp_dir))
    unlink(temp_dir, recursive=TRUE)
  dir.create(temp_dir)
  setwd(temp_dir)
  
  system('wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2011-03-28.tar.gz')
  system('tar xzvf gmap-gsnap-2011-03-28.tar.gz')
  orig_dir <- getwd()
  on.exit(setwd(orig_dir))
  setwd('gmap-2011-03-28/')
  
  gsnap_dir <- file.path(install_dir, "gmap")
  gsnap_genomes_dir <- file.path(install_dir, "gmap", "genomes")
  dir.create(gsnap_genomes_dir, recursive=TRUE)
  
  sys_call <- paste("./configure",
                    paste("prefix=", gsnap_dir, sep=''),
                    paste("with_gmapdb=", gsnap_genomes_dir, sep=''))
  system(sys_call)                  
  system('make')
  system('make install')

  unlink(temp_dir, recursive=TRUE)
  return(0)
}
