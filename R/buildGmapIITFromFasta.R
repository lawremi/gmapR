##' Given a FASTA file, will make the FASTA a reference sequence for GMAP/GSNAP
##'
##' 
##' @title Build GMAP/GSNAP Index from a FASTA file
##' @param genome the name of the genome you want to pass to gsnap to align against your FASTA
##' @param fasta A fasta of the reference sequence you wish to align against
##' @param gmap_data_dir The location to save the GMAP/GSNAP index into
##' @return 0 if successful
##' @author Cory Barr
buildGmapIITFromFasta <- function(genome, fasta, gmap_data_dir=NULL) {

  if(!file.exists(gmap_data_dir))
    dir.create(gmap_data_dir, recursive=TRUE)
  if(!is.null(gmap_data_dir))  gmap_data_dir <- file_path_as_absolute(gmap_data_dir)
  tmp_dir <- tempdir()
  dir.create(tmp_dir)
  cur_wd <- getwd()
  on.exit({unlink(tmp_dir, recursive=TRUE)
          setwd(cur_wd)})
  setwd(tmp_dir)
  
  gmap_setup <- file.path(globals()['gmap_bin'],
                          "gmap_setup")  
  sys_command <- paste(gmap_setup,
                       "-d",
                       genome)
  if(!is.null(gmap_data_dir))
    sys_command <- paste(sys_command, paste("-D",
                                            file.path(gmap_data_dir, genome)))
  sys_command <- paste(sys_command, fasta)
  system(sys_command)
  
  sys_command <- paste("make -f",
                       paste('Makefile.', genome, sep=''),
                       "coords")
  system(sys_command)
  
  sys_command <- paste("make -f",
                       paste('Makefile.', genome, sep=''),
                       "gmapdb")
  system(sys_command)
  
  sys_command <- paste("make -f",
                       paste('Makefile.', genome, sep=''),
                       "install")
  system(sys_command)
}
