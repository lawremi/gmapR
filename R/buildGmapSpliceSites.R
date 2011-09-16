buildGmapSpliceSites <- function(splice_sites,
                                 gsnap_data_dir,
                                 genome,
                                 splices_name) {
  scratch_dir <- file.path(tempdir(), "scratch_dir")
  cur_dir <- getwd()
  on.exit({unlink(scratch_dir, recursive=TRUE)
           setwd(cur_dir)})
  if(!file.exists(scratch_dir)) dir.create(scratch_dir)

  splice_file <- file.path(scratch_dir, "splice_file.txt")
  writeLines(splice_sites, con=splice_file)

  iit_store <- file.path(globals()$gsnap_bin_dir, "iit_store")
  splice_site_file <- paste(splices_name, ".splices", sep="")
  setwd(scratch_dir)
  sys_command <- paste("cat", splice_file, "|", iit_store, "-o", splice_site_file)
  system(sys_command)

  ##TODO: move this to the setup script for HTSeqGenieBase
  ##if (genome=='hg19') {
  ##  iit_dest_dir <- file.path(globals()$gsnap_save_dir, "hg19_ucsc", "hg19_ucsc.maps")
  ##} else if (genome=='mm9') {
  ##    iit_dest_dir <- file.path(globals()$gsnap_save_dir, "mm9", "mm9.maps") }
  ##else {
  ##  stop("genome not supported")
  ##}

  iit_dest_dir <- file.path(gsnap_data_dir,
                            genome,
                            paste(genome, ".maps", sep=""))

  iit_file <- paste(splices_name, '.splices.iit', sep="")
  
  sys_command <- paste("cp", iit_file, iit_dest_dir)
  system(sys_command)
  
  return(0)
}
