buildGmapSpliceSites <- function(splice_sites,
                                 gsnap_data_dir,
                                 genome,
                                 splices_name) {
  scratch_dir <- file.path(tempdir(), "scratch_dir")
  cur_dir <- getwd()
  on.exit(setwd(cur_dir))
  if(!file.exists(scratch_dir)) dir.create(scratch_dir)

  splice_file <- file.path(scratch_dir, "splice_file.txt")
  writeLines(splice_sites, con=splice_file)

  iit_store <- file.path(globals()$gsnap_bin_dir, "iit_store")
  splice_site_file <- paste(splices_name, ".splices", sep="")
  setwd(scratch_dir)
  sys_command <- paste("cat", splice_file, "|", iit_store, "-o", splice_site_file)
  if(system(sys_command) != 0)
    stop(paste("Error executing system command:", sys_command))

  iit_dest_dir <- file.path(gsnap_data_dir,
                            genome,
                            paste(genome, ".maps", sep=""))

  iit_file <- paste(splices_name, '.iit', sep="")
  
  sys_command <- paste("cp", iit_file, iit_dest_dir)
  res <- system(sys_command)
  if(res != 0)
    stop(paste("Error executing system command:", sys_command))
  return(0)
}
