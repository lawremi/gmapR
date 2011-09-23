library(gmapR)

buildConfigFile <- function() {
  gmap_bin <- system("which gmap", intern=TRUE)
  if(length(gmap_bin) ==0)
    stop("Did not find a gmap installation in $PATH")
  
  gsnap_bin_dir <- dirname(gmap_bin)
  
  samtools <- system("which samtools", intern=TRUE)
  ##DEVEL:
  warning("remove the next line..for devel only")
  samtools <- '/gne/home/coryba/bin/samtools'
  if (length(samtools) == 0)
    stop("Did not find samtools installation. Exiting.")
  
  ##check to see if gmap built with samtools
  if (!file.exists(file.path(gsnap_bin_dir, "bam_tally")))
    stop("Did not find bam_tally. Gmap likely not compiled with samtools source. Exiting.")
  
  gsnap_save_dir <- grep("^Default gmap directory \\(environment\\)",
                         system("gmap --version", intern=TRUE),
                         value=TRUE)
  gsnap_save_dir <- tail(unlist(strsplit(gsnap_save_dir, " ")), n=1)
  
  gmapR_config_file <- file.path(system.file("", package="gmapR"),
                                 "gmapR_globals-default.dcf")
  config_lines <- c(paste("samtools:", samtools),
                    paste("gsnap_bin_dir:", gsnap_bin_dir),
                    paste("gsnap_save_dir:", gsnap_save_dir))
  writeLines(text=config_lines, con=gmapR_config_file)
}

buildConfigFile()

##build hg19 IIT
buildGmapIndex('hg19', gmap_data_dir=globals()[['gsnap_save_dir']])
##build dnsnp131
buildGmapDbSNPIndex('snp131')
##build splice sites
buildRefseqSpliceSites('hg19')
