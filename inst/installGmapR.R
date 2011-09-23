library(gmapR)

gmap_bin <- system("which gmap", intern=TRUE)
if(length(gmap_bin) ==0)
  stop("Did not find a gmap installation in $PATH")

gsnap_bin_dir <- dirname(gmap_bin)




##check to see if gmap built with samtools
if (!file.exists(file.path(gsnap_bin_dir, "bam_tally")))
  stop("Did not find bam_tally. Gmap likely not compiled with samtools source. Exiting.")

gmapR_config_file <- 


##build gmap. Not needed if gmap already built at gmap_install_dir
##installGmap(gmap_install_dir, samtools_dir)

##build hg19 IIT
buildGmapIndex('hg19', gmap_data_dir=globals()[['gsnap_save_dir']])
##build dnsnp131
buildGmapDbSNPIndex('snp131')
##build splice sites
buildRefseqSpliceSites('hg19')
