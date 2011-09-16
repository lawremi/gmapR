library(gmapR)

##path to samtools source folder
samtools_dir <- '/gne/home/coryba/source/samtools/samtools-0.1.16/samtools'

##HTSeqGenieBase used to find location of gsnap binary. gmapR package should know its own location
library(HTSeqGenieBase)
gmap_install_dir <- dirname(HTSeqGenieBase::globals()$gsnap_bin_dir)

##for testing
##gmap_install_dir <- '/gne/home/coryba/sandbox/gmapTestInstallDir'

##build gmap. Not needed if gmap already built at gmap_install_dir
installGmap(gmap_install_dir, samtools_dir)

##NOTE: bin gets appended to gmap_install_dir, need to fix so it
##actually installs where it's supposed to. That's why dirname is
##called on globals()$gnap_bin_dir.

##build hg19 IIT
buildGmapIndex('hg19')
##build dnsnp131
buildGmapDbSNPIndex('hg19')
##build splice sites
buildRefseqSpliceSite('hg19')
