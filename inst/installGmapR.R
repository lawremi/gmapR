library(gmapR)

##warning: must have "~/.gmapR_globals.dcf" created and gsnap_save_dir
##along with gsnap_bin_dir set to location specified in this script
##(by gmap_install_dir variable).  This script will verify this
##constraint is met.  Specifically, if gmap_install_dir is
##gne/home/coryba/sandbox/gmapTestInstallDir, then in
##"~/.gmapR_globals.dcf" gsnap_bin_dir must be
##/gne/home/coryba/sandbox/gmapTestInstallDir/bin and gsnap_save_dir
##must be /gne/home/coryba/sandbox/gmapTestInstallDir/share
##
##Also requires samtools already installed and
##correctly pointed to in "~/.gmapR_globals.dcf".

##path to samtools source folder
samtools_dir <- '/gne/home/coryba/source/samtools/samtools-0.1.16/samtools'

##gmap_install_dir <- dirname(globals()$gsnap_bin_dir)
##for testing
gmap_install_dir <- '/gne/home/coryba/sandbox/gmapTestInstallDir'

if(gmap_install_dir != dirname(dirname(globals()[['gsnap']])))
  stop("settings in ~/.gmapR_globals.dcf do not match specified installation directory")

##build gmap. Not needed if gmap already built at gmap_install_dir
installGmap(gmap_install_dir, samtools_dir)

##build hg19 IIT
buildGmapIndex('hg19', gmap_data_dir=globals()[['gsnap_save_dir']])
##build dnsnp131
buildGmapDbSNPIndex('snp131')
##build splice sites
buildRefseqSpliceSite('hg19')
