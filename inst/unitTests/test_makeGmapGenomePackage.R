test_makeGmapGenomePackage_sacCer3 <- function() {
  library(BSgenome.Scerevisiae.UCSC.sacCer3)
     
  gmapGenomePath <- file.path(tempdir(), as.integer(runif(1) * 100000))
  on.exit(unlink(gmapGenomePath, recursive=TRUE))
  if (file.exists(gmapGenomePath)) unlink(gmapGenomePath, recursive=TRUE)
  ggd <- GmapGenomeDirectory(gmapGenomePath, create = TRUE)
  gmapGenome <- GmapGenome(genome=Scerevisiae,
                           directory = ggd,
                           name = "yeast",
                           create = TRUE)

  packageDestDir <- file.path(tempdir(), as.integer(runif(1) * 100000))
  on.exit(unlink(packageDestDir, recursive=TRUE), add=TRUE)  
  success <- makeGmapGenomePackage(gmapGenome=gmapGenome,
                                   version="0.1.0",
                                   maintainer="<your.name@somewhere.com>",
                                   author="Your Name",
                                   destDir=packageDestDir,
                                   license="Artistic-2.0",
                                   pkgName="GmapGenome.Scerevisiae.UCSC.sacCer3")
  checkTrue(success)
}
