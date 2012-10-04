##run this to update. The current working directory needs to be the
##top-level directory of the SVN checkout
updateGMAPRSrc <- function() {
  gmapSVNProj <- "http://resscm/bioinfo/projects/gmap/branches/internal-2011-12-28"
  gstructSVNProj <- "http://resscm/bioinfo/projects/gstruct/trunk"

  extractDirGmap <- file.path(getwd(), "src/gmap")
  extractDirGstruct <- file.path(getwd(), "src/gstruct")

  bootstrapAndExtract(projectSVNURL=gmapSVNProj, extractDir=extractDirGmap,
                      program="gmap")  
  bootstrapAndExtract(projectSVNURL=gstructSVNProj, extractDir=extractDirGstruct,
                      program="gstruct")
  .copySamflagsHeader(extractDirGstruct, extractDirGmap)
}
  
bootstrapAndExtract <- function(projectSVNURL, extractDir, program) {
  startingDir <- getwd()
  on.exit(setwd(startingDir))

  if (file.exists(extractDir)) {
    unlink(extractDir, recursive=TRUE)
  }
  dir.create(extractDir)
  
  svnCheckoutDir <- .getSVNProj(projectSVNURL, extractDir)
  on.exit(unlink(svnCheckoutDir, recursive=TRUE), add=TRUE)
  setwd(svnCheckoutDir)
  .bootstrapSVNCheckout()  
  .configureSrc(program)
  .makeDist()
  .extractDistTarballIntoSrcDirectory(extractDir)
  invisible(TRUE)
}

###################
###helper functions
###################

.copySamflagsHeader <- function(extractDirGstruct, extractDirGmap) {
  ##gstruct needs samflags.h. Copying from the gmap src
  gstructSamflagsLoc <- file.path(extractDirGstruct, "src/samflags.h")
  if (!file.exists(gstructSamflagsLoc)) {
    
    gmapSamflagsLoc <- file.path(extractDirGmap, "src/samflags.h")  
    if (!file.exists(gmapSamflagsLoc)) {
      stop("Could not find the samflags.h file in the gmap src code")
    }
    
    file.copy(gmapSamflagsLoc, gstructSamflagsLoc)
  }
}

.getSVNProj <- function(projectSVNURL, extractDir) {
  ##grab from svn
  tmpDir <- file.path(tempdir(), basename(extractDir))
  dir.create(tmpDir, recursive=TRUE)
  command <- paste("svn co",
                   projectSVNURL,
                   tmpDir)
  if (!system(command) == 0) {
    stop("Could not check out project from SVN")
  }

  return(tmpDir)
}

##assumes in the correct dir
.bootstrapSVNCheckout <- function() {
  if (!system("./bootstrap.Rdist") == 0) {
    stop("unable to bootstrap")
  }
}

.configureSrc <- function(program) {
  ##--with-gmapdb=${GMAPDB} --prefix=${PREFIX}
  ##configure. Set
  ##run a "make dist" to build a tarball
  command <- "./configure --disable-fulldist"
  ##if (program == "gstruct") command <- paste(command, "--disable-binaries")
  if (!system(command) == 0) {
    stop("unable to configure")
  }

  ##gstruct SVN repo doesn't have config.site (this code needs to work
  ##today)
  if (program == "gstruct") {
    if (!file.exists("config.site")) {
      system("touch config.site")
    }
  }
}

.makeDist <- function() {
  ##run a "make dist" to build a tarball
  if (!system("make distcheck") == 0) {
    stop("unable to 'make dist'")
  }  
}

.extractDistTarballIntoSrcDirectory <- function(extractDir) {
  ##extract tarball into 'src' directory of package
  distTarball <- dir(pattern="*\\.tar\\.gz$", full.names=TRUE)
  if (length(distTarball) > 1) {
    stop("Found more than one tarball in directory. ",
         "Not sure which is the distribution. Aborting.")
  }
  distTarball <- tools::file_path_as_absolute(distTarball)
  if (file.exists(extractDir)) {
    unlink(extractDir, recursive=TRUE)
  }
  dir.create(extractDir, recursive=TRUE)

  distName <- sub("\\..*", "", basename(distTarball))
  command <- paste("tar xvf", distTarball,
                   "-C",
                   extractDir,
                   distName)
  if (!system(command) == 0) {
    stop("Could not extract source from distribution into specified directory.")
  }

  ##move contents of untarred dir up one level
  ##TODO: How (or can) you do this w/ a tar command?
  whereTarExtracted <- file.path(extractDir, distName)
  command <- paste("mv",
                   paste0(whereTarExtracted, "/*"),
                   extractDir)
  if (!system(command) == 0) {
    stop("Could not move extracted directory down one level")
  }
  unlink(whereTarExtracted, recursive=TRUE)
}
