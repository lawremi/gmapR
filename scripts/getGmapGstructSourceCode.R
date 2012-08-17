bootstrapAndExtract <- function(projectSVNURL, extractDir, program) {

  startingDir <- getwd()
  on.exit(setwd(startingDir))
  
  ##grab from svn
  tmpDir <- file.path(tempdir(), basename(extractDir))
  dir.create(tmpDir, recursive=TRUE)
  on.exit(unlink(tmpDir, recursive=TRUE), add=TRUE)
  
  command <- paste("svn co",
                   projectSVNURL,
                   tmpDir)
  if (!system(command) == 0) {
    stop("Could not check out project from SVN")
  }
    
  ##bootstrap it
  setwd(tmpDir)
  if (!system("./bootstrap.Rdist") == 0) {
    stop("unable to bootstrap")
  }

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
  
  ##run a "make dist" to build a tarball
  if (!system("make distcheck") == 0) {
    stop("unable to 'make dist'")
  }  
  
  ##extract tarball into 'src' directory of package
  distTarball <- dir(pattern="*\\.tar\\.gz$")
  if (length(distTarball) > 1) {
    stop("Found more than one tarball in directory. ",
         "Not sure which is the distribution. Aborting.")
  }
  if (!file.exists(extractDir)) {
    dir.create(extractDir, recursive=TRUE)
  }
  distName <- sub("\\..*", "", distTarball)
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

  invisible(TRUE)
}

options(error=recover)

gmapSVNProj <- "http://resscm/bioinfo/projects/gmap/branches/internal-2011-12-28"
extractDirGmap <- file.path(getwd(), "src/gmap")
bootstrapAndExtract(projectSVNURL=gmapSVNProj, extractDir=extractDirGmap,
                    program="gmap")

gstructSVNProj <- "http://resscm/bioinfo/projects/gstruct/trunk"
extractDirGstruct <- file.path(getwd(), "src/gstruct")
bootstrapAndExtract(projectSVNURL=gstructSVNProj, extractDir=extractDirGstruct,
                    program="gstruct")

##gstruct needs samflags.h. Copying from the gmap src
gstructSamflagsLoc <- file.path(extractDirGstruct, "src/samflags.h")
if (!file.exists(gstructSamflagsLoc)) {
 
  gmapSamflagsLoc <- file.path(extractDirGmap, "src/samflags.h")  
  if (!file.exists(gmapSamflagsLoc)) {
    stop("Could not find the samflags.h file in the gmap src code")
  }

  file.copy(gmapSamflagsLoc, gstructSamflagsLoc)
}
    
