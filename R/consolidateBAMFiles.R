##' Helper utility to sort and merge a set of BAM files
##'
##' <details>
##' @title Sort and merge a set of BAM files
##' @param bam_files A character vector of the BAM files to be sorted
##' @param outfile The name of the resulting sorted and merged BAM
##' file (optional)
##' @return a string indicating the name of the outfile 
##' @author Cory Barr - Regular
consolidateBAMFiles <- function(bam_files, outfile) {
  if (! all(file.exists(bam_files))) {
    stop("Could not find one or more of the input BAM files")
  }

  bam_files <- sapply(bam_files,
                      file_path_as_absolute,
                      USE.NAMES=F)
  
  output_dir <- sub("/[^/]*$", "", bam_files[1])

  if(is.loaded("mc_fork", PACKAGE = "multicore")) {
    apply_func <- mclapply
  } else {
    apply_func <- lapply
  }  

  ##sort
  bam_files_sorted <- paste(bam_files,
                            ".samtools_sorted",
                            sep="")
  commands <- paste(globals()['samtools'],
                    "sort",
                    bam_files,
                    bam_files_sorted
                    )
  apply_func(commands, system)
  bam_files_sorted <- paste(bam_files_sorted,
                            ".bam",
                            sep="")
  
  ##merge
  if (missing(outfile)) {
    outfile <- paste(output_dir,
                     "merged.bam",
                     sep="/")
  } else {
    if(length(grep("^/", outfile, perl=T)) != 1) {
      outfile <- paste(output_dir,
                       outfile,
                       sep="/")
    }
  }

  if(length(bam_files_sorted) > 1) {
    command <- paste(globals()['samtools'],
                     "merge",
                     outfile,
                     paste(bam_files_sorted, collapse=" "))
    system(command)
  } else {
    file.copy(bam_files_sorted, outfile)
  }
  
  unlink(bam_files_sorted)
  return(outfile)
}
