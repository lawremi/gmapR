##' Consolidate SAM files created by gsnap
##'
##' This function takes uses samtools to sort and merge SAM files into a single BAM file
##' @title Consolidate parallelized gsnap SAM files
##' @param sam_files character vector of the SAM files 
##' @param outfile the name of the output BAM file. Will be place in the
##' same directory as the first SAM file 
##' @param remove_merged logical indication whether to delete the
##' contents of sam_files after merging
##' @return character vector indicating the name of the resulting BAM file
##' @author Cory Barr
consolidateSAMFiles <- function(sam_files, outfile, remove_merged) {

  if (! all(file.exists(sam_files))) {
    stop("Could not find one or more of the input SAM files")
  }

  sam_files <- sapply(sam_files,
                      file_path_as_absolute,
                      USE.NAMES=FALSE)

  output_dir <- dirname(sam_files[1])
  if (missing(outfile)) {
    outfile <- paste(output_dir,
                     "merged.bam",
                     sep="/")
  } else {
    if(length(grep("^/", outfile)) != 1) {
      outfile <- paste(output_dir,
                       outfile,
                       sep="/")
    }
  }
  
  if(is.loaded("mc_fork", PACKAGE = "multicore")) {
    apply_func <- mclapply
  } else {
    apply_func <- lapply
  }  

  ##convert to bam first
  ##Example sys call:
  ##samtools view -S gsnap_out.0.unpaired_uniq -b > junk.bam    
  sam_files_converted <- paste(sam_files,
                               ".bam.converted",
                               sep='')

  convert_commands <- paste(globals()['samtools'],
                            "view -S",
                            sam_files,
                            "-b >",
                            sam_files_converted
                         )
  apply_func(convert_commands, system)
  consolidated_BAM_file <-
    consolidateBAMFiles(sam_files_converted,
                        outfile=outfile)
  unlink(sam_files_converted)
  if(remove_merged)
    unlink(sam_files)
  return(consolidated_BAM_file)
}
