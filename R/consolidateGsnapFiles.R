##' Consolidates all pieces from parallelized gsnap output into
##' appropriately merged files
##'
##' If gsnap was run in single-end mode, this will consolidate all
##' parallelized into the 3 output files. If gsnap was run in
##' paired-end mode, consolidates to the 7 output files.
##' @title Consolidate gsnap's parallelized output files
##' @param sam_file_dir directory where gsnap's sam file output is stored
##' @param remove_merged remove the individual pieces after they are merged
##' @return list of the names of the files created from consolidation
##' @author Cory Barr
consolidateGsnapFiles <- function(sam_file_dir,
                                  remove_merged=FALSE,
                                  parallelized=TRUE) {

  if (! file.exists(sam_file_dir)) {
    stop(paste("Could not find the directory", sam_file_dir))
  } else {
    sam_file_dir <- file_path_as_absolute(sam_file_dir)
  }

  return_list <- list()
  
  dir_files <- dir(sam_file_dir, full.names=TRUE)

  .consolidate <- function(mapping_class,
                           remove_merged,
                           dir_files) {

    base_names <- basename(dir_files)

    unconsolidated_files <- dir_files[grep(paste("^gsnap_out\\.*",
                                                 "\\d+\\.",
                                                 mapping_class,
                                                 "$",
                                                 sep=''),
                                           base_names,
                                           perl=TRUE)]    

    consolidated_file <- consolidateSAMFiles(sam_files=unconsolidated_files,
                                             outfile=paste(
                                               "gsnap.merged",
                                               mapping_class,
                                               "bam",
                                               sep='.'),
                                             remove_merged=remove_merged)
  }

  mapping_classes <- basename(dir_files)
  mapping_classes <- mapping_classes[grep("gsnap_out\\.", mapping_classes)]
  mapping_classes <- unique(sub("^gsnap_out\\.\\d+\\.", "", mapping_classes))

  apply_func <- lapply
  if(parallelized)
    apply_func <- mclapply
  
  merged_files <- lapply(mapping_classes,
                         .consolidate,
                         remove_merged=remove_merged,
                         dir_files=dir_files)
  names(merged_files) <- mapping_classes
  return(merged_files)
}
