##' Consolidates all pieces from parallelized gsnap output into
##' appropriately merged files
##'
##' If gmap was run in single-end mode, this will consolidate all
##' parallelized into the 3 output files. If gmap was run in
##' paired-end mode, consolidates to the 7 output files.
##' @title Consolidate gmap's output files
##' @param sam_file_dir directory where gmap's sam file output is stored
##' @param paired_end indicated whether gmap was run in paired_end mode
##' @param remove_merged remove the individual pieces after they are merged
##' @return list of the names of the files created from consolidation
##' @author Cory Barr
##' @export
consolidateGmapFiles <- function(sam_file_dir,
                                  paired_end=FALSE,
                                  remove_merged=FALSE) {

  if (! file.exists(sam_file_dir)) {
    stop(paste("Could not find the directory", sam_file_dir))
  } else {
    sam_file_dir <- file_path_as_absolute(sam_file_dir)
  }

  return_list <- list()
  
  dir_files <- dir(sam_file_dir)

  no_mapping_files <- dir_files[grep("gmap_out.*nomapping$",
                                     dir_files,
                                     perl=T)]
  no_mapping_files <- paste(sam_file_dir, no_mapping_files, sep="/")
  consolidated_file <- consolidateSAMFiles(no_mapping_files,
                                           "gmap.merged.no_mapping.bam",
                                           remove_merged=remove_merged)
  return_list$no_mapping <- consolidated_file

  unpaired_uniq_files <- dir_files[grep("gmap_out.*uniq$",
                                        dir_files,
                                        perl=T)]
  unpaired_uniq_files <- paste(sam_file_dir, unpaired_uniq_files, sep="/")
  consolidated_file <- consolidateSAMFiles(unpaired_uniq_files,
                                           "gmap.merged.uniq.bam",
                                           remove_merged=remove_merged)
  return_list$unpaired_uniq <- consolidated_file

  unpaired_mult_files <- dir_files[grep("gmap_out.*mult$",
                                        dir_files,
                                        perl=T)]    
  unpaired_mult_files <- paste(sam_file_dir, unpaired_mult_files, sep="/")
  consolidated_file <- consolidateSAMFiles(unpaired_mult_files,
                                           "gmap.merged.mult.bam",
                                           remove_merged=remove_merged)
  return_list$unpaired_mult <- consolidated_file

  
  
  ##paired-end output adds an additional 4 files to those above
  if (paired_end) {
    ## not functional for 454 Gmap yet
  }
  
  return(return_list)
}
