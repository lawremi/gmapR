##' Return number of reads in each of gsnap's output files
##'
##' With the --split-output argument, gsnap can output multiple
##' files. If these files are BAM files, this function counts the
##' number of reads in each file.
##' @title Return number of reads in each of gsnap's output files
##' @param aligner_dir Directory containing output bam files from gsnap
##' @param parallelize Boolean indicating if counting should be done in parallel
##' @return named list of read counts per gsnap output file
##' @author Cory Barr
##' @export
getGsnapFileCounts <- function(aligner_dir, parallelize=TRUE) {

  if(!file.exists(aligner_dir))
    stop("aligner_dir does not exist")
  
  ##have to count gsnap output by chr, or can get overflow errors
  .getNumUniqueReadsInBam <- function(bam_file) {
    chr_granges <- getSeqlensFromBAM(bam_file)
      .getUniqueReadsByChr <- function(index) {
        sb_param <-
          ScanBamParam(what=c("qname"),
                       which=chr_granges[index])
        res <- scanBam(bam_file, param=sb_param)[[1]]
        unique(res$qname)
      }

    apply_func <- lapply
    if (parallelize) apply_func <- mclapply
    
    res <- unlist(apply_func(seq_len(length(chr_granges)),
                           .getUniqueReadsByChr))
    length(unique(res))
  }

  bam_files <- dir(aligner_dir, pattern="\\.bam$", full.names=TRUE)

  nomapping_index <- grep("\\.no_?mapping\\.bam$", bam_files)
  mult_mapping_index <- grep("_mult[_\\.]", bam_files)
  uniq_mapping_index <- seq_len(length(bam_files))[-c(mult_mapping_index, nomapping_index)]
  nomapping_bam <- bam_files[nomapping_index]
  mult_mapping_bams <- bam_files[mult_mapping_index]
  uniq_mapping_bams <- bam_files[uniq_mapping_index]

  named_names <- bam_files[-nomapping_index]
  names(named_names) <- named_names
  gsnap_counts <- lapply(named_names,
                         .getNumUniqueReadsInBam)
                         
  names(gsnap_counts) <- sub("^.*\\.gsnap\\.merged\\.", "", names(gsnap_counts))
  names(gsnap_counts) <- sub("\\.bam$", "", names(gsnap_counts))

  ##handling the 'nomapping' file differently, since nothing this (so
  ##nothing is returned if scanBam is given a ScanBamParam which arg
  nomapper_ids <- unlist(scanBam(nomapping_bam,
                          param=ScanBamParam(what=c("qname")))[[1]],
                         use.names=FALSE)
  num_nomappers <- sum(!duplicated(nomapper_ids))
  gsnap_counts[['nomapping']] <- num_nomappers

  names(gsnap_counts) <- paste("gsnap", names(gsnap_counts), sep=".")
  names(gsnap_counts) <- paste(names(gsnap_counts), "_reads", sep="")

  return(gsnap_counts)  
}
